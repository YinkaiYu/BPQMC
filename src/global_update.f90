module GlobalUpdate_mod
    use CalcBasic
    use DQMC_Model_mod
    use Dynamics_mod
    use MakeInitialState, only: Initial
    use Multiply_mod
    use ObserEqual_mod
    use OperatorHubbard_mod, only: AccCounter, OperatorHubbard
    use Stabilize_mod
    implicit none

    public
    private :: HMC_force_prop_L, HMC_set_overlap, HMC_log_overlap_abs2
    private :: HMC_draw_nfrog
    private :: HMC_rng_gaussian, HMC_measure_sweep_L, HMC_measure_sweep_R

    type :: GlobalUpdate
        type(Propagator), allocatable, private :: prop_work
        type(WrapList), allocatable, private :: wr_work
        real(kind=8), dimension(:,:,:), allocatable, private :: force_cur
        real(kind=8), dimension(:,:,:), allocatable, private :: momentum
        real(kind=8), dimension(:,:,:), allocatable, private :: phi_backup
        real(kind=8), private :: action_cur
        logical, private :: state_ready
    contains
        procedure :: init => Global_init
        procedure :: clear => Global_clear
        procedure :: debug_eval => Global_debug_eval
        procedure, private :: reset => Global_reset
        procedure, private :: prepare_work => Global_prepare_work
        procedure, private :: ensure_state => Global_ensure_state
        procedure, private :: eval_state => Global_eval_state
        procedure, private :: sample_momentum => Global_sample_momentum
        procedure, private :: leapfrog => Global_leapfrog
        procedure, private :: step => Global_step
        procedure, private :: measure_config => Global_measure_config
        procedure :: therm => Global_therm
        procedure :: sweep => Global_sweep
        procedure, nopass :: ctrl_print => Global_control_print
        procedure, nopass :: ctrl_print_t => Global_therm_print
    end type GlobalUpdate

    type(ObserTau), allocatable :: Obs_tau_hmc
    type(ObserEqual), allocatable :: Obs_equal_hmc
    type(Dynamics) :: Dyn
    type(AccCounter) :: Acc_HMC, Acc_HMC_warm
    real(kind=8), parameter :: overlap_floor = 1.d-300

contains
    subroutine Global_init(this, Init_obj)
        class(GlobalUpdate), intent(inout) :: this
        class(Initial), intent(in) :: Init_obj

        allocate(this%prop_work)
        call this%prop_work%make(Init_obj)
        allocate(this%wr_work)
        call this%wr_work%make()
        allocate(this%force_cur(Naux, Ndim, Ltrot))
        allocate(this%momentum(Naux, Ndim, Ltrot))
        allocate(this%phi_backup(Naux, Ndim, Ltrot))

        this%force_cur = 0.d0
        this%momentum = 0.d0
        this%phi_backup = 0.d0
        this%action_cur = 0.d0
        this%state_ready = .false.

        allocate(Obs_equal_hmc)
        call Obs_equal_hmc%make()
        if (is_tau) then
            allocate(Obs_tau_hmc)
            call Obs_tau_hmc%make()
            call Dyn%init(Init_obj)
        endif

        call Acc_HMC%init()
        call Acc_HMC_warm%init()
        return
    end subroutine Global_init

    subroutine Global_clear(this)
        class(GlobalUpdate), intent(inout) :: this

        if (allocated(Obs_equal_hmc)) deallocate(Obs_equal_hmc)
        if (is_tau) then
            if (allocated(Obs_tau_hmc)) deallocate(Obs_tau_hmc)
            call Dyn%clear()
        endif
        deallocate(this%prop_work)
        deallocate(this%wr_work)
        deallocate(this%force_cur, this%momentum, this%phi_backup)
        return
    end subroutine Global_clear

    subroutine Global_debug_eval(this, action, force, ok)
        class(GlobalUpdate), intent(inout) :: this
        real(kind=8), intent(out) :: action
        real(kind=8), dimension(Naux, Ndim, Ltrot), intent(out) :: force
        logical, intent(out) :: ok

        call this%eval_state(action, force)
        ok = .true.
        return
    end subroutine Global_debug_eval

    subroutine Global_reset(this, toggle)
        class(GlobalUpdate), intent(inout) :: this
        logical, intent(in) :: toggle

        call Acc_HMC%reset()
        call Obs_equal_hmc%reset()
        if (toggle) call Obs_tau_hmc%reset()
        return
    end subroutine Global_reset

    subroutine Global_prepare_work(this)
        class(GlobalUpdate), intent(inout) :: this

        call this%prop_work%reset(Init)
        call this%wr_work%reset()
        return
    end subroutine Global_prepare_work

    subroutine Global_ensure_state(this)
        class(GlobalUpdate), intent(inout) :: this

        if (this%state_ready) return
        call this%eval_state(this%action_cur, this%force_cur)
        this%state_ready = .true.
        return
    end subroutine Global_ensure_state

    subroutine Global_eval_state(this, action, force)
        class(GlobalUpdate), intent(inout) :: this
        real(kind=8), intent(out) :: action
        real(kind=8), dimension(Naux, Ndim, Ltrot), intent(out) :: force
        integer :: nt
        real(kind=8) :: log_overlap_abs2

        force = 0.d0
        action = 0.5d0 * sum(Conf%phi_list * Conf%phi_list)

        call this%prepare_work()
        if (Ltrot == 0) then
            call HMC_log_overlap_abs2(this%prop_work, log_overlap_abs2)
            action = action - dble(Nbos) * log_overlap_abs2
            return
        endif

        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call propU_pre(Op_U1, this%prop_work, 1, nt)
            if (abs(RU2) > Zero) call propU_pre(Op_U2, this%prop_work, 2, nt)
            call propT_pre(this%prop_work)
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_pre(this%prop_work, this%wr_work, nt)
        enddo

        call HMC_log_overlap_abs2(this%prop_work, log_overlap_abs2)
        action = action - dble(Nbos) * log_overlap_abs2

        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_L(this%prop_work, this%wr_work, nt, "H")
            call propT_L(this%prop_work)
            if (abs(RU2) > Zero) call HMC_force_prop_L(Op_U2, this%prop_work, 2, nt, force)
            if (abs(RU1) > Zero) call HMC_force_prop_L(Op_U1, this%prop_work, 1, nt, force)
        enddo
        return
    end subroutine Global_eval_state

    subroutine Global_sample_momentum(this, iseed)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer :: nt, ii, nf

        do nt = 1, Ltrot
            do ii = 1, Ndim
                do nf = 1, Naux
                    this%momentum(nf, ii, nt) = HMC_rng_gaussian(iseed)
                enddo
            enddo
        enddo
        return
    end subroutine Global_sample_momentum

    subroutine Global_leapfrog(this, action_new, nsteps)
        class(GlobalUpdate), intent(inout) :: this
        real(kind=8), intent(out) :: action_new
        integer, intent(in) :: nsteps
        integer :: nlf

        this%momentum = this%momentum + 0.5d0 * hmc_dt * this%force_cur
        do nlf = 1, nsteps
            Conf%phi_list = Conf%phi_list + hmc_dt * this%momentum
            call this%eval_state(action_new, this%force_cur)
            if (nlf == nsteps) then
                this%momentum = this%momentum + 0.5d0 * hmc_dt * this%force_cur
            else
                this%momentum = this%momentum + hmc_dt * this%force_cur
            endif
        enddo
        return
    end subroutine Global_leapfrog

    subroutine Global_step(this, iseed, Counter)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        class(AccCounter), intent(inout) :: Counter
        real(kind=8) :: action_new, ham_old, ham_new, delta_h, ratio, random
        integer :: nsteps
        logical :: accepted
        real(kind=8), external :: ranf

        call this%ensure_state()
        this%phi_backup = Conf%phi_list
        call this%sample_momentum(iseed)
        ham_old = 0.5d0 * sum(this%momentum * this%momentum) + this%action_cur

        nsteps = HMC_draw_nfrog(iseed)
        call this%leapfrog(action_new, nsteps)
        ham_new = 0.5d0 * sum(this%momentum * this%momentum) + action_new

        delta_h = ham_old - ham_new
        if (delta_h >= 0.d0) then
            ratio = 1.d0
        else
            ratio = exp(delta_h)
        endif

        random = ranf(iseed)
        accepted = (ratio > random)
        call Counter%count(accepted)
        if (accepted) then
            this%action_cur = action_new
            this%state_ready = .true.
        else
            Conf%phi_list = this%phi_backup
            call this%eval_state(this%action_cur, this%force_cur)
            this%state_ready = .true.
        endif
        return
    end subroutine Global_step

    integer function HMC_draw_nfrog(iseed) result(nsteps)
        integer, intent(inout) :: iseed
        integer :: lower, upper, span

        if (NfrogJitter <= 0) then
            nsteps = Nfrog
            return
        endif

        lower = max(1, Nfrog - NfrogJitter)
        upper = Nfrog + NfrogJitter
        span = upper - lower + 1
        nsteps = lower - 1 + nranf(iseed, span)
        return
    end function HMC_draw_nfrog

    subroutine Global_measure_config(this, toggle, Nobs, Nobst)
        class(GlobalUpdate), intent(inout) :: this
        logical, intent(in) :: toggle
        integer, intent(inout) :: Nobs, Nobst
        integer :: nt

        call this%prepare_work()
        if (Ltrot == 0) then
            call Obs_equal_hmc%calc(this%prop_work, 0)
            Nobs = Nobs + 1
            return
        endif

        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call propU_pre(Op_U1, this%prop_work, 1, nt)
            if (abs(RU2) > Zero) call propU_pre(Op_U2, this%prop_work, 2, nt)
            call propT_pre(this%prop_work)
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_pre(this%prop_work, this%wr_work, nt)
        enddo

        call HMC_measure_sweep_L(this%prop_work, this%wr_work, Nobs)
        if (toggle) then
            call HMC_set_overlap(this%prop_work)
            call Dyn%reset(this%prop_work)
            call Dyn%sweep_R(Obs_tau_hmc, this%wr_work)
            Nobst = Nobst + 1
        endif
        call HMC_measure_sweep_R(this%prop_work, this%wr_work, Nobs)
        return
    end subroutine Global_measure_config

    subroutine Global_therm(this, iseed)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed

        call Acc_HMC_warm%reset()
        call this%step(iseed, Acc_HMC_warm)
        call Acc_HMC_warm%ratio()
        return
    end subroutine Global_therm

    subroutine Global_sweep(this, Prop, WrList, iseed, is_beta, toggle)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: iseed
        logical, intent(inout) :: is_beta
        logical, intent(in) :: toggle
        integer :: Nobs, Nobst, nsw

        call this%reset(toggle)
        Nobs = 0
        Nobst = 0
        do nsw = 1, Nsweep
            call this%step(iseed, Acc_HMC)
            call this%measure_config(toggle, Nobs, Nobst)
        enddo

        call Obs_equal_hmc%ave(Nobs)
        if (toggle) call Obs_tau_hmc%ave(Nobst)
        call Acc_HMC%ratio()
        return
    end subroutine Global_sweep

    subroutine Global_control_print(toggle)
        include 'mpif.h'
        logical, intent(in) :: toggle
        real(kind=8) :: collect

        collect = 0.d0
        call MPI_Reduce(Acc_HMC%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            Acc_HMC%acc = collect / dble(ISIZE * Nbin)
            write(50,*) 'Accept_HMC                                     :', Acc_HMC%acc
        endif
        if (toggle) call Dyn%ctrl_print()
        return
    end subroutine Global_control_print

    subroutine Global_therm_print()
        include 'mpif.h'
        real(kind=8) :: collect

        collect = 0.d0
        call MPI_Reduce(Acc_HMC_warm%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            Acc_HMC_warm%acc = collect / dble(ISIZE * Nwarm)
            write(50,*) 'Thermalize HMC Accept Ratio                    :', Acc_HMC_warm%acc
        endif
        return
    end subroutine Global_therm_print

    subroutine HMC_force_prop_L(Op_U, Prop, nf, ntau, force)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(in) :: nf, ntau
        real(kind=8), dimension(Naux, Ndim, Ltrot), intent(inout) :: force
        complex(kind=8) :: alpha, gdiag, scale
        integer :: ii

        call HMC_set_overlap(Prop)
        alpha = Op_U%get_alpha()
        scale = dcmplx(dble(Nbos), 0.d0) / HMC_safe_overlap(Prop%overlap)
        do ii = Ndim, 1, -1
            gdiag = scale * Prop%UUR(ii, 1) * Prop%UUL(1, ii)
            force(nf, ii, ntau) = -Conf%phi_list(nf, ii, ntau) + 2.d0 * real(alpha * gdiag)
        enddo

        do ii = Ndim, 1, -1
            call Op_U%mmult_L(Prop%UUL, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
        enddo
        return
    end subroutine HMC_force_prop_L

    subroutine HMC_measure_sweep_L(Prop, WrList, Nobs)
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: Nobs
        integer :: nt

        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_L(Prop, WrList, nt, "M")
            if (nt == Ltrot/2) then
                call HMC_set_overlap(Prop)
                call Obs_equal_hmc%calc(Prop, nt)
                Nobs = Nobs + 1
            endif
            call propT_L(Prop)
            if (abs(RU2) > Zero) call propU_L(Op_U2, Prop, 2, nt)
            if (abs(RU1) > Zero) call propU_L(Op_U1, Prop, 1, nt)
        enddo
        return
    end subroutine HMC_measure_sweep_L

    subroutine HMC_measure_sweep_R(Prop, WrList, Nobs)
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: Nobs
        integer :: nt

        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call propU_R(Op_U1, Prop, 1, nt)
            if (abs(RU2) > Zero) call propU_R(Op_U2, Prop, 2, nt)
            call propT_R(Prop)
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_R(Prop, WrList, nt, "M")
            if (nt == Ltrot/2) then
                call HMC_set_overlap(Prop)
                call Obs_equal_hmc%calc(Prop, nt)
                Nobs = Nobs + 1
            endif
        enddo
        return
    end subroutine HMC_measure_sweep_R

    subroutine HMC_set_overlap(Prop)
        class(Propagator), intent(inout) :: Prop
        complex(kind=8), external :: ZDOTU

        Prop%overlap = ZDOTU(Ndim, Prop%UUL(1,1), 1, Prop%UUR(1,1), 1)
        return
    end subroutine HMC_set_overlap

    subroutine HMC_log_overlap_abs2(Prop, log_overlap_abs2)
        class(Propagator), intent(inout) :: Prop
        real(kind=8), intent(out) :: log_overlap_abs2
        real(kind=8) :: overlap_abs

        call HMC_set_overlap(Prop)
        overlap_abs = max(abs(Prop%overlap), overlap_floor)
        log_overlap_abs2 = 2.d0 * (Prop%log_norm_ur + Prop%log_norm_ul + log(overlap_abs))
        return
    end subroutine HMC_log_overlap_abs2

    pure complex(kind=8) function HMC_safe_overlap(overlap) result(safe_overlap)
        complex(kind=8), intent(in) :: overlap
        real(kind=8) :: overlap_abs

        overlap_abs = abs(overlap)
        if (overlap_abs >= overlap_floor) then
            safe_overlap = overlap
        elseif (overlap_abs > 0.d0) then
            safe_overlap = overlap * (overlap_floor / overlap_abs)
        else
            safe_overlap = dcmplx(overlap_floor, 0.d0)
        endif
        return
    end function HMC_safe_overlap

    real(kind=8) function HMC_rng_gaussian(iseed) result(X)
        integer, intent(inout) :: iseed
        real(kind=8) :: X1, X2
        real(kind=8), external :: ranf

        X1 = max(ranf(iseed), 1.d-12)
        X2 = ranf(iseed)
        X = sqrt(-2.d0 * log(X1)) * cos(2.d0 * PI * X2)
        return
    end function HMC_rng_gaussian
end module GlobalUpdate_mod
