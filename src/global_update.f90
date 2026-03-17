module GlobalUpdate_mod
    use CalcBasic
    use DQMC_Model_mod
    use Dynamics_mod
    use MakeInitialState, only: Initial
    use Multiply_mod
    use ObserEqual_mod
    use OperatorHubbard_mod, only: AccCounter, OperatorHubbard
    use Stabilize_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    public
    private :: HMC_force_prop_L, HMC_set_overlap, HMC_log_overlap_abs2
    private :: HMC_draw_nfrog
    private :: HMC_rng_gaussian, HMC_measure_sweep_L, HMC_measure_sweep_R
    private :: HMC_monitor_init, HMC_monitor_close, HMC_monitor_step
    private :: HMC_trace_init, HMC_trace_close, HMC_trace_begin, HMC_trace_step
    private :: HMC_read_env_int, HMC_trace_values

    type :: GlobalUpdate
        type(Propagator), allocatable, private :: prop_work
        type(WrapList), allocatable, private :: wr_work
        real(kind=8), dimension(:,:,:), allocatable, private :: force_cur
        real(kind=8), dimension(:,:,:), allocatable, private :: momentum
        real(kind=8), dimension(:,:,:), allocatable, private :: phi_backup
        real(kind=8), private :: action_cur
        logical, private :: state_ready
        integer, private :: ii_begin, ii_end
        integer, private :: tau_begin, tau_end
    contains
        procedure, public :: init => Global_init
        procedure, public :: clear => Global_clear
        procedure, public :: debug_eval => Global_debug_eval
        procedure, public :: debug_roundtrip => Global_debug_roundtrip
        procedure, public :: debug_measure_current => Global_debug_measure_current
        procedure, public :: reset_diag => Global_reset_diag
        procedure, private :: reset => Global_reset
        procedure, private :: prepare_work => Global_prepare_work
        procedure, private :: ensure_state => Global_ensure_state
        procedure, private :: eval_state => Global_eval_state
        procedure, private :: select_site_block => Global_select_site_block
        procedure, private :: select_tau_block => Global_select_tau_block
        procedure, private :: sample_momentum => Global_sample_momentum
        procedure, private :: leapfrog => Global_leapfrog
        procedure, private :: step => Global_step
        procedure, private :: measure_config => Global_measure_config
        procedure, private :: check_equal_obs => Global_check_equal_obs
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
    real(kind=8), save :: HMC_deltaH_sum = 0.d0
    real(kind=8), save :: HMC_deltaH_absmax = 0.d0
    integer, save :: HMC_deltaH_count = 0
    logical, save :: HMC_monitor_enabled = .false.
    integer, parameter :: HMC_monitor_unit = 86
    integer, save :: HMC_monitor_count = 0
    integer, save :: HMC_monitor_site = 1
    integer, save :: HMC_monitor_tau = 1
    logical, save :: HMC_trace_enabled = .false.
    integer, parameter :: HMC_trace_unit = 87
    integer, save :: HMC_trace_count = 0
    integer, save :: HMC_trace_proposal = 0
    integer, save :: HMC_trace_stage = 0
    integer, save :: HMC_trace_nsteps = 0
    integer, save :: HMC_debug_measure = 0

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
        this%ii_begin = 1
        this%ii_end = Ndim
        this%tau_begin = 1
        this%tau_end = Ltrot

        allocate(Obs_equal_hmc)
        call Obs_equal_hmc%make()
        if (is_tau) then
            allocate(Obs_tau_hmc)
            call Obs_tau_hmc%make()
            call Dyn%init(Init_obj)
        endif

        call Acc_HMC%init()
        call Acc_HMC_warm%init()
        call this%reset_diag()
        call HMC_read_env_int("BPQMC_HMC_DEBUG_MEASURE", HMC_debug_measure)
        call HMC_monitor_init()
        call HMC_trace_init()
        return
    end subroutine Global_init

    subroutine Global_clear(this)
        class(GlobalUpdate), intent(inout) :: this

        if (allocated(Obs_equal_hmc)) deallocate(Obs_equal_hmc)
        if (is_tau) then
            if (allocated(Obs_tau_hmc)) deallocate(Obs_tau_hmc)
            call Dyn%clear()
        endif
        call HMC_trace_close()
        call HMC_monitor_close()
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

    subroutine Global_debug_roundtrip(this, iseed, nsteps, diff_phi, diff_momentum, diff_action)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer, intent(in) :: nsteps
        real(kind=8), intent(out) :: diff_phi, diff_momentum, diff_action
        real(kind=8), allocatable :: phi_start(:,:,:), momentum_start(:,:,:)
        real(kind=8) :: action_trial, action_start

        call this%ensure_state()
        allocate(phi_start(Naux, Ndim, Ltrot))
        allocate(momentum_start(Naux, Ndim, Ltrot))

        phi_start = Conf%phi_list
        action_start = this%action_cur
        call this%select_site_block(iseed)
        call this%select_tau_block(iseed)
        call this%sample_momentum(iseed)
        momentum_start = this%momentum

        call this%leapfrog(action_trial, nsteps, -1)
        this%action_cur = action_trial
        this%state_ready = .true.
        this%momentum = -this%momentum

        call this%leapfrog(action_trial, nsteps, -2)
        diff_phi = maxval(abs(Conf%phi_list - phi_start))
        diff_momentum = maxval(abs(this%momentum + momentum_start))
        diff_action = abs(action_trial - action_start)

        Conf%phi_list = phi_start
        call this%eval_state(this%action_cur, this%force_cur)
        this%state_ready = .true.

        deallocate(phi_start)
        deallocate(momentum_start)
        return
    end subroutine Global_debug_roundtrip

    subroutine Global_debug_measure_current(this, toggle)
        class(GlobalUpdate), intent(inout) :: this
        logical, intent(in) :: toggle
        integer :: nobs, nobst

        call Obs_equal_hmc%reset()
        if (toggle .and. is_tau .and. allocated(Obs_tau_hmc)) call Obs_tau_hmc%reset()
        call this%measure_config(toggle .and. is_tau .and. allocated(Obs_tau_hmc), nobs, nobst)
        call Obs_equal_hmc%ave(nobs)
        call this%check_equal_obs('debug_measure_current', nobs, nobst, toggle)
        if (toggle .and. is_tau .and. allocated(Obs_tau_hmc)) call Obs_tau_hmc%ave(nobst)
        return
    end subroutine Global_debug_measure_current

    subroutine Global_reset(this, toggle)
        class(GlobalUpdate), intent(inout) :: this
        logical, intent(in) :: toggle

        call Acc_HMC%reset()
        call Obs_equal_hmc%reset()
        if (toggle .and. is_tau .and. allocated(Obs_tau_hmc)) call Obs_tau_hmc%reset()
        return
    end subroutine Global_reset

    subroutine Global_reset_diag(this)
        class(GlobalUpdate), intent(inout) :: this

        HMC_deltaH_sum = 0.d0
        HMC_deltaH_absmax = 0.d0
        HMC_deltaH_count = 0
        return
    end subroutine Global_reset_diag

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

    subroutine Global_select_site_block(this, iseed)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer :: span

        if (Ndim <= 0) then
            this%ii_begin = 1
            this%ii_end = 0
            return
        endif

        if (hmc_block_sites <= 0 .or. hmc_block_sites >= Ndim) then
            this%ii_begin = 1
            this%ii_end = Ndim
            return
        endif

        span = Ndim - hmc_block_sites + 1
        this%ii_begin = nranf(iseed, span)
        this%ii_end = this%ii_begin + hmc_block_sites - 1
        return
    end subroutine Global_select_site_block

    subroutine Global_select_tau_block(this, iseed)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer :: span

        if (Ltrot <= 0) then
            this%tau_begin = 1
            this%tau_end = 0
            return
        endif

        if (hmc_block_tau <= 0 .or. hmc_block_tau >= Ltrot) then
            this%tau_begin = 1
            this%tau_end = Ltrot
            return
        endif

        span = Ltrot - hmc_block_tau + 1
        this%tau_begin = nranf(iseed, span)
        this%tau_end = this%tau_begin + hmc_block_tau - 1
        return
    end subroutine Global_select_tau_block

    subroutine Global_sample_momentum(this, iseed)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer :: nt, ii, nf
        real(kind=8) :: sigma_p

        sigma_p = sqrt(hmc_mass)
        this%momentum = 0.d0

        do nt = this%tau_begin, this%tau_end
            do ii = this%ii_begin, this%ii_end
                do nf = 1, Naux
                    this%momentum(nf, ii, nt) = sigma_p * HMC_rng_gaussian(iseed)
                enddo
            enddo
        enddo
        return
    end subroutine Global_sample_momentum

    subroutine Global_leapfrog(this, action_new, nsteps, stage_id)
        class(GlobalUpdate), intent(inout) :: this
        real(kind=8), intent(out) :: action_new
        integer, intent(in) :: nsteps, stage_id
        integer :: nlf, nt, ii, nf

        if (this%ii_begin <= this%ii_end .and. this%tau_begin <= this%tau_end) then
            do nt = this%tau_begin, this%tau_end
                do ii = this%ii_begin, this%ii_end
                    do nf = 1, Naux
                        this%momentum(nf, ii, nt) = this%momentum(nf, ii, nt) + 0.5d0 * hmc_dt * this%force_cur(nf, ii, nt)
                    enddo
                enddo
            enddo
        endif
        call HMC_trace_begin(stage_id, nsteps, this%action_cur, this%momentum, this%force_cur)
        do nlf = 1, nsteps
            if (this%ii_begin <= this%ii_end .and. this%tau_begin <= this%tau_end) then
                do nt = this%tau_begin, this%tau_end
                    do ii = this%ii_begin, this%ii_end
                        do nf = 1, Naux
                            Conf%phi_list(nf, ii, nt) = Conf%phi_list(nf, ii, nt) + (hmc_dt / hmc_mass) * this%momentum(nf, ii, nt)
                        enddo
                    enddo
                enddo
            endif
            call this%eval_state(action_new, this%force_cur)
            if (nlf == nsteps) then
                if (this%ii_begin <= this%ii_end .and. this%tau_begin <= this%tau_end) then
                    do nt = this%tau_begin, this%tau_end
                        do ii = this%ii_begin, this%ii_end
                            do nf = 1, Naux
                                this%momentum(nf, ii, nt) = this%momentum(nf, ii, nt) + 0.5d0 * hmc_dt * this%force_cur(nf, ii, nt)
                            enddo
                        enddo
                    enddo
                endif
            else
                if (this%ii_begin <= this%ii_end .and. this%tau_begin <= this%tau_end) then
                    do nt = this%tau_begin, this%tau_end
                        do ii = this%ii_begin, this%ii_end
                            do nf = 1, Naux
                                this%momentum(nf, ii, nt) = this%momentum(nf, ii, nt) + hmc_dt * this%force_cur(nf, ii, nt)
                            enddo
                        enddo
                    enddo
                endif
            endif
            call HMC_trace_step(nlf, action_new, this%momentum, this%force_cur)
        enddo
        return
    end subroutine Global_leapfrog

    subroutine Global_step(this, iseed, Counter, stage_id)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer, intent(in) :: stage_id
        class(AccCounter), intent(inout) :: Counter
        real(kind=8) :: action_new, ham_old, ham_new, delta_h, ratio, random
        integer :: nsteps
        logical :: accepted
        real(kind=8), external :: ranf

        call this%ensure_state()
        this%phi_backup = Conf%phi_list
        call this%select_site_block(iseed)
        call this%select_tau_block(iseed)
        call this%sample_momentum(iseed)
        ham_old = 0.5d0 * sum(this%momentum * this%momentum) / hmc_mass + this%action_cur

        nsteps = HMC_draw_nfrog(iseed)
        call this%leapfrog(action_new, nsteps, stage_id)
        ham_new = 0.5d0 * sum(this%momentum * this%momentum) / hmc_mass + action_new
        if (.not. ieee_is_finite(this%action_cur) .or. .not. ieee_is_finite(action_new) .or. &
            .not. ieee_is_finite(ham_old) .or. .not. ieee_is_finite(ham_new)) then
            write(6,*) 'Non-finite HMC energy detected at rank=', IRANK
            write(6,*) '  action_cur=', this%action_cur, ' action_new=', action_new
            write(6,*) '  ham_old=', ham_old, ' ham_new=', ham_new
            write(6,*) '  block sites=', this%ii_begin, this%ii_end, ' block tau=', this%tau_begin, this%tau_end
            write(6,*) '  nsteps=', nsteps, ' hmc_dt=', hmc_dt, ' hmc_mass=', hmc_mass
            stop 1
        endif

        delta_h = ham_old - ham_new
        HMC_deltaH_sum = HMC_deltaH_sum + delta_h
        HMC_deltaH_absmax = max(HMC_deltaH_absmax, abs(delta_h))
        HMC_deltaH_count = HMC_deltaH_count + 1
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
        call HMC_monitor_step(stage_id, accepted, nsteps, delta_h, this%action_cur)
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
        integer, intent(out) :: Nobs, Nobst
        integer :: nt

        Nobs = 0
        Nobst = 0
        if (HMC_debug_measure > 0 .and. IRANK == 0) then
            write(6,*) 'Global_measure_config start: Nobs=', Nobs, ' Nobst=', Nobst, ' Ltrot=', Ltrot
        endif
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
        if (HMC_debug_measure > 0 .and. IRANK == 0) then
            write(6,*) 'Global_measure_config after sweep_L: Nobs=', Nobs
        endif
        if (toggle .and. is_tau .and. allocated(Obs_tau_hmc)) then
            call HMC_set_overlap(this%prop_work)
            call Dyn%reset(this%prop_work)
            call Dyn%sweep_R(Obs_tau_hmc, this%wr_work)
            Nobst = Nobst + 1
        endif
        call HMC_measure_sweep_R(this%prop_work, this%wr_work, Nobs)
        if (HMC_debug_measure > 0 .and. IRANK == 0) then
            write(6,*) 'Global_measure_config after sweep_R: Nobs=', Nobs, ' Nobst=', Nobst
        endif
        return
    end subroutine Global_measure_config

    subroutine Global_check_equal_obs(this, stage, Nobs, Nobst, toggle)
        class(GlobalUpdate), intent(in) :: this
        character(len=*), intent(in) :: stage
        integer, intent(in) :: Nobs, Nobst
        logical, intent(in) :: toggle

        if (Nobs <= 0) then
            write(6,*) 'Global_check_equal_obs: no equal-time samples in stage ', trim(stage)
            write(6,*) '  Nobs=', Nobs, ' Nobst=', Nobst, ' Nsweep=', Nsweep, ' toggle=', toggle
            write(6,*) '  block sites=', this%ii_begin, this%ii_end, ' block tau=', this%tau_begin, this%tau_end
            stop 1
        endif
        if (.not. ieee_is_finite(Obs_equal_hmc%density_up) .or. .not. ieee_is_finite(Obs_equal_hmc%density_do) .or. &
            .not. ieee_is_finite(Obs_equal_hmc%kinetic) .or. .not. ieee_is_finite(Obs_equal_hmc%doubleOcc) .or. &
            .not. ieee_is_finite(Obs_equal_hmc%squareOcc) .or. .not. ieee_is_finite(Obs_equal_hmc%IPR) .or. &
            .not. ieee_is_finite(Obs_equal_hmc%nearestOcc) .or. .not. ieee_is_finite(Obs_equal_hmc%num_up) .or. &
            .not. ieee_is_finite(Obs_equal_hmc%num_do) .or. .not. ieee_is_finite(Obs_equal_hmc%numsquare_up) .or. &
            .not. ieee_is_finite(Obs_equal_hmc%numsquare_do)) then
            write(6,*) 'Global_check_equal_obs: non-finite equal-time scalar in stage ', trim(stage)
            write(6,*) '  Nobs=', Nobs, ' Nobst=', Nobst, ' toggle=', toggle
            write(6,*) '  density_up=', Obs_equal_hmc%density_up, ' density_do=', Obs_equal_hmc%density_do
            write(6,*) '  kinetic=', Obs_equal_hmc%kinetic, ' doubleOcc=', Obs_equal_hmc%doubleOcc
            write(6,*) '  squareOcc=', Obs_equal_hmc%squareOcc, ' IPR=', Obs_equal_hmc%IPR
            write(6,*) '  nearestOcc=', Obs_equal_hmc%nearestOcc
            write(6,*) '  block sites=', this%ii_begin, this%ii_end, ' block tau=', this%tau_begin, this%tau_end
            stop 1
        endif
        return
    end subroutine Global_check_equal_obs

    subroutine Global_therm(this, iseed)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed

        call Acc_HMC_warm%reset()
        call this%step(iseed, Acc_HMC_warm, 0)
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
        integer :: Nobs, Nobst, nobs_cfg, nobst_cfg, nsw

        call this%reset(toggle)
        Nobs = 0
        Nobst = 0
        do nsw = 1, Nsweep
            if (HMC_debug_measure > 0 .and. IRANK == 0) then
                write(6,*) 'Global_sweep before step: sweep=', nsw, ' Nobs=', Nobs, ' Nobst=', Nobst
            endif
            call this%step(iseed, Acc_HMC, 1)
            if (HMC_debug_measure > 0 .and. IRANK == 0) then
                write(6,*) 'Global_sweep after step: sweep=', nsw, ' Nobs=', Nobs, ' Nobst=', Nobst
            endif
            call this%measure_config(toggle, nobs_cfg, nobst_cfg)
            Nobs = Nobs + nobs_cfg
            Nobst = Nobst + nobst_cfg
            if (HMC_debug_measure > 0 .and. IRANK == 0) then
                write(6,*) 'Global_sweep after measure: sweep=', nsw, ' Nobs=', Nobs, ' Nobst=', Nobst
            endif
        enddo

        call Obs_equal_hmc%ave(Nobs)
        call this%check_equal_obs('global_sweep', Nobs, Nobst, toggle)
        if (toggle .and. is_tau .and. allocated(Obs_tau_hmc)) call Obs_tau_hmc%ave(Nobst)
        call Acc_HMC%ratio()
        return
    end subroutine Global_sweep

    subroutine Global_control_print(toggle)
        include 'mpif.h'
        logical, intent(in) :: toggle
        real(kind=8) :: collect, collect_sum, collect_absmax
        integer :: collect_count

        collect = 0.d0
        call MPI_Reduce(Acc_HMC%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        collect_sum = 0.d0
        call MPI_Reduce(HMC_deltaH_sum, collect_sum, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        collect_absmax = 0.d0
        call MPI_Reduce(HMC_deltaH_absmax, collect_absmax, 1, MPI_Real8, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
        collect_count = 0
        call MPI_Reduce(HMC_deltaH_count, collect_count, 1, MPI_Integer, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            Acc_HMC%acc = collect / dble(ISIZE * Nbin)
            write(50,*) 'Accept_HMC                                     :', Acc_HMC%acc
            if (collect_count > 0) then
                write(50,*) 'HMC DeltaH mean                                :', collect_sum / dble(collect_count)
                write(50,*) 'HMC DeltaH abs max                             :', collect_absmax
            endif
        endif
        if (toggle .and. is_tau) call Dyn%ctrl_print()
        return
    end subroutine Global_control_print

    subroutine Global_therm_print()
        include 'mpif.h'
        real(kind=8) :: collect

        collect = 0.d0
        call MPI_Reduce(Acc_HMC_warm%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            if (Nwarm > 0) then
                Acc_HMC_warm%acc = collect / dble(ISIZE * Nwarm)
            else
                Acc_HMC_warm%acc = 0.d0
            endif
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
                if (HMC_debug_measure > 0 .and. IRANK == 0) then
                    write(6,*) 'HMC_measure_sweep_L hit ntau=', nt, ' Nobs before=', Nobs
                endif
                call HMC_set_overlap(Prop)
                call Obs_equal_hmc%calc(Prop, nt)
                Nobs = Nobs + 1
                if (HMC_debug_measure > 0 .and. IRANK == 0) then
                    write(6,*) 'HMC_measure_sweep_L after calc: Nobs=', Nobs, ' doubleOcc=', Obs_equal_hmc%doubleOcc
                endif
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
                if (HMC_debug_measure > 0 .and. IRANK == 0) then
                    write(6,*) 'HMC_measure_sweep_R hit ntau=', nt, ' Nobs before=', Nobs
                endif
                call HMC_set_overlap(Prop)
                call Obs_equal_hmc%calc(Prop, nt)
                Nobs = Nobs + 1
                if (HMC_debug_measure > 0 .and. IRANK == 0) then
                    write(6,*) 'HMC_measure_sweep_R after calc: Nobs=', Nobs, ' doubleOcc=', Obs_equal_hmc%doubleOcc
                endif
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

    subroutine HMC_monitor_init()
        character(len=32) :: env_value
        integer :: env_status

        HMC_monitor_enabled = .false.
        HMC_monitor_count = 0
        HMC_monitor_site = 1
        HMC_monitor_tau = max(1, Ltrot / 2)

        call get_environment_variable("BPQMC_HMC_MONITOR", env_value, status=env_status)
        if (env_status /= 0) return
        if (len_trim(env_value) == 0) return
        if (trim(env_value) == "0") return

        call HMC_read_env_int("BPQMC_HMC_MONITOR_SITE", HMC_monitor_site)
        call HMC_read_env_int("BPQMC_HMC_MONITOR_TAU", HMC_monitor_tau)
        HMC_monitor_site = min(max(HMC_monitor_site, 1), Ndim)
        HMC_monitor_tau = min(max(HMC_monitor_tau, 1), max(Ltrot, 1))
        HMC_monitor_enabled = .true.
        open(unit=HMC_monitor_unit, file='hmc_monitor.dat', status='replace', action='write')
        write(HMC_monitor_unit, '(A)') '# step stage accepted nsteps delta_h action phi_f1 phi_f2 phi_rms_f1 phi_rms_f2'
        return
    end subroutine HMC_monitor_init

    subroutine HMC_monitor_close()
        if (.not. HMC_monitor_enabled) return
        close(unit=HMC_monitor_unit)
        HMC_monitor_enabled = .false.
        return
    end subroutine HMC_monitor_close

    subroutine HMC_monitor_step(stage_id, accepted, nsteps, delta_h, action_cur)
        integer, intent(in) :: stage_id, nsteps
        logical, intent(in) :: accepted
        real(kind=8), intent(in) :: delta_h, action_cur
        integer :: accepted_flag
        real(kind=8) :: phi_f1, phi_f2, phi_rms_f1, phi_rms_f2

        if (.not. HMC_monitor_enabled) return
        if (Ltrot <= 0) return

        HMC_monitor_count = HMC_monitor_count + 1
        accepted_flag = merge(1, 0, accepted)
        call HMC_trace_values(Conf%phi_list, phi_f1, phi_f2, phi_rms_f1, phi_rms_f2)
        write(HMC_monitor_unit, '(I0,1X,I0,1X,I0,1X,I0,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)') &
            HMC_monitor_count, stage_id, accepted_flag, nsteps, delta_h, action_cur, phi_f1, phi_f2, phi_rms_f1, phi_rms_f2
        return
    end subroutine HMC_monitor_step

    subroutine HMC_trace_init()
        character(len=32) :: env_value
        integer :: env_status

        HMC_trace_enabled = .false.
        HMC_trace_count = 0
        HMC_trace_proposal = 0
        HMC_trace_stage = 0
        HMC_trace_nsteps = 0

        call get_environment_variable("BPQMC_HMC_TRACE", env_value, status=env_status)
        if (env_status /= 0) return
        if (len_trim(env_value) == 0) return
        if (trim(env_value) == "0") return

        call HMC_read_env_int("BPQMC_HMC_MONITOR_SITE", HMC_monitor_site)
        call HMC_read_env_int("BPQMC_HMC_MONITOR_TAU", HMC_monitor_tau)
        HMC_monitor_site = min(max(HMC_monitor_site, 1), Ndim)
        HMC_monitor_tau = min(max(HMC_monitor_tau, 1), max(Ltrot, 1))
        HMC_trace_enabled = .true.
        open(unit=HMC_trace_unit, file='hmc_trace.dat', status='replace', action='write')
        write(HMC_trace_unit, '(A)') '# proposal stage lf_step nsteps md_time action phi_f1 phi_f2 mom_f1 mom_f2 force_f1 force_f2 phi_rms_f1 phi_rms_f2'
        return
    end subroutine HMC_trace_init

    subroutine HMC_trace_close()
        if (.not. HMC_trace_enabled) return
        close(unit=HMC_trace_unit)
        HMC_trace_enabled = .false.
        return
    end subroutine HMC_trace_close

    subroutine HMC_trace_begin(stage_id, nsteps, action_cur, momentum, force)
        integer, intent(in) :: stage_id, nsteps
        real(kind=8), intent(in) :: action_cur
        real(kind=8), dimension(Naux, Ndim, Ltrot), intent(in) :: momentum, force

        if (.not. HMC_trace_enabled) return
        if (Ltrot <= 0) return

        HMC_trace_proposal = HMC_trace_proposal + 1
        HMC_trace_stage = stage_id
        HMC_trace_nsteps = nsteps
        call HMC_trace_step(0, action_cur, momentum, force)
        return
    end subroutine HMC_trace_begin

    subroutine HMC_trace_step(lf_step, action_cur, momentum, force)
        integer, intent(in) :: lf_step
        real(kind=8), intent(in) :: action_cur
        real(kind=8), dimension(Naux, Ndim, Ltrot), intent(in) :: momentum, force
        integer :: ntau_idx
        real(kind=8) :: md_time
        real(kind=8) :: phi_f1, phi_f2, phi_rms_f1, phi_rms_f2
        real(kind=8) :: mom_f1, mom_f2, force_f1, force_f2

        if (.not. HMC_trace_enabled) return
        if (Ltrot <= 0) return

        ntau_idx = min(max(HMC_monitor_tau, 1), Ltrot)
        md_time = dble(lf_step) * hmc_dt
        call HMC_trace_values(Conf%phi_list, phi_f1, phi_f2, phi_rms_f1, phi_rms_f2)
        mom_f1 = momentum(1, HMC_monitor_site, ntau_idx)
        force_f1 = force(1, HMC_monitor_site, ntau_idx)
        if (Naux >= 2) then
            mom_f2 = momentum(2, HMC_monitor_site, ntau_idx)
            force_f2 = force(2, HMC_monitor_site, ntau_idx)
        else
            mom_f2 = 0.d0
            force_f2 = 0.d0
        endif
        HMC_trace_count = HMC_trace_count + 1
        write(HMC_trace_unit, '(I0,1X,I0,1X,I0,1X,I0,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)') &
            HMC_trace_proposal, HMC_trace_stage, lf_step, HMC_trace_nsteps, md_time, action_cur, phi_f1, phi_f2, mom_f1, mom_f2, force_f1, force_f2, phi_rms_f1, phi_rms_f2
        return
    end subroutine HMC_trace_step

    subroutine HMC_trace_values(phi_list, phi_f1, phi_f2, phi_rms_f1, phi_rms_f2)
        real(kind=8), dimension(Naux, Ndim, Ltrot), intent(in) :: phi_list
        real(kind=8), intent(out) :: phi_f1, phi_f2, phi_rms_f1, phi_rms_f2
        integer :: ntau_idx

        ntau_idx = min(max(HMC_monitor_tau, 1), Ltrot)
        phi_f1 = phi_list(1, HMC_monitor_site, ntau_idx)
        phi_rms_f1 = sqrt(sum(phi_list(1,:,:) * phi_list(1,:,:)) / dble(Ndim * Ltrot))
        if (Naux >= 2) then
            phi_f2 = phi_list(2, HMC_monitor_site, ntau_idx)
            phi_rms_f2 = sqrt(sum(phi_list(2,:,:) * phi_list(2,:,:)) / dble(Ndim * Ltrot))
        else
            phi_f2 = 0.d0
            phi_rms_f2 = 0.d0
        endif
        return
    end subroutine HMC_trace_values

    subroutine HMC_read_env_int(name, value)
        character(len=*), intent(in) :: name
        integer, intent(inout) :: value
        character(len=32) :: env_value
        integer :: env_status, ios

        call get_environment_variable(name, env_value, status=env_status)
        if (env_status /= 0) return
        if (len_trim(env_value) == 0) return
        read(env_value, *, iostat=ios) value
        if (ios /= 0) return
        return
    end subroutine HMC_read_env_int
end module GlobalUpdate_mod
