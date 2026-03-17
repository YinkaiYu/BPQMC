program check_hmc_ratio_scan
    use CalcBasic
    use DQMC_Model_mod
    use GlobalUpdate_mod
    use Multiply_mod
    use ProcessMatrix
    use Stabilize_mod
    implicit none
    include 'mpif.h'

    type(GlobalUpdate) :: sweep_global
    type(Propagator), allocatable :: prop
    type(WrapList), allocatable :: wr_list
    real(kind=8), allocatable :: force(:,:,:)
    real(kind=8) :: action_old, action_new
    real(kind=8) :: phi_old, phi_new, ratio_action, ratio_local
    real(kind=8) :: rel_diff, max_rel_diff, sum_rel_diff, shift_scale, shift_delta
    complex(kind=8) :: overlap, delta_term
    integer :: iseed, ntrial, itrial, ii, ntau, nf, active_flavor_count
    integer :: flavor_list(Naux), worst_nf, worst_ii, worst_ntau
    logical :: ok
    complex(kind=8), external :: ZDOTU
    real(kind=8), external :: ranf

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)

    call Model_init(iseed)
    call Stabilize_init()
    call sweep_global%init(Init)
    allocate(force(Naux, Ndim, Ltrot))
    allocate(prop)
    call prop%make(Init)
    allocate(wr_list)
    call wr_list%make()

    ntrial = 32
    shift_scale = 0.2d0
    call read_optional_int("BPQMC_RATIO_SCAN_TRIALS", ntrial)
    call read_optional_real("BPQMC_RATIO_SCAN_SHIFT", shift_scale)

    call sweep_global%debug_eval(action_old, force, ok)
    if (.not. ok) then
        if (IRANK == 0) write(6,*) 'debug_eval failed on original configuration'
        stop 1
    endif

    active_flavor_count = 0
    if (abs(RU1) > Zero) then
        active_flavor_count = active_flavor_count + 1
        flavor_list(active_flavor_count) = 1
    endif
    if (abs(RU2) > Zero) then
        active_flavor_count = active_flavor_count + 1
        flavor_list(active_flavor_count) = 2
    endif
    if (active_flavor_count == 0) then
        if (IRANK == 0) write(6,*) 'No active Hubbard flavor for ratio scan.'
        stop 1
    endif

    max_rel_diff = 0.d0
    sum_rel_diff = 0.d0
    worst_nf = 1
    worst_ii = 1
    worst_ntau = 1

    do itrial = 1, ntrial
        nf = flavor_list(nranf(iseed, active_flavor_count))
        ii = nranf(iseed, Ndim)
        ntau = nranf(iseed, Ltrot)
        shift_delta = shift_scale * (ranf(iseed) - 0.5d0)
        if (abs(shift_delta) < 1.d-12) cycle

        call build_local_context(prop, wr_list, nf, ntau)
        overlap = ZDOTU(Ndim, prop%UUL(1,1), 1, prop%UUR(1,1), 1)

        phi_old = Conf%phi_list(nf, ii, ntau)
        phi_new = phi_old + shift_delta
        if (nf == 1) then
            call Op_U1%get_delta(phi_old, phi_new)
            delta_term = Op_U1%Delta * prop%UUR(ii,1) * prop%UUL(1,ii) / overlap
        else
            call Op_U2%get_delta(phi_old, phi_new)
            delta_term = Op_U2%Delta * prop%UUR(ii,1) * prop%UUL(1,ii) / overlap
        endif
        ratio_local = exp(-0.5d0 * (phi_new * phi_new - phi_old * phi_old)) * abs(1.d0 + delta_term) ** (2 * Nbos)

        Conf%phi_list(nf, ii, ntau) = phi_new
        call sweep_global%debug_eval(action_new, force, ok)
        if (.not. ok) then
            if (IRANK == 0) write(6,*) 'debug_eval failed on shifted configuration'
            stop 1
        endif
        ratio_action = exp(-(action_new - action_old))
        rel_diff = abs(ratio_local - ratio_action) / max(abs(ratio_action), 1.d-300)
        sum_rel_diff = sum_rel_diff + rel_diff
        if (rel_diff > max_rel_diff) then
            max_rel_diff = rel_diff
            worst_nf = nf
            worst_ii = ii
            worst_ntau = ntau
        endif

        Conf%phi_list(nf, ii, ntau) = phi_old
    enddo

    if (IRANK == 0) then
        write(6,'(A,I0)') 'ratio_scan_trials = ', ntrial
        write(6,'(A,ES24.16)') 'ratio_scan_shift = ', shift_scale
        write(6,'(A,ES24.16)') 'ratio_scan_mean_reldiff = ', sum_rel_diff / dble(max(ntrial, 1))
        write(6,'(A,ES24.16)') 'ratio_scan_max_reldiff = ', max_rel_diff
        write(6,'(A,I0,A,I0,A,I0)') 'ratio_scan_worst nf=', worst_nf, ' ii=', worst_ii, ' ntau=', worst_ntau
    endif

    call sweep_global%clear()
    call Stabilize_clear()
    deallocate(force)
    deallocate(prop)
    deallocate(wr_list)
    call Model_clear(iseed)
    call MPI_FINALIZE(IERR)

contains
    subroutine build_local_context(prop_ctx, wr_ctx, nf_target, ntau_target)
        class(Propagator), intent(inout) :: prop_ctx
        class(WrapList), intent(inout) :: wr_ctx
        integer, intent(in) :: nf_target, ntau_target
        integer :: nt

        call prop_ctx%reset(Init)
        call wr_ctx%reset()
        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call propU_pre(Op_U1, prop_ctx, 1, nt)
            if (abs(RU2) > Zero) call propU_pre(Op_U2, prop_ctx, 2, nt)
            call propT_pre(prop_ctx)
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_pre(prop_ctx, wr_ctx, nt)
        enddo

        do nt = Ltrot, ntau_target + 1, -1
            if (mod(nt, Nwrap) == 0 .or. nt == Ltrot) call Wrap_L(prop_ctx, wr_ctx, nt, "H")
            call propT_L(prop_ctx)
            if (abs(RU2) > Zero) call propU_L(Op_U2, prop_ctx, 2, nt)
            if (abs(RU1) > Zero) call propU_L(Op_U1, prop_ctx, 1, nt)
        enddo

        if (mod(ntau_target, Nwrap) == 0 .or. ntau_target == Ltrot) call Wrap_L(prop_ctx, wr_ctx, ntau_target, "H")
        call propT_L(prop_ctx)
        if (nf_target == 1 .and. abs(RU2) > Zero) call propU_L(Op_U2, prop_ctx, 2, ntau_target)
        return
    end subroutine build_local_context

    subroutine read_optional_int(name, value)
        character(len=*), intent(in) :: name
        integer, intent(inout) :: value
        character(len=32) :: env_value
        integer :: env_status, ios

        call get_environment_variable(name, env_value, status=env_status)
        if (env_status /= 0) return
        if (len_trim(env_value) == 0) return
        read(env_value, *, iostat=ios) value
        if (ios /= 0) return
    end subroutine read_optional_int

    subroutine read_optional_real(name, value)
        character(len=*), intent(in) :: name
        real(kind=8), intent(inout) :: value
        character(len=32) :: env_value
        integer :: env_status, ios

        call get_environment_variable(name, env_value, status=env_status)
        if (env_status /= 0) return
        if (len_trim(env_value) == 0) return
        read(env_value, *, iostat=ios) value
        if (ios /= 0) return
    end subroutine read_optional_real
end program check_hmc_ratio_scan
