program check_hmc_force_fd
    use CalcBasic
    use DQMC_Model_mod
    use GlobalUpdate_mod
    use Stabilize_mod
    implicit none
    include 'mpif.h'

    type(GlobalUpdate) :: sweep_global
    real(kind=8), allocatable :: force_ref(:,:,:), force_work(:,:,:)
    real(kind=8) :: action_ref, action_plus, action_minus
    real(kind=8) :: eps_scale, eps, phi_old, force_analytic, force_fd
    real(kind=8) :: abs_diff, rel_diff
    real(kind=8) :: max_abs_diff, max_rel_diff, sum_abs_diff, sum_rel_diff
    real(kind=8) :: worst_force_analytic, worst_force_fd, worst_eps
    integer :: iseed, ntrial, itrial, ii, ntau, nf, active_flavor_count
    integer :: flavor_list(Naux), worst_nf, worst_ii, worst_ntau
    logical :: ok

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)

    call Model_init(iseed)
    call Stabilize_init()
    call sweep_global%init(Init)

    allocate(force_ref(Naux, Ndim, Ltrot))
    allocate(force_work(Naux, Ndim, Ltrot))

    ntrial = 32
    eps_scale = 1.d-6
    call read_optional_int("BPQMC_FORCE_FD_TRIALS", ntrial)
    call read_optional_real("BPQMC_FORCE_FD_EPS", eps_scale)

    call sweep_global%debug_eval(action_ref, force_ref, ok)
    if (.not. ok) then
        if (IRANK == 0) write(6,*) 'debug_eval failed on reference configuration'
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
        if (IRANK == 0) write(6,*) 'No active Hubbard flavor for force finite-difference check.'
        stop 1
    endif

    max_abs_diff = 0.d0
    max_rel_diff = 0.d0
    sum_abs_diff = 0.d0
    sum_rel_diff = 0.d0
    worst_nf = 1
    worst_ii = 1
    worst_ntau = 1
    worst_force_analytic = 0.d0
    worst_force_fd = 0.d0
    worst_eps = 0.d0

    do itrial = 1, ntrial
        nf = flavor_list(nranf(iseed, active_flavor_count))
        ii = nranf(iseed, Ndim)
        ntau = nranf(iseed, Ltrot)
        phi_old = Conf%phi_list(nf, ii, ntau)
        eps = eps_scale * max(1.d0, abs(phi_old))

        Conf%phi_list(nf, ii, ntau) = phi_old + eps
        call sweep_global%debug_eval(action_plus, force_work, ok)
        if (.not. ok) then
            if (IRANK == 0) write(6,*) 'debug_eval failed on plus-shift configuration'
            stop 1
        endif

        Conf%phi_list(nf, ii, ntau) = phi_old - eps
        call sweep_global%debug_eval(action_minus, force_work, ok)
        if (.not. ok) then
            if (IRANK == 0) write(6,*) 'debug_eval failed on minus-shift configuration'
            stop 1
        endif

        Conf%phi_list(nf, ii, ntau) = phi_old

        force_analytic = force_ref(nf, ii, ntau)
        force_fd = -(action_plus - action_minus) / (2.d0 * eps)
        abs_diff = abs(force_analytic - force_fd)
        rel_diff = abs_diff / max(abs(force_fd), 1.d-300)
        sum_abs_diff = sum_abs_diff + abs_diff
        sum_rel_diff = sum_rel_diff + rel_diff
        if (abs_diff > max_abs_diff) then
            max_abs_diff = abs_diff
            worst_nf = nf
            worst_ii = ii
            worst_ntau = ntau
            worst_force_analytic = force_analytic
            worst_force_fd = force_fd
            worst_eps = eps
        endif
        if (rel_diff > max_rel_diff) max_rel_diff = rel_diff
    enddo

    call sweep_global%clear()
    call Stabilize_clear()
    deallocate(force_ref)
    deallocate(force_work)
    call Model_clear(iseed)

    if (IRANK == 0) then
        write(6,'(A,I0)') 'force_fd_trials = ', ntrial
        write(6,'(A,ES24.16)') 'force_fd_eps_scale = ', eps_scale
        write(6,'(A,ES24.16)') 'force_fd_mean_absdiff = ', sum_abs_diff / dble(max(ntrial, 1))
        write(6,'(A,ES24.16)') 'force_fd_max_absdiff = ', max_abs_diff
        write(6,'(A,ES24.16)') 'force_fd_mean_reldiff = ', sum_rel_diff / dble(max(ntrial, 1))
        write(6,'(A,ES24.16)') 'force_fd_max_reldiff = ', max_rel_diff
        write(6,'(A,I0,A,I0,A,I0)') 'force_fd_worst nf=', worst_nf, ' ii=', worst_ii, ' ntau=', worst_ntau
        write(6,'(A,ES24.16)') 'force_fd_worst_analytic = ', worst_force_analytic
        write(6,'(A,ES24.16)') 'force_fd_worst_fd = ', worst_force_fd
        write(6,'(A,ES24.16)') 'force_fd_worst_eps = ', worst_eps
    endif

    call MPI_FINALIZE(IERR)

contains
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
end program check_hmc_force_fd
