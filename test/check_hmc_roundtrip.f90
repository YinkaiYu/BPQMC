program check_hmc_roundtrip
    use CalcBasic
    use DQMC_Model_mod
    use GlobalUpdate_mod
    use Stabilize_mod
    implicit none
    include 'mpif.h'

    type(GlobalUpdate) :: sweep_global
    real(kind=8) :: diff_phi, diff_momentum, diff_action
    integer :: iseed

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)

    call Model_init(iseed)
    call Stabilize_init()
    call sweep_global%init(Init)
    call sweep_global%debug_roundtrip(iseed, Nfrog, diff_phi, diff_momentum, diff_action)

    if (IRANK == 0) then
        write(6,'(A,I0)') 'nsteps       = ', Nfrog
        write(6,'(A,ES24.16)') 'diff_phi     = ', diff_phi
        write(6,'(A,ES24.16)') 'diff_momentum= ', diff_momentum
        write(6,'(A,ES24.16)') 'diff_action  = ', diff_action
    endif

    call sweep_global%clear()
    call Stabilize_clear()
    call Model_clear(iseed)
    call MPI_FINALIZE(IERR)
end program check_hmc_roundtrip
