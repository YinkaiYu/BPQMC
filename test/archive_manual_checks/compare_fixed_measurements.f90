program compare_fixed_measurements
    use CalcBasic
    use DQMC_Model_mod
    use GlobalUpdate_mod
    use LocalSweep_mod
    use ProcessMatrix
    use Stabilize_mod
    implicit none
    include 'mpif.h'

    type(GlobalUpdate) :: sweep_global
    type(LocalSweep) :: sweep_local
    type(Propagator), allocatable :: prop
    type(WrapList), allocatable :: wr_list
    integer :: iseed
    logical :: is_beta
    real(kind=8) :: shift_save

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)

    call Model_init(iseed)
    call Stabilize_init()
    allocate(prop)
    call prop%make(Init)
    allocate(wr_list)
    call wr_list%make()
    call sweep_local%init(Init)
    call sweep_global%init(Init)

    shift_save = shiftLoc
    shiftLoc = 0.d0
    is_beta = .true.
    call sweep_local%pre(prop, wr_list)
    call sweep_local%sweep(prop, wr_list, iseed, is_beta, .false.)
    shiftLoc = shift_save

    call sweep_global%debug_measure_current(.false.)

    if (IRANK == 0) then
        write(6,'(A)') 'Fixed-configuration local vs HMC measurement comparison'
        call print_scalar('density_up', Obs_equal%density_up, Obs_equal_hmc%density_up)
        call print_scalar('density_do', Obs_equal%density_do, Obs_equal_hmc%density_do)
        call print_scalar('kinetic', Obs_equal%kinetic, Obs_equal_hmc%kinetic)
        call print_scalar('doubleOcc', Obs_equal%doubleOcc, Obs_equal_hmc%doubleOcc)
        call print_scalar('squareOcc', Obs_equal%squareOcc, Obs_equal_hmc%squareOcc)
        call print_scalar('nearestOcc', Obs_equal%nearestOcc, Obs_equal_hmc%nearestOcc)
        call print_scalar('IPR', Obs_equal%IPR, Obs_equal_hmc%IPR)
        call print_scalar('num_up', Obs_equal%num_up, Obs_equal_hmc%num_up)
        call print_scalar('num_do', Obs_equal%num_do, Obs_equal_hmc%num_do)
        call print_scalar('numsquare_up', Obs_equal%numsquare_up, Obs_equal_hmc%numsquare_up)
        call print_scalar('numsquare_do', Obs_equal%numsquare_do, Obs_equal_hmc%numsquare_do)
        call print_complex_max1('SF_corr_up', Obs_equal%SF_corr_up, Obs_equal_hmc%SF_corr_up)
        call print_complex_max1('SF_corr_do', Obs_equal%SF_corr_do, Obs_equal_hmc%SF_corr_do)
        call print_complex_max1('PF_corr', Obs_equal%PF_corr, Obs_equal_hmc%PF_corr)
        call print_complex_max1('C3_corr', Obs_equal%C3_corr, Obs_equal_hmc%C3_corr)
        call print_complex_max3('den_corr_up', Obs_equal%den_corr_up, Obs_equal_hmc%den_corr_up)
        call print_complex_max3('den_corr_do', Obs_equal%den_corr_do, Obs_equal_hmc%den_corr_do)
        call print_complex_max3('single_corr', Obs_equal%single_corr, Obs_equal_hmc%single_corr)
    endif

    call sweep_global%clear()
    call sweep_local%clear()
    call Stabilize_clear()
    deallocate(prop)
    deallocate(wr_list)
    call Model_clear(iseed)
    call MPI_FINALIZE(IERR)

contains
    subroutine print_scalar(name, local_value, hmc_value)
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: local_value, hmc_value
        write(6,'(A,1X,ES24.16,1X,ES24.16,1X,ES24.16)') trim(name), local_value, hmc_value, abs(local_value - hmc_value)
        return
    end subroutine print_scalar

    subroutine print_complex_max1(name, local_value, hmc_value)
        character(len=*), intent(in) :: name
        complex(kind=8), intent(in) :: local_value(:), hmc_value(:)
        write(6,'(A,1X,ES24.16)') trim(name), maxval(abs(local_value - hmc_value))
    end subroutine print_complex_max1

    subroutine print_complex_max3(name, local_value, hmc_value)
        character(len=*), intent(in) :: name
        complex(kind=8), intent(in) :: local_value(:,:,:), hmc_value(:,:,:)
        write(6,'(A,1X,ES24.16)') trim(name), maxval(abs(local_value - hmc_value))
    end subroutine print_complex_max3
end program compare_fixed_measurements
