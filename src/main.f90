program BPQMC
    use LocalSweep_mod
    use GlobalUpdate_mod
    use FourierTrans_mod
    implicit none
    include 'mpif.h'
    
    integer :: status(MPI_STATUS_SIZE)
    integer:: iseed, nth, nbc, N
    logical :: is_beta, istau_tmp
    real(kind=8) :: collect, CPUT, wall_1, wall_2
    
    type(GlobalUpdate) :: Sweep_global
    type(LocalSweep) :: Sweep_local
    type(FourierTrans) :: Fourier
    type(Propagator), allocatable :: Prop
    type(WrapList), allocatable :: WrList

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
    
    wall_1 = MPI_WTIME()
! initiate
    call Model_init(iseed)
    allocate(Prop)
    call Prop%make(Init)
    allocate(WrList)
    call WrList%make()
    call Stabilize_init()
    if (is_global) then
        call Sweep_global%init(Init)
    else
        call Sweep_local%init(Init)
    endif
! boson warm-up
    if (is_warm) then
        do nth = 1, Nwarm
            if (is_global) then
                call Sweep_global%therm(iseed)
            else
                call Sweep_local%therm(iseed)
            endif
        enddo
        if (is_global) then
            call Sweep_global%ctrl_print_t()
        else
            call Sweep_local%ctrl_print_t()
        endif
    else
        if (IRANK == 0) write(50,*) "Skipping Bosonic warm-up"
    endif
! Sweep
    if (.not. is_global) call Sweep_local%pre(Prop, WrList)
    is_beta = .true.; istau_tmp = .false.
    do nbc = 1, Nbin
        if (nbc .gt. Nthermal) istau_tmp = is_tau
        if (is_global) then
            call Sweep_global%sweep(Prop, WrList, iseed, is_beta, istau_tmp)
            call Fourier%preq(Obs_equal_hmc)
        else
            call Sweep_local%sweep(Prop, WrList, iseed, is_beta, istau_tmp)
            call Fourier%preq(Obs_equal)
        endif
        if (istau_tmp) then
            if (is_global) then
                call Fourier%prtau(Obs_tau_hmc)
            else
                call Fourier%prtau(Obs_tau)
            endif
        endif
    enddo
! control print
    if (is_global) then
        call Sweep_global%ctrl_print(istau_tmp)
    else
        call Sweep_local%ctrl_print_l(istau_tmp)
        collect = 0.d0
        call MPI_Reduce(Prop%Xmaxm, collect, 1, MPI_Real8, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Prop%Xmaxm = collect
        N = 2* ISIZE * Nbin * Nst * Nsweep
        collect = 0.d0
        call MPI_Reduce(Prop%Xmeanm, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Prop%Xmeanm = collect / dble(N)
    endif
    wall_2 = MPI_WTIME()
    CPUT = wall_2 - wall_1
    collect = 0.d0
    call MPI_Reduce(CPUT, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) CPUT = collect / dble(ISIZE)
    if (IRANK == 0) then
        if (.not. is_global) then
            write(50,*) 'Max diff Matrix                                :', Prop%Xmaxm
            write(50,*) 'Mean diff Matrix                               :', Prop%Xmeanm
        endif
        write(50,*) 'Tot CPU time                                   :', CPUT
    endif
! deallocate
    if (is_global) then
        call Sweep_global%clear()
    else
        call Sweep_local%clear()
    endif
    call Stabilize_clear()
    deallocate(Prop)
    deallocate(WrList)
    call Model_clear(iseed) ! conf-out
    
    call MPI_FINALIZE(IERR)
end program BPQMC
