module FourierTrans_mod
    use DQMC_Model_mod
    use ObserEqual_mod
    use ObserTau_mod
    implicit none
    
    type, public :: FourierTrans
    contains
        procedure, private, nopass :: m_write_real_1
        procedure, private, nopass :: m_write_real_2
        generic :: write_real => m_write_real_1
        generic :: write_real => m_write_real_2
        
        ! procedure, private, nopass :: write_cmplx => m_write_cmplx_3

        procedure, private, nopass :: m_write_cmplx_1
        procedure, private, nopass :: m_write_cmplx_3
        generic :: write_cmplx => m_write_cmplx_1
        generic :: write_cmplx => m_write_cmplx_3
        
        procedure, private, nopass :: m_write_reciprocal_1
        procedure, private, nopass :: m_write_reciprocal_2
        procedure, private, nopass :: m_write_reciprocal_3
        generic :: write_reciprocal => m_write_reciprocal_1
        generic :: write_reciprocal => m_write_reciprocal_2
        generic :: write_reciprocal => m_write_reciprocal_3
        
        procedure, private, nopass :: m_write_k_1
        procedure, private, nopass :: m_write_k_2
        procedure, private, nopass :: m_write_k_3
        generic :: write_k => m_write_k_1
        generic :: write_k => m_write_k_2
        generic :: write_k => m_write_k_3
        
        procedure, private, nopass :: integrate_susc => m_integrate_susc_2_mom
        procedure, private, nopass :: integrate_susc_freq => m_integrate_susc_2_freq
        
        procedure, private, nopass :: m_write_k_tau_2
        procedure, private, nopass :: m_write_k_tau_4
        generic :: write_k_tau => m_write_k_tau_2
        generic :: write_k_tau => m_write_k_tau_4
        
        procedure, private, nopass :: write_r_tau => m_write_r_tau_2
        
        procedure, private, nopass :: write_w => m_write_w
        
        procedure, private :: write_obs_equal => m_write_obs_equal
        procedure, private :: write_obs_tau => m_write_obs_tau
        
        procedure, public :: preq => m_process_obs_equal
        procedure, public :: prtau => m_process_obs_tau
    end type FourierTrans

contains
    subroutine m_write_real_1(gr, filek) ! overloading routine for other correlations
! Arguments:
        real(kind=8), dimension(Lq), intent(in) :: gr
        character(len=*), intent(in) :: filek
! Local: 
        integer :: nr
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            write(20,*) gr(nr)
        enddo
        close(20)
        return
    end subroutine m_write_real_1
    
    subroutine m_write_real_2(gr, filek) ! overloading routine
! Arguments:
        real(kind=8), dimension(:,:), intent(in) :: gr
        character(len=*), intent(in) :: filek
! Local: 
        integer :: nr, nf
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            if (size(gr, 1) == Lq) then ! (Lq, Nbond)
                do nf = 1, Nbond
                    write(20,*) gr(nr, nf)
                enddo
            elseif (size(gr, 2) == Lq) then ! (Naux, Lq)
                do nf = 1, Naux
                    write(20,*) gr(nf, nr)
                enddo
            else
                write(6,*) "ERROR: incorrect input size in write_real_2"; stop
            endif
        enddo
        close(20)
        return
    end subroutine m_write_real_2
    
    subroutine m_write_cmplx_1(gr, filek)
        complex(kind=8), dimension(Lq), intent(in) :: gr
        character(len=*), intent(in) :: filek
        integer :: nr, no1, no2
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            write(20,*) gr(nr)
        enddo
        close(20)
        return
    end subroutine m_write_cmplx_1
    
    subroutine m_write_cmplx_3(gr, filek)
        complex(kind=8), dimension(Lq, Norb, Norb), intent(in) :: gr
        character(len=*), intent(in) :: filek
        integer :: nr, no1, no2
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do no2 = 1, Norb
            do no1 = 1, Norb
                do nr = 1, Lq
                    write(20, '(1X,E16.8)', advance='no') real(gr(nr, no1, no2))
                enddo
            enddo
        enddo
        do no2 = 1, Norb
            do no1 = 1, Norb
                do nr = 1, Lq
                    write(20, '(1X,E16.8)', advance='no') aimag(gr(nr, no1, no2))
                enddo
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_cmplx_3
   
    subroutine m_write_reciprocal_1(gk, filek)
        complex(kind=8), dimension(Lq), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            write(20,*) gk(nk)
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_1
    
    subroutine m_write_reciprocal_2(gk, filek)
        complex(kind=8), dimension(Lq, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk, no
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            do no = 1, Nbond
                write(20,*) gk(nk, no)
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_2
    
    subroutine m_write_reciprocal_3(gk, filek)
        complex(kind=8), dimension(Lq, Nbond, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk, no1, no2
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            do no2 = 1, Nbond
                do no1 = 1, Nbond
                    write(20,*) gk(nk, no1, no2)
                enddo
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_3

    subroutine m_write_k_1(gk, filek, momindex)
        complex(kind=8), dimension(Lq), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30, *) real(gk(momindex)), imag(gk(momindex))
        close(30)
        return
    end subroutine m_write_k_1
    
    subroutine m_write_k_2(gk, filek, momindex, nf)
        complex(kind=8), dimension(Lq, Norb), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, nf
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30,*) real(gk(momindex, nf)), imag(gk(momindex, nf))
        close(30)
        return
    end subroutine m_write_k_2
    
    subroutine m_write_k_3(gk, filek, momindex, no1, no2)
        complex(kind=8), dimension(Lq, Norb, Norb), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, no1, no2
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30,*) real(gk(momindex, no1, no2)), imag(gk(momindex, no1, no2))
        close(30)
        return
    end subroutine m_write_k_3

    subroutine m_integrate_susc_2_mom(gr, gk)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        complex(kind=8), dimension(Lq), intent(out) :: gk
        integer :: nt, nk
        gk = dcmplx(0.d0, 0.d0)
        do nt = 1, Ltrot
            do nk = 1, Lq
                gk(nk) = gk(nk) + gr(nk, nt)
            enddo
        enddo
        gk = gk * dcmplx(Dtau, 0.d0)
        return
    end subroutine m_integrate_susc_2_mom
    
    subroutine m_integrate_susc_2_freq(gr, gk)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        complex(kind=8), dimension(Ltrot), intent(out) :: gk
        integer :: nt, nw
        gk = dcmplx(0.d0, 0.d0)
        do nw = 1, Ltrot
            do nt = 1, Ltrot
                gk(nw) = gk(nw) + exp( dcmplx(0.d0, 2.d0*Pi*dble((nt-1)*(nw-1))/dble(Ltrot)) ) * gr(1, nt)
            enddo
        enddo
        gk = gk * dcmplx(Dtau, 0.d0)
        return
    end subroutine m_integrate_susc_2_freq
    
    subroutine m_write_r_tau_2(gr, filek, momindex)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        integer :: nt
        complex(kind=8) :: tmp
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%aimj_v(momindex, 1), Latt%aimj_v(momindex, 2)
        do nt = 1, Ltrot
            tmp = gr(momindex, nt) / dble(Lq)
            write(20,*) real(tmp), imag(tmp)
        enddo
        close(20)
        return
    end subroutine m_write_r_tau_2
    
    subroutine m_write_w(gk, filek)
        complex(kind=8), dimension(Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        real(kind=8) :: omega
        integer :: nw
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nw = 1, Ltrot
            omega = 2.d0 * Pi * dble(nw-1) / beta ! Matsubara frequency
            write(20,*) omega
            write(20,*) real(gk(nw)), imag(gk(nw))
        enddo
        close(20)
        return
    end subroutine m_write_w
    
    subroutine m_write_k_tau_2(gk, filek, momindex)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        integer :: nt
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%xk_v(momindex, 1), Latt%xk_v(momindex, 2)
        do nt = 1, Ltrot
            write(20,*) real(gk(momindex, nt)), imag(gk(momindex, nt))
        enddo
        close(20)
        return
    end subroutine m_write_k_tau_2
    
    subroutine m_write_k_tau_4(gk, filek, momindex, no1, no2)
        complex(kind=8), dimension(Lq, Nbond, Nbond, Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, no1, no2
        integer :: nt
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%xk_v(momindex, 1), Latt%xk_v(momindex, 2)
        do nt = 1, Ltrot
            write(20,*) real(gk(momindex, no1, no2, nt))
        enddo
        close(20)
        return
    end subroutine m_write_k_tau_4
    
    subroutine m_write_obs_equal(this, Obs)
        class(FourierTrans), intent(inout) :: this
        class(ObserEqual), intent(in) :: Obs
        complex(kind=8) :: correlation_up(Lq, Norb, Norb), correlation_do(Lq, Norb, Norb), correlation_updo(Lq)
        complex(kind=8) :: dentot_corr(Lq), dentot_corr_up(Lq), dentot_corr_do(Lq), SF_corr(Lq)
        complex(kind=8) :: SF_structure(Lq), SF_structure_up(Lq), SF_structure_do(Lq),PF_structure(Lq), C3_structure(Lq), dentot_structure(Lq), dentot_structure_up(Lq), dentot_structure_do(Lq)
        character(len=25) :: filek
        integer :: indexzero, no1, no2, i
        real(kind=8) temp
        
        indexzero = Latt%inv_cell_list(1, 1)
            
        open(unit=80, file='density_up', status='unknown', action="write", position="append")
        write(80,*) Obs%density_up
        close(80)
            
        open(unit=80, file='density_do', status='unknown', action="write", position="append")
        write(80,*) Obs%density_do
        close(80)
            
        open(unit=80, file='kinetic', status='unknown', action="write", position="append")
        write(80,*) Obs%kinetic
        close(80)
            
        open(unit=80, file='doubleOcc', status='unknown', action="write", position="append")
        write(80,*) Obs%doubleOcc
        close(80)
            
        open(unit=80, file='squareOcc', status='unknown', action="write", position="append")
        write(80,*) Obs%squareOcc
        close(80)
            
        open(unit=80, file='num_up', status='unknown', action="write", position="append")
        write(80,*) Obs%num_up
        close(80)
            
        open(unit=80, file='num_do', status='unknown', action="write", position="append")
        write(80,*) Obs%num_do
        close(80)
            
        open(unit=80, file='numsquare_up', status='unknown', action="write", position="append")
        write(80,*) Obs%numsquare_up
        close(80)
            
        open(unit=80, file='numsquare_do', status='unknown', action="write", position="append")
        write(80,*) Obs%numsquare_do
        close(80)
            
        open(unit=80, file='C3breaking', status='unknown', action="write", position="append")
        write(80,*) Obs%C3breaking
        close(80)

        open(unit=80, file='C3breaking_up', status='unknown', action="write", position="append")
        write(80,*) Obs%C3breaking_up
        close(80)

        open(unit=80, file='C3breaking_do', status='unknown', action="write", position="append")
        write(80,*) Obs%C3breaking_do
        close(80)

        dentot_corr_up = dcmplx(0.d0, 0.d0)
        dentot_corr_do = dcmplx(0.d0, 0.d0)
        do no1 = 1, Norb
            do no2 = 1, Norb
                dentot_corr_up = dentot_corr_up + Obs%den_corr_up(:, no1, no2)
                dentot_corr_do = dentot_corr_do + Obs%den_corr_do(:, no1, no2)
            enddo
        enddo

        dentot_corr = dentot_corr_up + dentot_corr_do + Obs%den_corr_updo + Obs%den_corr_updo
        SF_corr = Obs%SF_corr_up + Obs%SF_corr_do

        call Fourier_R_to_K(Obs%den_corr_up, correlation_up, Latt)
        call Fourier_R_to_K(Obs%den_corr_do, correlation_do, Latt)
        call Fourier_R_to_K(Obs%den_corr_updo, correlation_updo, Latt)

        C3_structure = dcmplx(0.d0, 0.d0)

        call Fourier_R_to_K(Obs%SF_corr_up, SF_structure_up, Latt)
        call Fourier_R_to_K(Obs%SF_corr_do, SF_structure_do, Latt)
        call Fourier_R_to_K(Obs%PF_corr, PF_structure, Latt)
        call Fourier_R_to_K(Obs%C3_corr, C3_structure, Latt)
        call Fourier_R_to_K(dentot_corr, dentot_structure, Latt)
        call Fourier_R_to_K(dentot_corr_up, dentot_structure_up, Latt)
        call Fourier_R_to_K(dentot_corr_do, dentot_structure_do, Latt)
        SF_structure = SF_structure_up + SF_structure_do

        do no1 = 1, Norb
            do no2 = 1, Norb
                write(filek, "('den_upup_sub',I0,'',I0)") no1, no2
                call this%write_k(correlation_up, filek, indexzero, no1, no2 )
                write(filek, "('den_dodo_sub',I0,'',I0)") no1, no2
                call this%write_k(correlation_do, filek, indexzero, no1, no2 )
            enddo
        enddo

        filek = 'den_updo'
        call this%write_k(correlation_updo, filek, indexzero )

        filek = 'SF_Gamma'
        call this%write_k(SF_structure, filek, indexzero )
        filek = 'SF_up_Gamma'
        call this%write_k(SF_structure_up, filek, indexzero )
        filek = 'SF_do_Gamma'
        call this%write_k(SF_structure_do, filek, indexzero )
        filek = 'PF_Gamma'
        call this%write_k(PF_structure, filek, indexzero )
        filek = 'C3_Gamma'
        call this%write_k(C3_structure, filek, indexzero )
        filek = 'dentot_Gamma'
        call this%write_k(dentot_structure, filek, indexzero )
        filek = 'dentot_up_Gamma'
        call this%write_k(dentot_structure_up, filek, indexzero )
        filek = 'dentot_do_Gamma'
        call this%write_k(dentot_structure_do, filek, indexzero )

        filek = 'SF_corr'
        call this%write_cmplx(SF_corr, filek)
        filek = 'SF_corr_up'
        call this%write_cmplx(Obs%SF_corr_up, filek)
        filek = 'SF_corr_do'
        call this%write_cmplx(Obs%SF_corr_do, filek)
        filek = 'PF_corr'
        call this%write_cmplx(Obs%PF_corr, filek)
        filek = 'C3_corr'
        call this%write_cmplx(Obs%C3_corr, filek)
        filek = 'dentot_corr'
        call this%write_cmplx(dentot_corr, filek)
        filek = 'dentot_corr_up'
        call this%write_cmplx(dentot_corr_up, filek)
        filek = 'dentot_corr_do'
        call this%write_cmplx(dentot_corr_do, filek)

        filek = 'SF_structure_up'
        call this%write_reciprocal(SF_structure_up, filek)
        filek = 'SF_structure_do'
        call this%write_reciprocal(SF_structure_do, filek)
        filek = 'PF_structure'
        call this%write_reciprocal(PF_structure, filek)
        filek = 'C3_structure'
        call this%write_reciprocal(C3_structure, filek)
        filek = 'dentot_structure'
        call this%write_reciprocal(dentot_structure, filek)
        filek = 'dentot_structure_up'
        call this%write_reciprocal(dentot_structure_up, filek)
        filek = 'dentot_structure_do'
        call this%write_reciprocal(dentot_structure_do, filek)

        ! filek = 'green'
        ! call this%write_cmplx(Obs%single_corr, filek)

        return
    end subroutine m_write_obs_equal

    subroutine m_process_obs_equal(this, Obs)
!#define DEC
        include 'mpif.h'
! Arguments:
        class(FourierTrans), intent(inout) :: this
        class(ObserEqual), intent(inout) :: Obs
! Local:
!        complex(kind=8), dimension(Lq, Nbond, Nbond) :: Collect3
        real(kind=8), dimension(Lq) :: Collect1
        real(kind=8), dimension(Lq, Nbond) :: Collect2
        real(kind=8), dimension(Naux, Lq) :: Collect2prime
        real(kind=8) :: Collect0, Collect1prime(Nbond)
        complex(kind=8), dimension(Lq) :: Collect1cmplx
        complex(kind=8) :: Collect3(Lq, Norb, Norb)
        integer :: N
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density_do = Collect0/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%kinetic, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%kinetic = Collect0/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%doubleOcc, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%doubleOcc = Collect0/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%squareOcc, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%squareOcc = Collect0/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%num_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%num_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%num_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%num_do = Collect0/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%numsquare_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%numsquare_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%numsquare_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%numsquare_do = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%C3breaking, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%C3breaking = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%C3breaking_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%C3breaking_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%C3breaking_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%C3breaking_do = Collect0/dble(ISIZE)

        N = Lq * Norb * Norb

        Collect3 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%den_corr_up, Collect3, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_corr_up = Collect3/dcmplx(dble(ISIZE),0.d0)

        Collect3 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%den_corr_do, Collect3, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_corr_do = Collect3/dcmplx(dble(ISIZE),0.d0)

        Collect3 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%single_corr, Collect3, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%single_corr = Collect3/dcmplx(dble(ISIZE),0.d0)

        N = Lq

        Collect1cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%den_corr_updo, Collect1cmplx, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_corr_updo = Collect1cmplx/dcmplx(dble(ISIZE),0.d0)

        Collect1cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%SF_corr_up, Collect1cmplx, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%SF_corr_up = Collect1cmplx/dcmplx(dble(ISIZE),0.d0)

        Collect1cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%SF_corr_do, Collect1cmplx, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%SF_corr_do = Collect1cmplx/dcmplx(dble(ISIZE),0.d0)
        
        Collect1cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%PF_corr, Collect1cmplx, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%PF_corr = Collect1cmplx/dcmplx(dble(ISIZE),0.d0)

        Collect1cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%C3_corr, Collect1cmplx, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%C3_corr = Collect1cmplx/dcmplx(dble(ISIZE),0.d0)

        if (IRANK == 0) call this%write_obs_equal(Obs)
        return
    end subroutine m_process_obs_equal
    
    subroutine m_write_obs_tau(this, Obs)
        class(FourierTrans), intent(inout) :: this
        class(ObserTau), intent(in) :: Obs
        character(len=25) :: filek
        return
    end subroutine m_write_obs_tau
    
    subroutine m_process_obs_tau(this, Obs)
!#define DEC
        include 'mpif.h'
! Arguments:
        class(FourierTrans), intent(inout) :: this
        class(ObserTau), intent(inout) :: Obs
! Local:
        complex(kind=8), dimension(Lq, Ltrot) :: Collect2
!        complex(kind=8), dimension(Lq, Nbond, Nbond, Ltrot) :: Collect4
        integer :: N

        N = Lq * Ltrot

        if (IRANK == 0) call this%write_obs_tau(Obs)
        return
    end subroutine m_process_obs_tau
end module FourierTrans_mod