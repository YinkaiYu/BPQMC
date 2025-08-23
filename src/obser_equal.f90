module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserEqual
        complex(kind=8), dimension(:,:,:), allocatable  :: den_corr_up, den_corr_do, single_corr
        complex(kind=8), dimension(:), allocatable      :: den_corr_updo
        real(kind=8)                                    :: density_up,  density_do
        real(kind=8)                                    :: kinetic, doubleOcc, squareOcc
        real(kind=8)                                    :: num_up, num_do, numsquare_up, numsquare_do
    contains
        procedure :: make   => Obs_equal_make
        procedure :: reset  => Obs_equal_reset
        procedure :: ave    => Obs_equal_ave
        procedure :: calc   => Obs_equal_calc
        final :: Obs_equal_clear
    end type ObserEqual
    
contains
    subroutine Obs_equal_make(this)
        class(ObserEqual), intent(inout) :: this
        allocate( this%den_corr_up(Lq, Norb, Norb), this%den_corr_do(Lq, Norb, Norb), this%single_corr(Lq, Norb, Norb), this%den_corr_updo(Lq) )
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        deallocate( this%den_corr_up, this%den_corr_do, this%single_corr, this%den_corr_updo )
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%den_corr_up   = dcmplx(0.d0,0.d0)
        this%den_corr_do   = dcmplx(0.d0,0.d0)
        this%single_corr   = dcmplx(0.d0,0.d0)
        this%den_corr_updo = dcmplx(0.d0,0.d0)
        this%density_up  = 0.d0
        this%density_do  = 0.d0
        this%kinetic     = 0.d0
        this%doubleOcc   = 0.d0
        this%squareOcc   = 0.d0
        this%num_up      = 0.d0
        this%num_do      = 0.d0
        this%numsquare_up = 0.d0
        this%numsquare_do = 0.d0
        return
    end subroutine Obs_equal_reset
    
    subroutine Obs_equal_ave(this, Nobs)
        class(ObserEqual), intent(inout) :: this
        integer, intent(in) :: Nobs
        real(kind=8) :: znorm
        znorm = 1.d0 / dble(Nobs)
        this%den_corr_up = this%den_corr_up * znorm
        this%den_corr_do = this%den_corr_do * znorm
        this%single_corr = this%single_corr * znorm
        this%density_up  = this%density_up  * znorm
        this%density_do  = this%density_do  * znorm
        this%kinetic     = this%kinetic     * znorm
        this%doubleOcc   = this%doubleOcc   * znorm
        this%squareOcc   = this%squareOcc   * znorm
        this%den_corr_updo = this%den_corr_updo * znorm
        this%num_up      = this%num_up * znorm
        this%num_do      = this%num_do * znorm
        this%numsquare_up = this%numsquare_up * znorm
        this%numsquare_do = this%numsquare_do * znorm
        return
    end subroutine Obs_equal_ave
    
    subroutine Obs_equal_calc(this, Prop, ntau)
!   Arguments: 
        class(ObserEqual), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Grupc, Grup
        complex(kind=8), dimension(Ndim, Ndim) :: Grdoc, Grdo
        integer :: i, j, no1, no2, ii, jj, imj, nb, no
        
        Grup    = Prop%Gbar + ZKRON                 !   Gr(i, j)    = <b_i b^+_j > = Gbar + I  
        Grupc   = transpose(Prop%Gbar)              !   Grc(i, j)   = <b^+_i b_j > = transpose(Gbar)
        
        Grdo    = dconjg(Prop%Gbar + ZKRON)         !   Gr(i, j)    = <c_i c^+_j > = conjg(Gbar + I)
        Grdoc   = dconjg(transpose(Prop%Gbar))      !   Grc(i, j)   = <c^+_i c_j > = conjg(transpose(Gbar))
        
        do ii = 1, Ndim
            this%density_up = this%density_up + real( Grupc(ii,ii) ) / dble(Lq)
            this%density_do = this%density_do + real( Grdoc(ii,ii) ) / dble(Lq)
            this%doubleOcc  = this%doubleOcc  + real( Grupc(ii,ii) * Grdoc(ii,ii) ) / dble(Lq)
            this%squareOcc  = this%squareOcc  + real( Grupc(ii,ii) * Grupc(ii,ii) + Grdoc(ii,ii) * Grdoc(ii,ii) ) / dble(Lq)
            this%num_up = this%num_up + real( Grupc(ii,ii) ) 
            this%num_do = this%num_do + real( Grdoc(ii,ii) ) 
        enddo

        do ii = 1, Ndim
            do jj = 1, Ndim
                this%numsquare_up = this%numsquare_up + real( Grupc(ii,ii) * Grupc(jj,jj) + Grupc(ii,jj) * Grup(ii,jj) )
                this%numsquare_do = this%numsquare_do + real( Grdoc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,jj) * Grdo(ii,jj) )
            enddo
        enddo

        do i = 1, Lq
            do j = 1, Lq
                imj = Latt%imj(i, j)
                do no1 = 1, Norb
                    do no2 = 1, Norb
                        ii = Latt%inv_dim_list(i, no1)
                        jj = Latt%inv_dim_list(j, no2)
                        this%den_corr_up(imj, no1, no2) = this%den_corr_up(imj, no1, no2) + ( Grupc(ii,ii) * Grupc(jj,jj) + Grupc(ii,jj) * Grup(ii,jj) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_do(imj, no1, no2) = this%den_corr_do(imj, no1, no2) + ( Grdoc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,jj) * Grdo(ii,jj) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_updo(imj) = this%den_corr_updo(imj) + ( Grupc(ii,ii) * Grdoc(jj,jj) ) / dcmplx(dble(Lq), 0.d0)
                        this%single_corr(imj, no1, no2) = this%single_corr(imj, no1, no2) + Grupc(ii,jj) / dcmplx(dble(Lq), 0.d0)
                    enddo
                enddo
            enddo
        enddo

        do ii = 1, Ndim
            do nb = 1, Nbond
                jj = Latt%L_bonds(ii, nb)
                this%kinetic = this%kinetic + RT * real( Grupc(ii,jj) + Grupc(jj,ii) + Grdoc(ii,jj) + Grdoc(jj,ii) ) / dble(Lq)
            enddo
        enddo

        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod