module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserEqual
        complex(kind=8), dimension(:,:,:), allocatable  :: den_corr_up, den_corr_do, single_corr
        complex(kind=8), dimension(:), allocatable      :: den_corr_updo
        complex(kind=8), dimension(:), allocatable      :: SF_corr_up, SF_corr_do
        complex(kind=8), dimension(:), allocatable      :: PF_corr
        complex(kind=8), dimension(:), allocatable      :: C3_corr, C3_corr_up
        real(kind=8)                                    :: density_up,  density_do
        real(kind=8)                                    :: kinetic, doubleOcc, squareOcc, nearestOcc, IPR
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
        allocate( this%SF_corr_up(Lq), this%SF_corr_do(Lq), this%PF_corr(Lq), this%C3_corr(Lq), this%C3_corr_up(Lq) )
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        deallocate( this%den_corr_up, this%den_corr_do, this%single_corr, this%den_corr_updo )
        deallocate( this%SF_corr_up, this%SF_corr_do, this%PF_corr, this%C3_corr, this%C3_corr_up )
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%den_corr_up   = dcmplx(0.d0,0.d0)
        this%den_corr_do   = dcmplx(0.d0,0.d0)
        this%single_corr   = dcmplx(0.d0,0.d0)
        this%den_corr_updo = dcmplx(0.d0,0.d0)
        this%SF_corr_up   = dcmplx(0.d0,0.d0)
        this%SF_corr_do   = dcmplx(0.d0,0.d0)
        this%PF_corr      = dcmplx(0.d0,0.d0)
        this%C3_corr     = dcmplx(0.d0,0.d0)
        this%C3_corr_up  = dcmplx(0.d0,0.d0)
        this%density_up  = 0.d0
        this%density_do  = 0.d0
        this%kinetic     = 0.d0
        this%doubleOcc   = 0.d0
        this%squareOcc   = 0.d0
        this%IPR         = 0.d0
        this%nearestOcc  = 0.d0
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
        this%SF_corr_up  = this%SF_corr_up  * znorm
        this%SF_corr_do  = this%SF_corr_do  * znorm
        this%PF_corr     = this%PF_corr     * znorm
        this%C3_corr     = this%C3_corr     * znorm
        this%C3_corr_up  = this%C3_corr_up  * znorm
        this%density_up  = this%density_up  * znorm
        this%density_do  = this%density_do  * znorm
        this%kinetic     = this%kinetic     * znorm
        this%doubleOcc   = this%doubleOcc   * znorm
        this%squareOcc   = this%squareOcc   * znorm
        this%IPR         = this%IPR         * znorm
        this%nearestOcc  = this%nearestOcc  * znorm
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
        complex(kind=8), dimension(Ndim, Ndim) :: Gbar
        complex(kind=8) :: alpha, overlap
        integer :: i, j, no1, no2, ii, jj, imj, nb, no
        real(kind=8) :: GGfactor, num_per_site
        complex(kind=8) :: den_temp(Lq, Norb, Norb), den_temp_up(Lq, Norb, Norb), den_temp_do(Lq, Norb, Norb)
        complex(kind=8) :: temp_upup, temp_dodo, temp_updo, temp_doup, temp_up, temp_do
        complex(kind=8), external :: ZDOTU
        external :: ZGERU

        ! Reconstruct rank-1 Gbar = (N_b / overlap) * (UUR * UUL)
        overlap = ZDOTU(Ndim, Prop%UUL(1,1), 1, Prop%UUR(1,1), 1)
        alpha = dcmplx(dble(Nbos), 0.d0) / overlap
        Gbar = dcmplx(0.d0, 0.d0)
        call ZGERU(Ndim, Ndim, alpha, Prop%UUR(1,1), 1, Prop%UUL(1,1), 1, Gbar, Ndim)

        Grup    = Gbar + ZKRON                 !   Gr(i, j)    = <b_i b^+_j > = Gbar + I  
        Grupc   = transpose(Gbar)              !   Grc(i, j)   = <b^+_i b_j > = transpose(Gbar)
        
        Grdo    = dconjg(Grup)                       !   Gr(i, j)    = <c_i c^+_j > = conjg(Gbar + I)
        Grdoc   = dconjg(Grupc)                      !   Grc(i, j)   = <c^+_i c_j > = conjg(transpose(Gbar))

        GGfactor = dble(Nbos-1) / dble(Nbos)
        num_per_site = dble(Nbos) / dble(Norb*Lq)
        den_temp = dcmplx(0.d0, 0.d0)
        den_temp_up = dcmplx(0.d0, 0.d0)
        den_temp_do = dcmplx(0.d0, 0.d0)
        
        do ii = 1, Ndim
            this%density_up = this%density_up + real( Grupc(ii,ii) ) / dble(Lq)
            this%density_do = this%density_do + real( Grdoc(ii,ii) ) / dble(Lq)
            this%doubleOcc  = this%doubleOcc  + real( Grupc(ii,ii) * Grdoc(ii,ii) ) / dble(Lq)
            this%squareOcc  = this%squareOcc  + ( GGfactor * real( Grupc(ii,ii) * Grupc(ii,ii) + Grdoc(ii,ii) * Grdoc(ii,ii) ) + real( Grupc(ii,ii) + Grdoc(ii,ii) ) ) / dble(Lq)
            this%IPR        = this%IPR        + real( Grupc(ii,ii) * Grupc(ii,ii) ) / dble(Nbos**2)
            this%num_up = this%num_up + real( Grupc(ii,ii) ) 
            this%num_do = this%num_do + real( Grdoc(ii,ii) ) 
        enddo

        do ii = 1, Ndim
            do jj = 1, Ndim
                this%numsquare_up = this%numsquare_up + GGfactor * real( Grupc(ii,ii) * Grupc(jj,jj) )
                this%numsquare_do = this%numsquare_do + GGfactor * real( Grdoc(ii,ii) * Grdoc(jj,jj) )
                if (ii==jj) then
                    this%numsquare_up = this%numsquare_up + real( Grupc(ii,ii) )
                    this%numsquare_do = this%numsquare_do + real( Grdoc(ii,ii) )
                endif
            enddo
        enddo

        do i = 1, Lq
            do j = 1, Lq
                imj = Latt%imj(i, j)
                do no = 1, Norb
                    ii = Latt%inv_dim_list(i, no)
                    jj = Latt%inv_dim_list(j, no)
                    this%SF_corr_up(imj) = this%SF_corr_up(imj) + Grupc(ii,jj) / dcmplx(dble(Norb*Lq), 0.d0)
                    this%SF_corr_do(imj) = this%SF_corr_do(imj) + Grdoc(ii,jj) / dcmplx(dble(Norb*Lq), 0.d0)
                    this%PF_corr(imj)    = this%PF_corr(imj)    + ( Grupc(ii,jj) * Grdoc(ii,jj) ) / dcmplx(dble(Norb*Lq), 0.d0)
                enddo
            enddo
        enddo

        do i = 1, Lq
            do j = 1, Lq
                imj = Latt%imj(i, j)
                do no1 = 1, Norb
                    do no2 = 1, Norb
                        ii = Latt%inv_dim_list(i, no1)
                        jj = Latt%inv_dim_list(j, no2)
                        temp_upup = ( GGfactor * Grupc(ii,ii) * Grupc(jj,jj) - num_per_site * (Grupc(ii,ii)+Grupc(jj,jj)-num_per_site) ) / dcmplx(dble(Lq), 0.d0)
                        temp_dodo = ( GGfactor * Grdoc(ii,ii) * Grdoc(jj,jj) - num_per_site * (Grdoc(ii,ii)+Grdoc(jj,jj)-num_per_site) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_up(imj, no1, no2) = this%den_corr_up(imj, no1, no2) + temp_upup
                        this%den_corr_do(imj, no1, no2) = this%den_corr_do(imj, no1, no2) + temp_dodo
                        den_temp_up(imj, no1, no2) = den_temp_up(imj, no1, no2) + temp_upup
                        den_temp_do(imj, no1, no2) = den_temp_do(imj, no1, no2) + temp_dodo
                        den_temp(imj, no1, no2) = den_temp(imj, no1, no2) + temp_upup
                        den_temp(imj, no1, no2) = den_temp(imj, no1, no2) + temp_dodo
                        if (ii==jj) then
                            temp_up = Grupc(ii,ii) / dcmplx(dble(Lq), 0.d0)
                            temp_do = Grdoc(ii,ii) / dcmplx(dble(Lq), 0.d0)
                            this%den_corr_up(imj, no1, no2) = this%den_corr_up(imj, no1, no2) + temp_up
                            this%den_corr_do(imj, no1, no2) = this%den_corr_do(imj, no1, no2) + temp_do
                            den_temp_up(imj, no1, no2) = den_temp_up(imj, no1, no2) + temp_up
                            den_temp_do(imj, no1, no2) = den_temp_do(imj, no1, no2) + temp_do
                            den_temp(imj, no1, no2) = den_temp(imj, no1, no2) + temp_up
                            den_temp(imj, no1, no2) = den_temp(imj, no1, no2) + temp_do
                        endif
                        temp_updo = ( Grupc(ii,ii) * Grdoc(jj,jj) - num_per_site * (Grupc(ii,ii)+Grdoc(jj,jj)-num_per_site) ) / dcmplx(dble(Lq), 0.d0)
                        temp_doup = ( Grdoc(ii,ii) * Grupc(jj,jj) - num_per_site * (Grdoc(ii,ii)+Grupc(jj,jj)-num_per_site) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_updo(imj) = this%den_corr_updo(imj) + temp_updo
                        den_temp(imj, no1, no2) = den_temp(imj, no1, no2) + temp_updo
                        den_temp(imj, no1, no2) = den_temp(imj, no1, no2) + temp_doup
                        this%single_corr(imj, no1, no2) = this%single_corr(imj, no1, no2) + Grupc(ii,jj) / dcmplx(dble(Lq), 0.d0)
                    enddo
                enddo
            enddo
        enddo

        ! do imj = 1, Lq
        !     this%C3_corr(imj) = this%C3_corr(imj) + ( &
        !         & 4*den_temp(imj,1,1) + 4*den_temp(imj,2,2) + 4*den_temp(imj,3,3) &
        !         & - 2*den_temp(imj,1,2) - 2*den_temp(imj,2,1) &
        !         & - 2*den_temp(imj,1,3) - 2*den_temp(imj,3,1) &
        !         & - 2*den_temp(imj,2,3) - 2*den_temp(imj,3,2) ) / dcmplx(6.d0,0.d0)
        !     this%C3_corr_up(imj) = this%C3_corr_up(imj) + ( &
        !         & 4*den_temp_up(imj,1,1) + 4*den_temp_up(imj,2,2) + 4*den_temp_up(imj,3,3) &
        !         & - 2*den_temp_up(imj,1,2) - 2*den_temp_up(imj,2,1) &
        !         & - 2*den_temp_up(imj,1,3) - 2*den_temp_up(imj,3,1) &
        !         & - 2*den_temp_up(imj,2,3) - 2*den_temp_up(imj,3,2) ) / dcmplx(6.d0,0.d0)
        ! enddo

        do ii = 1, Ndim
            do nb = 1, Nbond
                jj = Latt%L_bonds(ii, nb)
                this%kinetic = this%kinetic + RT * real( Grupc(ii,jj) + Grupc(jj,ii) + Grdoc(ii,jj) + Grdoc(jj,ii) ) / dble(Lq)
                this%nearestOcc = this%nearestOcc + ( GGfactor * real( Grupc(ii,ii) * Grupc(jj,jj) + Grdoc(ii,ii) * Grdoc(jj,jj) ) ) / dble(Lq)
            enddo
        enddo

        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod
