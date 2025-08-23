module LocalU_mod
    use Multiply_mod
    use CalcBasic, only: Nbos
    implicit none
    
    public
    private :: LocalU_metro, LocalU_metro_therm
    
    
contains
    subroutine LocalU_init(Op_U)
        type(OperatorHubbard), intent(inout) :: Op_U
        call Op_U%Acc_U_local%init()
        call Op_U%Acc_U_therm%init()
        return
    end subroutine LocalU_init
    
    subroutine LocalU_clear(Op_U)
        type(OperatorHubbard), intent(inout) :: Op_U
        return
    end subroutine LocalU_clear
    
    subroutine LocalU_reset(Op_U)
        type(OperatorHubbard), intent(inout) :: Op_U
        call Op_U%Acc_U_local%reset()
        call Op_U%Acc_U_therm%reset()
        return
    end subroutine LocalU_reset
    
    subroutine LocalU_metro(Op_U, Gbar, iseed, nf, ii, ntau)
        use MyMats
! Arguments:
        type(OperatorHubbard), intent(inout) :: Op_U
	    complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gbar
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau, nf
!   Local: 
        real(kind=8), external :: ranf
        external :: ZSCAL
        real(kind=8) :: phi_old, phi_new
        real(kind=8) :: xflip, Xdif, random
        real(kind=8) :: ratio_abs
        complex(kind=8) :: ratio_Pfa, ratio_exp, r_b
        complex(kind=8) :: one_plus_delta

! Local update on space-time (ii, ntau) for auxiliary field flavor (nf)
        phi_old = Conf%phi_list(nf, ii, ntau)
        xflip = ranf(iseed)
        Xdif = dble((xflip - 0.5) * abs(shiftLoc))
        phi_new = phi_old + Xdif
! Calculate PQMC Metropolis ratio   
        call Op_U%get_delta(phi_old, phi_new)
        ratio_exp = Op_U%ratio_gaussian
        ! Calculate r_b = 1 + Δ_i/N_b * Gbar_{ii}
        r_b = dcmplx(1.d0,0.d0) + Op_U%Delta * Gbar(ii,ii) / dble(Nbos)
        ! Calculate ratio_Pfa = r_b^{N_b}
        ratio_Pfa = r_b ** Nbos
        ratio_abs = abs(ratio_exp * ratio_Pfa * dconjg(ratio_Pfa))
! Update Gbar and phi if accepted
        random = ranf(iseed)
        if (ratio_abs .gt. random) then
            call Op_U%Acc_U_local%count(.true.)

            ! PQMC Gbar update: Gbar' = (1+Δ)/r_b * Gbar
            ! Step 1: Apply (1+Δ_i) to row ii only
            one_plus_delta = dcmplx(1.d0,0.d0) + Op_U%Delta
            call ZSCAL(Ndim, one_plus_delta, Gbar(ii,1), Ndim)
            
            ! Step 2: Apply factor 1/r_b to entire matrix
            call ZSCAL(Ndim*Ndim, dcmplx(1.d0,0.d0)/r_b, Gbar, 1)

            Conf%phi_list(nf, ii, ntau) = phi_new
        else
            call Op_U%Acc_U_local%count(.false.)
        endif
        return
    end subroutine LocalU_metro
    
    subroutine LocalU_prop_L(Op_U, Prop, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = Ndim, 1, -1
            call LocalU_metro(Op_U, Prop%Gbar, iseed, nf, ii, ntau)
            call Op_U%mmult_L(Prop%Gbar, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_R(Prop%Gbar, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
        enddo
        ! wrap the left
        do ii = Ndim, 1, -1
            call Op_U%mmult_L(Prop%UUL, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
        enddo
        return
    end subroutine LocalU_prop_L
    
    subroutine LocalU_prop_R(Op_U, Prop, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%Gbar, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_L(Prop%Gbar, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
            call LocalU_metro(Op_U, Prop%Gbar, iseed, nf, ii, ntau)
        enddo
        ! wrap the right
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
        enddo
        return
    end subroutine LocalU_prop_R
    
    subroutine LocalU_metro_therm(Op_U, nf, ii, ntau, iseed)
! Arguments: 
        type(OperatorHubbard), intent(inout) :: Op_U
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau, nf
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: phi_old, phi_new
        real(kind=8) :: xflip, Xdif, random
        real(kind=8) :: ratio_abs
        complex(kind=8) :: ratio_exp

! Local update on space-time (ii, ntau) for auxiliary field flavor (nf)
        phi_old = Conf%phi_list(nf, ii, ntau)
        xflip = ranf(iseed)
        Xdif = dble((xflip - 0.5) * abs(shiftLoc))
        phi_new = phi_old + Xdif
! Calculate auxiliary Gaussian ratio   
        call Op_U%get_delta(phi_old, phi_new)
        ratio_exp = Op_U%ratio_gaussian
        ratio_abs = abs(ratio_exp)
! Upgrade phi
        random = ranf(iseed)
        if (ratio_abs .gt. random) then
            call Op_U%Acc_U_therm%count(.true.)
            Conf%phi_list(nf, ii, ntau) = phi_new
        else
            call Op_U%Acc_U_therm%count(.false.)
        endif
        return
    end subroutine LocalU_metro_therm
    
    subroutine LocalU_prop_therm(Op_U, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = 1, Ndim
            call LocalU_metro_therm(Op_U, nf, ii, ntau, iseed)
        enddo
        return
    end subroutine LocalU_prop_therm
end module LocalU_mod