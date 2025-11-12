module Stabilize_mod
    use ProcessMatrix
    use MyMats
    implicit none
    
    private
    public :: Wrap_pre, Wrap_L, Wrap_R, Wrap_tau, Stabilize_init, Stabilize_clear, stab_green
    complex(kind=8) ::  Z_one
    type(PropGreen), allocatable :: Gr_tmp
    
contains
    subroutine Stabilize_init()
        Z_one = dcmplx(1.d0, 0.d0)
        allocate(Gr_tmp)
        call Gr_tmp%make()
        return
    end subroutine Stabilize_init
    
    subroutine Stabilize_clear()
        deallocate(Gr_tmp)
        return
    end subroutine Stabilize_clear
    
    ! QDRP_decompose function removed - no longer needed for PQMC vector operations
    
    subroutine stab_UR(Prop)
        ! PQMC version: normalize right vector B(tau,0)P using LAPACK
        class(Propagator), intent(inout) :: Prop
        real(kind=8) :: norm_val
        
        ! LAPACK function interfaces
        real(kind=8), external :: DZNRM2
        external :: ZDSCAL
        
        ! Calculate ||B(tau,0)P|| using LAPACK DZNRM2
        norm_val = DZNRM2(Ndim, Prop%UUR(1,1), 1)
        
        ! Normalize: P_R = B(tau,0)P / ||B(tau,0)P||
        if (norm_val > 1.d-14) then
            call ZDSCAL(Ndim, 1.d0/norm_val, Prop%UUR(1,1), 1)
        else
            write(6,*) "Warning: very small norm in stab_UR, norm=", norm_val
        endif
        
        return
    end subroutine stab_UR
    
    subroutine  stab_UL(Prop)
        ! PQMC version: normalize left vector P^dagger B(2theta,tau) using LAPACK
        class(Propagator), intent(inout) :: Prop
        real(kind=8) :: norm_val
        
        ! LAPACK function interfaces
        real(kind=8), external :: DZNRM2
        external :: ZDSCAL
        
        ! Calculate ||P^dagger B(2theta,tau)|| using LAPACK DZNRM2
        norm_val = DZNRM2(Ndim, Prop%UUL(1,1), 1)
        
        ! Normalize: P_L^dagger = P^dagger B(2theta,tau) / ||P^dagger B(2theta,tau)||
        if (norm_val > 1.d-14) then
            call ZDSCAL(Ndim, 1.d0/norm_val, Prop%UUL(1,1), 1)
        else
            write(6,*) "Warning: very small norm in stab_UL, norm=", norm_val
        endif
        
        return
    end subroutine stab_UL
    
    subroutine stab_green(Gbar, Prop, nt)
        ! PQMC version: calculate Gbar using optimized LAPACK operations
        ! Arguments: 
        complex(kind=8), dimension(Ndim, Ndim), intent(out) :: Gbar
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: nt
        ! Local: 
        complex(kind=8) :: denominator, alpha
        
        ! LAPACK function interfaces
        complex(kind=8), external :: ZDOTU
        real(kind=8), external :: DZNRM2
        external :: ZDSCAL, ZGERU
        
        ! Calculate denominator: P_L^dagger * P_R using LAPACK ZDOTU
        ! Since UUL is already P_L^dagger (1,Ndim) and UUR is P_R (Ndim,1):
        ! denominator = sum_i UUL(1,i) * UUR(i,1) (no conjugation needed)
        denominator = ZDOTU(Ndim, Prop%UUL(1,1), 1, Prop%UUR(1,1), 1)
        
        ! Check for numerical stability
        if (abs(denominator) < 1.d-14) then
            write(6,*) "Warning: small denominator in stab_green, |denom|=", abs(denominator)
        endif
        
        ! Calculate Gbar = N_b * (P_R * P_L^dagger) / denominator
        ! Using ZGERU: A := alpha*x*y^T + A (no conjugation)
        ! Here: Gbar := (N_b/denominator) * UUR * UUL + 0
        ! Since UUL is already P_L^dagger, no conjugate needed
        alpha = dcmplx(dble(Nbos), 0.d0) / denominator
        Gbar = dcmplx(0.d0, 0.d0)
        call ZGERU(Ndim, Ndim, alpha, Prop%UUR(1,1), 1, Prop%UUL(1,1), 1, Gbar, Ndim)
        
        return
    end subroutine stab_green
        
    subroutine stab_green_big(Prop)
        ! PQMC version: simplified version for big Green's function
        ! Arguments:
        class(Propagator), intent(in) :: Prop
        
        ! For PQMC, this function may not be needed in the same way
        ! as the original matrix-based algorithm. 
        ! Placeholder implementation - may need future adaptation
        write(6,*) "stab_green_big: PQMC version not yet implemented"
        
        return
    end subroutine stab_green_big
    
    real(kind=8) function compare_mat(Gr, Gr2) result(dif)
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: Gr, Gr2
        integer :: nl, nr
        dif = 0.d0
        do nr = 1, Ndim
            do nl = 1, Ndim
                dif = dif + real( abs(Gr(nl, nr) - Gr2(nl, nr)) )
            enddo
        enddo
        return
    end function compare_mat
    
    subroutine Wrap_pre(Prop, WrList, nt)
        ! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        ! Local: 
        integer :: nt_st
        ! Allow wrapping at nt=Ltrot even if not divisible by Nwrap
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. Ltrot) then
            write(6,*) "incorrect preortho time slice, NT = ", nt; stop
        endif
        ! Calculate correct nt_st for non-divisible Ltrot/Nwrap cases
        if (nt == 0) then
            nt_st = 0
        else if (mod(nt, Nwrap) == 0) then
            nt_st = nt / Nwrap
        else
            nt_st = nt / Nwrap + 1  ! Next interval for partial wrap
        endif
        if (nt .ne. 0) call stab_UR(Prop)
        WrList%URlist(1:Ndim, 1, nt_st) = Prop%UUR(1:Ndim, 1)
        if (nt == Ltrot) then
            Prop%Gbar = dcmplx(0.d0, 0.d0)
            call stab_green(Prop%Gbar, Prop, nt)
        endif
        return
    end subroutine Wrap_pre
    
    subroutine Wrap_L(Prop, WrList, nt, flag)
        ! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        character(len=*), optional, intent(in) :: flag
        ! Local: 
        integer :: nt_st
        real(kind=8) :: dif
        ! Allow wrapping at nt=Ltrot even if not divisible by Nwrap
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0 .and. nt .ne. Ltrot) then
            write(6,*) "incorrect ortholeft time slice, NT = ", nt; stop
        endif
        ! Calculate correct nt_st for non-divisible Ltrot/Nwrap cases
        if (nt == 0) then
            nt_st = 0
        else if (mod(nt, Nwrap) == 0) then
            nt_st = nt / Nwrap
        else
            nt_st = nt / Nwrap + 1  ! Next interval for partial wrap
        endif
        ! Bounds check for WrapList array access
        if (nt_st > Nst) then
            write(6,*) "Warning: nt_st exceeds Nst in Wrap_L, nt=", nt, "nt_st=", nt_st
            nt_st = Nst  ! Use last valid index
        endif
        Prop%UUR(1:Ndim, 1) = WrList%URlist(1:Ndim, 1, nt_st)
        if (nt == 0) then ! clear URlist
            WrList%URlist = dcmplx(0.d0, 0.d0)
        endif
        if (nt .ne. Ltrot) then
            call stab_UL(Prop)
            Gr_tmp%Gr00 = dcmplx(0.d0, 0.d0)
            call stab_green(Gr_tmp%Gr00, Prop, nt)
            dif = compare_mat(Gr_tmp%Gr00, Prop%Gbar)
            if (dif > Prop%Xmaxm) Prop%Xmaxm = dif
            if (dif .ge. 5.5d-5) write(6,*) nt, dif, "left ortho unstable in RANK ", IRANK
            if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
            Prop%Gbar = Gr_tmp%Gr00
        endif
        WrList%ULlist(1, 1:Ndim, nt_st) = Prop%UUL(1, 1:Ndim)
        return
    end subroutine Wrap_L
    
    subroutine Wrap_R(Prop, WrList, nt, flag)
        ! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        character(len=*), optional, intent(in) :: flag
        ! Local: 
        integer :: nt_st
        real(kind=8) :: dif
        ! Allow wrapping at nt=Ltrot even if not divisible by Nwrap
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0 .and. nt .ne. Ltrot) then
            write(6,*) "incorrect orthoright time slice, NT = ", nt; stop
        endif
        ! Calculate correct nt_st for non-divisible Ltrot/Nwrap cases
        if (nt == 0) then
            nt_st = 0
        else if (mod(nt, Nwrap) == 0) then
            nt_st = nt / Nwrap
        else
            nt_st = nt / Nwrap + 1  ! Next interval for partial wrap
        endif
        ! Bounds check for WrapList array access
        if (nt_st > Nst) then
            write(6,*) "Warning: nt_st exceeds Nst in Wrap_R, nt=", nt, "nt_st=", nt_st
            nt_st = Nst  ! Use last valid index
        endif
        Prop%UUL(1, 1:Ndim) = WrList%ULlist(1, 1:Ndim, nt_st)
        if (nt == Ltrot) then
            WrList%ULlist = dcmplx(0.d0, 0.d0)
        endif
        if (nt .ne. 0) then
            call stab_UR(Prop)
            Gr_tmp%Gr00 = dcmplx(0.d0, 0.d0)
            call stab_green(Gr_tmp%Gr00, Prop, nt)
            dif = compare_mat(Gr_tmp%Gr00, Prop%Gbar)
            if (dif > Prop%Xmaxm) Prop%Xmaxm = dif
            if (dif .ge. 5.5d-5) write(6,*) nt, dif, "right ortho unstable in RANK ", IRANK
            if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
            Prop%Gbar = Gr_tmp%Gr00
        endif
        WrList%URlist(1:Ndim, 1, nt_st) = Prop%UUR(1:Ndim, 1)
        return
    end subroutine Wrap_R

    subroutine Wrap_tau(Prop, PropGr, WrList, nt)
        ! PQMC version: simplified time-dependent wrapping
        ! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: PropGr
        class(WrapList), intent(in) :: WrList
        integer, intent(in) :: nt
        ! Local: 
        integer :: nt_st
        
        ! Allow wrapping at nt=Ltrot even if not divisible by Nwrap
        if ((mod(nt, Nwrap) .ne. 0 .and. nt .ne. Ltrot) .or. nt == 0) then
            write(6,*) "incorrect orthobig time slice, NT = ", nt; stop
        endif
        ! Calculate correct nt_st for non-divisible Ltrot/Nwrap cases
        if (mod(nt, Nwrap) == 0) then
            nt_st = nt / Nwrap
        else
            nt_st = nt / Nwrap + 1  ! Next interval for partial wrap
        endif
        Prop%UUL(1, 1:Ndim) = WrList%ULlist(1, 1:Ndim, nt_st)
        call stab_UR(Prop)
        
        ! For PQMC, the time-dependent Green's function calculation
        ! needs to be adapted. For now, keep minimal functionality.
        write(6,*) "Wrap_tau: PQMC time-dependent Green's function not fully implemented"
        
        return
    end subroutine Wrap_tau
end module Stabilize_mod