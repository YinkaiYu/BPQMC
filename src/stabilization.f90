module Stabilize_mod
    use ProcessMatrix
    use MyMats
    use CalcBasic, only: norm_diff_vec, norm_threshold
    implicit none
    
    private
    public :: Wrap_pre, Wrap_L, Wrap_R, Wrap_tau
    public :: Stabilize_init, Stabilize_clear, stab_green
    complex(kind=8) ::  Z_one
    type(PropGreen), allocatable :: Gr_tmp
    logical :: warned_wrap_tau
    
contains
    subroutine Stabilize_init()
        Z_one = dcmplx(1.d0, 0.d0)
        warned_wrap_tau = .false.
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
            Prop%overlap = Prop%overlap / norm_val
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
            Prop%overlap = Prop%overlap / norm_val
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
        if (Nst <= 0) return
        if (nt < 1 .or. nt > Ltrot) then
            write(6,*) "incorrect preortho time slice, NT = ", nt; stop
        endif
        nt_st = (nt - 1) / Nwrap + 1  ! ceil(nt / Nwrap)
        call stab_UR(Prop)
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
        real(kind=8) :: dif, dif_wr
        if (Nst <= 0) return
        if (nt < 1 .or. nt > Ltrot) then
            write(6,*) "incorrect ortholeft time slice, NT = ", nt; stop
        endif
        nt_st = (nt - 1) / Nwrap + 1
        call stab_UR(Prop)
        dif_wr = norm_diff_vec(Prop%UUR(:,1), WrList%URlist(:,1,nt_st), Ndim)
        if (dif_wr > norm_threshold) write(6,*) "wrap_L UR diff at nt=", nt, " diff=", dif_wr, " rank=", IRANK
        Prop%UUR(1:Ndim, 1) = WrList%URlist(1:Ndim, 1, nt_st)
        call stab_UL(Prop)
        if (nt .ne. Ltrot) then
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
        real(kind=8) :: dif, dif_wr
        if (Nst <= 0) return
        if (nt < 1 .or. nt > Ltrot) then
            write(6,*) "incorrect orthoright time slice, NT = ", nt; stop
        endif
        nt_st = (nt - 1) / Nwrap + 1
        call stab_UL(Prop)
        dif_wr = norm_diff_vec(Prop%UUL(1,:), WrList%ULlist(1,1:Ndim,nt_st), Ndim)
        if (dif_wr > norm_threshold) write(6,*) "wrap_R UL diff at nt=", nt, " diff=", dif_wr, " rank=", IRANK
        Prop%UUL(1, 1:Ndim) = WrList%ULlist(1, 1:Ndim, nt_st)
        call stab_UR(Prop)
        Gr_tmp%Gr00 = dcmplx(0.d0, 0.d0)
        call stab_green(Gr_tmp%Gr00, Prop, nt)
        dif = compare_mat(Gr_tmp%Gr00, Prop%Gbar)
        if (dif > Prop%Xmaxm) Prop%Xmaxm = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "right ortho unstable in RANK ", IRANK
        if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
        Prop%Gbar = Gr_tmp%Gr00
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
        
        if (Nst <= 0) return
        if (nt <= 0 .or. nt > Ltrot) then
            write(6,*) "incorrect orthobig time slice, NT = ", nt; stop
        endif
        nt_st = (nt - 1) / Nwrap + 1
        Prop%UUL(1, 1:Ndim) = WrList%ULlist(1, 1:Ndim, nt_st)
        call stab_UR(Prop)
        
        if (.not. warned_wrap_tau) then
            write(6,*) "Wrap_tau: PQMC time-dependent Green's function not fully implemented"
            warned_wrap_tau = .true.
        endif
        
        return
    end subroutine Wrap_tau
end module Stabilize_mod
