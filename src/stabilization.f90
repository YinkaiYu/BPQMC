module Stabilize_mod
    use ProcessMatrix
    use MyMats
    use CalcBasic, only: norm_diff_vec, norm_threshold
    implicit none
    
    private
    public :: Wrap_pre, Wrap_L, Wrap_R, Wrap_tau
    public :: Stabilize_init, Stabilize_clear
    logical :: warned_wrap_tau
    
contains
    subroutine Stabilize_init()
        warned_wrap_tau = .false.
        return
    end subroutine Stabilize_init
    
    subroutine Stabilize_clear()
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
        real(kind=8) :: dif_wr
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
        real(kind=8) :: dif_wr
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
