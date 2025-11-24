module ProcessMatrix
    use MyLattice
    use MakeInitialState
    implicit none
    public
    
    type :: Propagator ! allocated in one spin-orbital sector for PQMC
        complex(kind=8), dimension(:,:), allocatable :: UUR  ! (Ndim, 1) - B(tau,0)P
        complex(kind=8), dimension(:,:), allocatable :: UUL  ! (1, Ndim) - P^dagger B(2theta,tau)
        complex(kind=8) :: overlap                        ! scalar P_L^dagger P_R
        real(kind=8) :: Xmaxm, Xmeanm
    contains
        procedure :: make => Prop_make
        procedure :: asgn => Prop_assign
        final :: Prop_clear
    end type Propagator
    
    type :: PropGreen
        complex(kind=8), dimension(:,:), allocatable :: Gr00, Grtt, Grt0, Gr0t
        real(kind=8) :: Xmaxm(3), Xmeanm(3)
    contains
        procedure :: make => Propgr_make
        procedure :: reset => Propgr_reset
        final :: Propgr_clear
    end type PropGreen
    
    type, public :: WrapList
        complex(kind=8), dimension(:,:,:), allocatable :: URlist, ULlist  ! (Ndim,1,1:Nst), (1,Ndim,1:Nst)
        ! Note: VRlist, VLlist, DRlist, DLlist removed for PQMC algorithm
    contains
        procedure :: make => Wrlist_make
        procedure :: asgn => Wrlist_assign
        final :: Wrlist_clear
    end type WrapList
    
contains
    subroutine Prop_make(this, Init_obj)
        class(Propagator), intent(inout) :: this
        class(Initial), intent(in) :: Init_obj
        allocate(this%UUR(Ndim, 1), this%UUL(1, Ndim))
        this%UUR = Init_obj%PR  ! Initialize with trial wave function
        this%UUL = Init_obj%PL  ! Initialize with trial wave function
        this%overlap = sum(this%UUL(1, 1:Ndim) * this%UUR(1:Ndim, 1))
        this%Xmaxm = 0.d0; this%Xmeanm = 0.d0
        return
    end subroutine Prop_make
    
    subroutine Prop_assign(this, that)
        class(Propagator), intent(inout) :: this
        class(Propagator), intent(in) :: that
        this%UUL = that%UUL
        this%UUR = that%UUR
        this%overlap = that%overlap
        this%Xmaxm = that%Xmaxm; this%Xmeanm = that%Xmeanm
        return
    end subroutine Prop_assign
    
    subroutine Prop_clear(this)
        type(Propagator), intent(inout) :: this
        deallocate(this%UUR, this%UUL)
        return
    end subroutine Prop_clear
    
    subroutine Propgr_make(this)
        class(PropGreen), intent(inout) :: this
        allocate(this%Gr00(Ndim, Ndim), this%Grtt(Ndim, Ndim))
        allocate(this%Grt0(Ndim, Ndim), this%Gr0t(Ndim, Ndim))
        this%Gr00 = ZKRON
        this%Grtt = ZKRON
        this%Gr0t = dcmplx(1.d0, 0.d0)
        this%Grt0 = ZKRON
        this%Xmaxm = 0.d0
        this%Xmeanm = 0.d0
        return
    end subroutine Propgr_make

    subroutine Propgr_reset(this, Gr)
        class(PropGreen), intent(inout) :: this
        complex(kind=8), dimension(:,:), intent(in) :: Gr
        this%Gr00 = Gr
        this%Grtt = Gr
        this%Grt0 = Gr
        this%Gr0t = - ( ZKRON - Gr )
        return
    end subroutine Propgr_reset
    
    subroutine Propgr_clear(this)
        type(PropGreen), intent(inout) :: this
        deallocate(this%Gr00, this%Grtt, this%Grt0, this%Gr0t)
        return
    end subroutine Propgr_clear
    
    subroutine Wrlist_make(this)
        class(WrapList), intent(inout) :: this
        allocate(this%URlist(Ndim, 1, max(Nst, 1)), this%ULlist(1, Ndim, max(Nst, 1)))
        this%URlist = dcmplx(0.d0, 0.d0); this%ULlist = dcmplx(0.d0, 0.d0)
        return
    end subroutine Wrlist_make

    subroutine Wrlist_assign(this, that)
        class(WrapList), intent(inout) :: this
        class(WrapList), intent(in) :: that
        this%URlist = that%URlist; this%ULlist = that%ULlist
        return
    end subroutine Wrlist_assign
    
    subroutine Wrlist_clear(this)
        type(WrapList), intent(inout) :: this
        deallocate(this%URlist, this%ULlist)
        return
    end subroutine Wrlist_clear
end module ProcessMatrix
