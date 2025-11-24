module GlobalK_mod
    use DQMC_Model_mod
    use ProcessMatrix
    implicit none
    
    public
    private :: ratioK_fermion
    logical :: warned_global = .false.
    
contains
    real(kind=8) function ratioK_fermion(Gr, phi_new, ii, ntau) result(ratio_fermion)
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi_new    
        integer, intent(in) :: ii, ntau

        ratio_fermion = 1.d0
        return
    end function ratioK_fermion
    
    subroutine GlobalK_prop_L(Prop, ratio_fermion, phi_new, nt) ! placeholder: global updates disabled
        class(Propagator), intent(inout) :: Prop
        real(kind=8), intent(inout) :: ratio_fermion
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi_new
        integer, intent(in) :: nt
        if (.not. warned_global) then
            write(6,*) "GlobalK_prop_L: global update path disabled in rank-1 flow"
            warned_global = .true.
        endif
        return
    end subroutine GlobalK_prop_L
    
    subroutine GlobalK_prop_R(Prop, ratio_fermion, phi_new, nt)
        class(Propagator), intent(inout) :: Prop
        real(kind=8), intent(inout) :: ratio_fermion
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi_new
        integer, intent(in) :: nt
        if (.not. warned_global) then
            write(6,*) "GlobalK_prop_R: global update path disabled in rank-1 flow"
            warned_global = .true.
        endif
        return
    end subroutine GlobalK_prop_R
end module GlobalK_mod
