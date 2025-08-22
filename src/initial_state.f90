module MakeInitialState
    use MyLattice
    use MyMats
    use NonInteract
    implicit none
    
    public
    
    type :: Initial
        complex(kind=8), dimension(:,:), allocatable :: PR  ! Right trial wave function (Ndim, 1)
        complex(kind=8), dimension(:,:), allocatable :: PL  ! Left trial wave function (1, Ndim)
        real(kind=8) :: energy_ground                       ! Ground state energy
        real(kind=8) :: energy_gap                          ! Energy gap to first excited state
        real(kind=8) :: norm_factor                         ! Normalization factor
    contains
        procedure :: make        => initial_make
        procedure :: set         => initial_set
        procedure :: get_PR      => initial_get_PR
        procedure :: get_PL      => initial_get_PL
        procedure :: get_gap     => initial_get_gap
        final :: initial_clear
    end type Initial
    
contains
    
    subroutine initial_make(this)
        ! Allocate memory for Initial structure
        class(Initial), intent(inout) :: this
        
        allocate(this%PR(Ndim, 1))
        allocate(this%PL(1, Ndim))
        this%PR = dcmplx(0.d0, 0.d0)
        this%PL = dcmplx(0.d0, 0.d0)
        this%energy_ground = 0.d0
        this%energy_gap = 0.d0
        this%norm_factor = 1.d0
        
        return
    end subroutine initial_make
    
    subroutine initial_clear(this)
        ! Deallocate memory for Initial structure
        type(Initial), intent(inout) :: this
        
        if (allocated(this%PR)) deallocate(this%PR)
        if (allocated(this%PL)) deallocate(this%PL)
        
        return
    end subroutine initial_clear
    
    subroutine initial_set(this, Latt)
        ! Set up initial state wave function by diagonalizing non-interacting Hamiltonian
        class(Initial), intent(inout) :: this
        class(kagomeLattice), intent(in) :: Latt
        
        ! Local variables
        complex(kind=8), dimension(Ndim, Ndim) :: HamT, eigvecs
        real(kind=8), dimension(Ndim) :: eigvals
        integer :: i
        
        ! Build non-interacting Hamiltonian using existing function
        call def_hamT_initial(HamT, Latt)
        
        ! Diagonalize Hamiltonian to get eigenvalues and eigenvectors
        ! Note: diag returns eigenvalues sorted from smallest to largest
        ! and eigenvectors are orthonormal by default
        call diag(HamT, eigvecs, eigvals)
        
        ! Ground state energy (first eigenvalue)
        this%energy_ground = eigvals(1)
        
        ! Energy gap (difference between first excited and ground state)
        if (Ndim > 1) then
            this%energy_gap = eigvals(2) - eigvals(1)
        else
            this%energy_gap = 0.d0
        endif
        
        ! Set up PR (ground state wave function as column vector)
        do i = 1, Ndim
            this%PR(i, 1) = eigvecs(i, 1)
        enddo
        
        ! Set up PL (conjugate transpose of PR as row vector)
        do i = 1, Ndim
            this%PL(1, i) = dconjg(this%PR(i, 1))
        enddo
        
        ! Calculate normalization factor to check if it's 1
        this%norm_factor = sqrt(real(sum(dconjg(this%PR(:,1)) * this%PR(:,1))))
        
        ! Only renormalize if norm factor is significantly different from 1
        if (abs(this%norm_factor - 1.d0) > 1.d-12) then
            this%PR(:,1) = this%PR(:,1) / this%norm_factor
            this%PL(1,:) = this%PL(1,:) / this%norm_factor
        endif
        
        return
    end subroutine initial_set
    
    subroutine def_hamT_initial(HamT, Latt)
        ! Define non-interacting Hamiltonian for initial state
        ! This calls the existing def_hamT function from NonInteract module
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: HamT
        class(kagomeLattice), intent(in) :: Latt
        
        call def_hamT(HamT, Latt)
        
        return
    end subroutine def_hamT_initial
    
    function initial_get_PR(this) result(PR_out)
        ! Get PR wave function
        class(Initial), intent(in) :: this
        complex(kind=8), dimension(Ndim, 1) :: PR_out
        
        PR_out = this%PR
        return
    end function initial_get_PR
    
    function initial_get_PL(this) result(PL_out)
        ! Get PL wave function
        class(Initial), intent(in) :: this
        complex(kind=8), dimension(1, Ndim) :: PL_out
        
        PL_out = this%PL
        return
    end function initial_get_PL
    
    function initial_get_gap(this) result(gap)
        ! Get energy gap
        class(Initial), intent(in) :: this
        real(kind=8) :: gap
        
        gap = this%energy_gap
        return
    end function initial_get_gap
    
end module MakeInitialState