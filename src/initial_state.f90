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
    
    subroutine initial_set(this, Latt, iseed)
        ! Set up initial state wave function by diagonalizing non-interacting Hamiltonian
        class(Initial), intent(inout) :: this
        class(kagomeLattice), intent(in) :: Latt
        integer, intent(inout) :: iseed
        
    ! Local variables
        complex(kind=8), allocatable :: HamT_initial(:,:), eigvecs(:,:)
        real(kind=8), allocatable :: eigvals(:)
        integer :: i
        
        allocate(HamT_initial(Ndim, Ndim), eigvecs(Ndim, Ndim))
        allocate(eigvals(Ndim))

        ! Build non-interacting Hamiltonian using existing function
        call def_hamT_initial(HamT_initial, Latt, iseed)
        
        ! Diagonalize Hamiltonian to get eigenvalues and eigenvectors
        ! Note: diag returns eigenvalues sorted from smallest to largest
        ! and eigenvectors are orthonormal by default
        call diag(HamT_initial, eigvecs, eigvals)
        
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
        
        deallocate(eigvals)
        deallocate(eigvecs, HamT_initial)
        return
    end subroutine initial_set
    
    subroutine def_hamT_initial(HamT_initial, Latt, iseed)
        ! Define non-interacting Hamiltonian for initial state
        ! This calls the existing def_hamT function from NonInteract module
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: HamT_initial
        class(kagomeLattice), intent(in) :: Latt
        integer, intent(inout) :: iseed
        real(kind=8), external :: ranf
        real(kind=8) :: twist_val
        complex(kind=8) :: Z
        integer :: ii, jj, nb, no, i, j
        
        select case (iniHam)
        case (3) ! twist with chemical potential bias
            HamT_initial = dcmplx(0.d0, 0.d0)
            Z = dcmplx( RT, 0.d0) 
            do ii = 1, Ndim
                do nb = 1, Nbond
                    jj = Latt%L_bonds(ii, nb)
                    HamT_initial(ii,jj) = HamT_initial(ii,jj) + Z
                    HamT_initial(jj,ii) = HamT_initial(jj,ii) + dconjg(Z)
                enddo
            enddo
            do ii = 1, Ndim
                no = Latt%dim_list(ii,2)
                if (no==1) HamT_initial(ii,ii) = HamT_initial(ii,ii) + dcmplx(-iniTwist, 0.d0)
                ! if (no==2) HamT_initial(ii,ii) = HamT_initial(ii,ii) + dcmplx( 0.d0, 0.d0)
                if (no==3) HamT_initial(ii,ii) = HamT_initial(ii,ii) + dcmplx( iniTwist, 0.d0)
            enddo
        case (2) ! twist with random hopping
            call def_hamT(HamT_initial, Latt)
            do ii = 1, Ndim
                do nb = 1, Nbond
                    jj = Latt%L_bonds(ii, nb)
                    twist_val = (ranf(iseed) - 0.5d0) * 2.d0 * iniTwist
                    HamT_initial(ii,jj) = HamT_initial(ii,jj) + dcmplx(twist_val, 0.d0)
                    HamT_initial(jj,ii) = HamT_initial(jj,ii) + dcmplx(twist_val, 0.d0)
                enddo
            enddo
        case (1) ! twist with random chemical potential
            call def_hamT(HamT_initial, Latt)
            do ii = 1, Ndim
                twist_val = (ranf(iseed) - 0.5d0) * 2.d0 * iniTwist
                HamT_initial(ii, ii) = HamT_initial(ii, ii) + dcmplx(twist_val, 0.d0)
            enddo
        case (0) ! no twist
            call def_hamT(HamT_initial, Latt)
        case (-1) ! all on A sites
            HamT_initial = dcmplx(0.d0, 0.d0)
            do ii = 1, Ndim
                if (Latt%dim_list(ii,2)==1) HamT_initial(ii, ii) = dcmplx(-1.d0, 0.d0)
            enddo
        case default
            call def_hamT(HamT_initial, Latt)
        end select
        
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
