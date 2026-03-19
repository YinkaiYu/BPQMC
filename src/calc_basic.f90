module CalcBasic ! Global parameters
    implicit none

! constants
    real(kind=8),           parameter           :: Zero = 1.0d-10
    real(kind=8),           parameter           :: PI = acos(-1.d0)
    real(kind=8),           parameter           :: upbound = 1.0d+200
! lattice parameters
    integer,                parameter           :: Naux  = 2  ! flavor number of auxiliary field, respectively for U1 term and U2 term
    integer,                public,     save    :: Norb = 3
    integer,                public,     save    :: Nsub = 3
    integer,                public,     save    :: Nbond = 2
    integer,                public              :: Nlx, Nly, NlxTherm, NlyTherm
    integer,                public              :: Lq, LqTherm
    integer,                public              :: Ndim, NdimTherm
    real(kind=8),           public              :: Dtau
    real(kind=8),           public              :: Beta
    integer,                public              :: Ltrot, LtrotTherm
    character(len=16),      public,     save    :: lattice_type = 'kagome'
! Hamiltonian parameters
    real(kind=8),           public,     save    :: RT
    real(kind=8),           public,     save    :: RU1, RU2
    integer,                public,     save    :: Nbos
    real(kind=8),           public,     save    :: imbalance ! chemical potential imbalance 
! update parameters
    real(kind=8),           public              :: shiftLoc
    logical,                public              :: is_global ! HMC switch; false => local update
    integer,                public              :: Nfrog
    integer,                public              :: NfrogJitter
    real(kind=8),           public              :: hmc_dt
    real(kind=8),           public              :: hmc_mass
    real(kind=8),           public              :: hmc_mass_spatial_uniform
    real(kind=8),           public              :: hmc_mass_spatial_shell1
    integer,                public              :: hmc_block_tau
    integer,                public              :: hmc_block_sites
! initial state parameters
    integer,                public              :: iniType ! type of initial phonon field configuration
    real(kind=8),           public              :: iniAmpl ! Gaussian amplitude of initial phonon fields
    real(kind=8),           public              :: iniBias(Naux) ! initial phonon field balance position
    integer,                public              :: iniHam  ! initial Hamiltonian twist mode
    real(kind=8),           public              :: iniTwist ! twist amplitude for initial Hamiltonian
! process control parameters
    logical,                    public          :: is_tau ! whether to calculate time-sliced Green function
    integer,                    public          :: Nthermal ! calculate time-sliced Green function from the (Nthermal + 1)-th bin
    logical,                    public          :: is_warm ! bosonic warm-up switch
    integer,                    public          :: Nwarm
    real(kind=8),               public          :: shiftWarm(Naux)
    integer,                    public          :: Nst  ! number of stored imaginary-time slices (0..Ltrot)
    integer,                    public          :: Nwrap
    integer,                    public          :: Nbin
    integer,                    public          :: Nsweep
    integer,                    public          :: ISIZE, IRANK, IERR
! utilities
    real(kind=8),               public          :: norm_threshold
    
contains
    subroutine read_input()
        include 'mpif.h'
        character(len=256) :: first_line
        character(len=256) :: hmc_line
        integer :: ios, ios_hmc, ios_ini
        if (IRANK == 0) then
            open(unit=20, file='paramC_sets.txt', status='unknown')
            read(20,'(A)') first_line
            call read_lattice_header(first_line)
            read(20,*) Nlx, Nly, Ltrot, Beta
            read(20,*) NlxTherm, NlyTherm, LtrotTherm
            read(20,*) Nwrap, Nbin, Nsweep, shiftLoc
            read(20,*) is_tau, Nthermal
            read(20,*) is_warm, Nwarm, shiftWarm(1), shiftWarm(2)
            read(20,'(A)') hmc_line
            is_global = .false.
            Nfrog = 0
            hmc_dt = 0.d0
            NfrogJitter = 0
            hmc_mass = 1.d0
            hmc_mass_spatial_uniform = 0.d0
            hmc_mass_spatial_shell1 = 0.d0
            hmc_block_tau = 0
            hmc_block_sites = 0
            read(hmc_line, *, iostat=ios_hmc) is_global, Nfrog, hmc_dt, NfrogJitter, hmc_mass, hmc_block_tau, hmc_block_sites
            if (ios_hmc /= 0) then
                read(hmc_line, *, iostat=ios_hmc) is_global, Nfrog, hmc_dt, NfrogJitter, hmc_mass, hmc_block_tau
                if (ios_hmc == 0) hmc_block_sites = 0
            endif
            if (ios_hmc /= 0) then
                read(hmc_line, *, iostat=ios_hmc) is_global, Nfrog, hmc_dt, NfrogJitter, hmc_mass
                if (ios_hmc == 0) then
                    hmc_block_tau = 0
                    hmc_block_sites = 0
                endif
            endif
            if (ios_hmc /= 0) then
                read(hmc_line, *, iostat=ios_hmc) is_global, Nfrog, hmc_dt, NfrogJitter
                if (ios_hmc == 0) then
                    hmc_mass = 1.d0
                    hmc_block_tau = 0
                    hmc_block_sites = 0
                endif
            endif
            if (ios_hmc /= 0) then
                read(hmc_line, *, iostat=ios_hmc) is_global, Nfrog, hmc_dt
                if (ios_hmc == 0) then
                    NfrogJitter = 0
                    hmc_mass = 1.d0
                    hmc_block_tau = 0
                    hmc_block_sites = 0
                endif
            endif
            if (ios_hmc /= 0) then
                read(hmc_line, *, iostat=ios_ini) iniType, iniAmpl, iniBias(1), iniBias(2)
                if (ios_ini /= 0) then
                    write(6,*) "Failed to parse HMC/input line:", trim(hmc_line)
                    stop
                endif
            else
                read(20,*) iniType, iniAmpl, iniBias(1), iniBias(2)
            endif
            read(20,*) iniHam, iniTwist, imbalance
            close(20)
            call read_env_real("BPQMC_HMC_MASS_SPATIAL_UNIFORM", hmc_mass_spatial_uniform)
            call read_env_real("BPQMC_HMC_MASS_SPATIAL_SHELL1", hmc_mass_spatial_shell1)
        endif 
!   MPI process: parallelization
        call MPI_BCAST(Beta, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RT, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RU1, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RU2, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nbos, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(shiftLoc, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hmc_dt, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hmc_mass, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hmc_mass_spatial_uniform, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hmc_mass_spatial_shell1, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hmc_block_tau, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hmc_block_sites, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(shiftWarm, Naux, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniAmpl, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniBias, Naux, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniType, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniHam, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniTwist, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(imbalance, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(lattice_type, len(lattice_type), MPI_CHARACTER, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nlx, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nly, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Ltrot, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlxTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlyTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(LtrotTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwrap, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nbin, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwarm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nfrog, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NfrogJitter, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nthermal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nsweep, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_tau, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_warm, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_global, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        norm_threshold = 1.d-6
        return
    end subroutine read_input
    
    subroutine Params_set()
        call set_lattice_layout()
        Lq = Nlx * Nly
        LqTherm = NlxTherm * NlyTherm
        Ndim = Lq * Norb
        NdimTherm = LqTherm * Norb
        if (Ltrot > 0) then
            Dtau = Beta / dble(Ltrot)
        else
            Dtau = 0.d0
        endif
        if (Nwrap <= 0) then
            write(6,*) "Nwrap must be positive"; stop
        endif
        if (is_global) then
            if (Ltrot <= 0) then
                write(6,*) "Ltrot must be positive in HMC mode"; stop
            endif
            if (Nfrog <= 0) then
                write(6,*) "Nfrog must be positive in HMC mode"; stop
            endif
            if (NfrogJitter < 0) then
                write(6,*) "NfrogJitter must be non-negative in HMC mode"; stop
            endif
            if (hmc_dt <= 0.d0) then
                write(6,*) "hmc_dt must be positive in HMC mode"; stop
            endif
            if (hmc_mass <= 0.d0) then
                write(6,*) "hmc_mass must be positive in HMC mode"; stop
            endif
            if (hmc_mass_spatial_uniform < 0.d0) then
                write(6,*) "hmc_mass_spatial_uniform must be non-negative in HMC mode"; stop
            endif
            if (hmc_mass_spatial_shell1 < 0.d0) then
                write(6,*) "hmc_mass_spatial_shell1 must be non-negative in HMC mode"; stop
            endif
            if (hmc_block_tau < 0) then
                write(6,*) "hmc_block_tau must be non-negative in HMC mode"; stop
            endif
            if (hmc_block_tau > Ltrot) then
                write(6,*) "hmc_block_tau cannot exceed Ltrot in HMC mode"; stop
            endif
            if (hmc_block_sites < 0) then
                write(6,*) "hmc_block_sites must be non-negative in HMC mode"; stop
            endif
            if (hmc_block_sites > Ndim) then
                write(6,*) "hmc_block_sites cannot exceed Ndim in HMC mode"; stop
            endif
        endif
        if (Ltrot <= 0) then
            Nst = 0
        else
            Nst = (Ltrot + Nwrap - 1) / Nwrap  ! ceil(Ltrot / Nwrap)
        endif
        return
    end subroutine Params_set

    subroutine read_lattice_header(first_line)
        character(len=*), intent(in) :: first_line
        integer :: ios

        read(first_line, *, iostat=ios) RT, RU1, RU2, Nbos
        if (ios == 0) then
            lattice_type = 'kagome'
        else
            lattice_type = normalize_lattice_name(first_line)
            read(20,*) RT, RU1, RU2, Nbos
        endif
        return
    end subroutine read_lattice_header

    subroutine set_lattice_layout()
        select case (trim(lattice_type))
        case ('kagome')
            Norb = 3
            Nsub = 3
            Nbond = 2
        case ('triangular')
            Norb = 1
            Nsub = 1
            Nbond = 3
        case default
            write(6,*) "Unsupported lattice type:", trim(lattice_type)
            stop
        end select
        return
    end subroutine set_lattice_layout

    pure function normalize_lattice_name(text) result(name)
        character(len=*), intent(in) :: text
        character(len=16) :: name
        integer :: idx, code

        name = adjustl(text)
        do idx = 1, len_trim(name)
            code = ichar(name(idx:idx))
            if (code >= ichar('A') .and. code <= ichar('Z')) then
                name(idx:idx) = achar(code + 32)
            endif
        enddo
        select case (trim(name))
        case ('tri', 'triangle', 'triangular')
            name = 'triangular'
        case ('kag', 'kagome')
            name = 'kagome'
        end select
        return
    end function normalize_lattice_name

    subroutine read_env_real(name, value)
        character(len=*), intent(in) :: name
        real(kind=8), intent(inout) :: value
        character(len=64) :: env_value
        integer :: env_status, ios

        call get_environment_variable(name, env_value, status=env_status)
        if (env_status /= 0) return
        if (len_trim(env_value) == 0) return
        read(env_value, *, iostat=ios) value
        if (ios /= 0) return
        return
    end subroutine read_env_real

    pure logical function is_triangular_lattice()
        is_triangular_lattice = trim(lattice_type) == 'triangular'
        return
    end function is_triangular_lattice

    pure logical function is_kagome_lattice()
        is_kagome_lattice = trim(lattice_type) == 'kagome'
        return
    end function is_kagome_lattice

    real(kind=8) function norm_diff_vec(vec1, vec2, n)
        complex(kind=8), intent(in) :: vec1(:), vec2(:)
        integer, intent(in) :: n
        integer :: i
        real(kind=8) :: acc
        acc = 0.d0
        do i = 1, n
            acc = acc + abs(vec1(i) - vec2(i))**2
        enddo
        norm_diff_vec = sqrt(acc)
        return
    end function norm_diff_vec
    
    integer function nranf(iseed, N)
        integer, intent(inout) :: iseed
        integer, intent(in) :: N
        real (kind=8), external :: ranf
        nranf  = nint(ranf(iseed)*dble(N) + 0.5)
        if (nranf .lt. 1 ) nranf = 1
        if (nranf .gt. N ) nranf = N 
        return
    end function nranf
    
    ! define the periodic boundary condition in x,y-directions.
    integer function npbc(nr, L)
        integer, intent(in) :: nr, L
        npbc = nr
        if (nr .gt. L) npbc = nr - L
        if (nr .lt. 1) npbc = nr + L
        return
    end function npbc
    
    pure function sqr_vec(vec)
        real(kind=8) :: sqr_vec
        real(kind=8), dimension(Naux), intent(in) :: vec
        integer :: ns
        sqr_vec = 0.d0
        do ns = 1, Naux
            sqr_vec = sqr_vec + vec(ns) * vec(ns)
        enddo
        return
    end function sqr_vec
    
    subroutine write_info()
        character(len=64) :: title
        if (IRANK == 0) then
            open (unit=50, file='info.txt', status='unknown', action="write")
            write(title, '(A,A,A)') 'DQMC for boson Hubbard on ', trim(lattice_type), ' lattice'
            write(50,*) '========================='
            write(50,*) trim(title)
            write(50,*) 'Lattice type                                   :', trim(lattice_type)
            write(50,*) 'Number of orbitals per cell                    :', Norb
            write(50,*) 'Number of bonds per site                       :', Nbond
            write(50,*) 'Linear lengh Lx                                :', Nlx
            write(50,*) 'Linear lengh Ly                                :', Nly
            write(50,*) 'Hopping t                                      :', RT
            write(50,*) 'Hubbard U1                                     :', RU1
            write(50,*) 'Hubbard U2                                     :', RU2
            write(50,*) 'boson particle number                          :', Nbos
            if (is_global) then
                write(50,*) 'Sampler                                        :', 'HMC'
                write(50,*) 'Leapfrog steps                                :', Nfrog
                write(50,*) 'Leapfrog jitter                               :', NfrogJitter
                write(50,*) 'Leapfrog step size                             :', hmc_dt
                write(50,*) 'Leapfrog mass                                  :', hmc_mass
                write(50,*) 'Spatial-uniform leapfrog mass                  :', hmc_mass_spatial_uniform
                write(50,*) 'Lowest-shell leapfrog mass                     :', hmc_mass_spatial_shell1
                write(50,*) 'HMC tau block size                             :', hmc_block_tau
                write(50,*) 'HMC site block size                            :', hmc_block_sites
            else
                write(50,*) 'Sampler                                        :', 'Local'
                write(50,*) 'Local update auxiliary field magnitude Shift   :', shiftLoc
            endif
            if (is_warm) then
            write(50,*) '# Warm                                         :', Nwarm
            if (.not. is_global) then
                write(50,*) 'Thermalize auxiliary field magnitude for U1    :', shiftWarm(1)
                write(50,*) 'Thermalize auxiliary field magnitude for U2    :', shiftWarm(2)
            endif
            endif
            write(50,*) 'Choosing initial distribution type             :', iniType
            write(50,*) 'Magnitude of Gaussian distribution             :', iniAmpl
            write(50,*) 'Initial Hamiltonian mode (0=none,1=twist,...)  :', iniHam
            write(50,*) 'Initial Hamiltonian twist amplitude            :', iniTwist
            write(50,*) 'Hamiltonian chemical potential imbalance       :', imbalance
            write(50,*) 'Beta                                           :', Beta
            write(50,*) 'Trotter number                                 :', Ltrot
            write(50,*) '=>Dtau                                         :', Dtau
            write(50,*) 'N_Ortho                                        :', Nwrap
            write(50,*) '# Bins                                         :', Nbin
            write(50,*) '# Nsweep in one bin                            :', Nsweep
            write(50,*) '# Cores                                        :', ISIZE
        endif
        return
    end subroutine write_info
    
end module CalcBasic
