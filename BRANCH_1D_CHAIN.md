# Branch: feature/1d-chain

## Overview

This branch modifies the main branch's kagome lattice implementation to a 1D chain lattice for simplified testing and validation of the PQMC algorithm.

## Key Differences from Main Branch

### 1. Lattice Geometry (`src/calc_basic.f90`)

**Main Branch (Kagome Lattice):**
- `Norb = 3` (three orbitals: A, B, C sublattices)
- `Nsub = 3` (three sublattices for observables)  
- `Nbond = 2` (two bonds per site: A→B/C, B→C/A, C→A/B)

**1D Chain Branch:**
- `Norb = 1` (single orbital per site)
- `Nsub = 1` (single sublattice)
- `Nbond = 1` (each site connects to one neighbor)

### 2. Neighbor Bond Structure (`src/lattice.f90`)

**Main Branch:** Complex kagome connectivity with orbital-dependent neighbor mapping:
```fortran
if (no==1) then      ! A --> B, C
if (no==2) then      ! B --> C, A  
if (no==3) then      ! C --> A, B
```

**1D Chain Branch:** Simple linear connectivity:
```fortran
! x --> x+1
n1  = Latt%inv_cell_list( npbc(ix+1, Nlx), npbc(iy, Nly) )
nn1 = Latt%inv_dim_list(n1, 1) ! only 1 orb
Latt%L_bonds(ii, 1) = nn1
```

### 3. Space-Time Bonds (`src/lattice.f90`)

**Main Branch:** 
```fortran
Latt%LT_bonds(iit, 1) = ... ! spatial neighbor 1
Latt%LT_bonds(iit, 2) = ... ! spatial neighbor 2  
Latt%LT_bonds(iit, 3) = ... ! time forward
Latt%LT_bonds(iit, 4) = ... ! time backward
```

**1D Chain Branch:**
```fortran
Latt%LT_bonds(iit, 1) = ... ! single spatial neighbor
Latt%LT_bonds(iit, 2) = ... ! time forward
Latt%LT_bonds(iit, 3) = ... ! time backward
```

### 4. Test Configuration (`test/paramC_sets.txt`)

**Main Branch:**
- `Nly = 6` (6×6 kagome lattice)
- `NlyTherm = 6` (thermalization lattice size)

**1D Chain Branch:** 
- `Nly = 1` (1D chain)
- `NlyTherm = 1` (1D thermalization)

## Physical Implications

- **Reduced complexity:** 1D chain eliminates triangular plaquette physics of kagome lattice
- **Simplified testing:** Easier to validate PQMC algorithm behavior and debug
- **Maintained algorithm:** All PQMC modifications (Stages 1-8) remain intact
- **Preserved functionality:** Observable calculations, measurements, and propagation algorithms unchanged

## Usage Notes

1. This branch maintains full PQMC algorithm compatibility
2. All build commands remain the same (`cd src && make`)
3. Configuration files automatically use 1D chain parameters
4. Results should be interpreted as 1D Hubbard model physics

## Merge Considerations

When merging back to main:
- Revert lattice parameters to kagome values
- Update test configurations to 2D geometry
- Verify all kagome-specific neighbor mappings are restored
- Test compilation and basic functionality

## Maintenance

- Keep PQMC algorithm changes synchronized with main branch
- Update this document when making geometry-specific modifications
- Test both 1D and kagome configurations before major releases