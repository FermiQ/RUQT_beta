# `Build_G_SD_Invert.f90` - Documentation

## Overview

The `Build_G_SD_Invert.f90` file contains the `Build_G_SD_Invert` subroutine. This subroutine is central to a quantum chemistry calculation, likely related to Green's function methods or many-body perturbation theory. Its primary role involves constructing and inverting a matrix (`G_C`) that represents a component of the Green's function or a related quantity. The calculation depends on whether singles, doubles, or CISD contributions are included, and it can interface with outputs from different quantum chemistry packages like GAMESS and Maple. It also handles reading molecular orbital (MO) energies, coefficients, and one- and two-electron integrals (`B1data`, `B2data`).

## Key Components

### `Build_G_SD_Invert`
*   **Description:** This subroutine calculates and inverts a central matrix `G_C`. The construction of this matrix involves several steps:
    1.  Reading MO energies, coefficients, and integral data (T1 and T2 amplitudes, which are stored in `B1data` and `B2data`) from files if it's the first iteration (`iter.eq.1`). The file reading logic differs based on whether `gamess` or `maple` flags are true.
    2.  If `doubles` is true and `use_b0` is true with `b0_type` as "cisd", it calls `Build_B0_CISD` to calculate the `cisd_b0` term.
    3.  It constructs a temporary Green's function matrix (`temp_gf`) in the MO basis. The terms included in `temp_gf` depend on whether `singles` or `doubles` are active, and the value of `use_b0` and `b0_type` for certain diagonal elements.
    4.  The `temp_gf` matrix is transformed into the atomic orbital (AO) basis to populate `G_S%en(ener_val)%gf`.
    5.  Finally, a specific block of this AO basis Green's function (`G_S%en(ener_val)%gf`) is extracted into `G_C`, adjusted by self-energy terms (`Sigma_l`, `Sigma_r`), and then inverted.
*   **Type:** Subroutine
*   **Key Arguments/Parameters:**
    *   `G_C` (output, `complex(8), dimension(:,:)`): The final inverted matrix.
    *   `Sigma_l`, `Sigma_r` (input, `complex(8), dimension(:,:)`): Left and Right self-energy matrices.
    *   `energy` (input, `real(8)`): Energy value at which calculations are performed.
    *   `size_l`, `size_c`, `size_lc`, `size_lcr` (input, `integer`): Integers likely defining matrix sizes or partitioning for block matrices.
    *   `norb` (input, `integer`): Total number of orbitals.
    *   `inputfile` (input, `character(len=40)`): Base name for input files.
    *   `numocc` (input, `integer`): Number of occupied orbitals.
    *   `numvirt` (input, `integer`): Number of virtual orbitals.
    *   `iter` (input, `integer`): Iteration number, controls initialization steps.
    *   `B1data` (inout, `type(B1)`): Derived type instance for one-electron integrals/amplitudes. Initialized if `iter.eq.1`.
    *   `B2data` (inout, `type(B2)`): Derived type instance for two-electron integrals/amplitudes. Initialized if `iter.eq.1`.
    *   `mo_ener` (inout, `real(8), dimension(:)`): Molecular orbital energies. Read from file if `iter.eq.1`.
    *   `mo_coeff`, `mo_coeff2` (inout, `real(8), dimension(:,:)`): Molecular orbital coefficients. Read from file if `iter.eq.1`.
    *   `doubles` (input, `logical`): Flag to include doubles contributions.
    *   `currentflag` (input, `logical`): Flag to control parts of the calculation flow (e.g., allocation of `temp_gf`).
    *   `energy_values` (input, `integer`): Number of energy points for `G_S`.
    *   `ener_val` (input, `integer`): Current energy index for `G_S`.
    *   `G_S` (inout, `type(energr)`): Derived type instance likely holding Green's function data at different energies.
    *   `corr_ener` (inout, `real(8)`): Correlation energy. Read from file or modified.
    *   `numatomic` (input, `integer`): Number of atomic orbitals/basis functions.
    *   `B0_Coeff` (input, `real(8), dimension(:)`): Coefficients for B0 term if `b0_type` is "rdm".
    *   `use_b0` (input, `logical`): Flag to use B0 term.
    *   `gamess` (input, `logical`): Flag indicating GAMESS is the source of input data.
    *   `maple` (input, `logical`): Flag indicating Maple is the source of input data.
    *   `numfcore` (input, `integer`): Number of frozen core orbitals.
    *   `numfvirt` (input, `integer`): Number of frozen virtual orbitals.
    *   `b0_type` (input, `character(len=40)`): Specifies the type of B0 term ("cisd" or "rdm").
*   **Return Value/Outputs:** The main output is the `G_C` matrix. The `G_S`, `B1data`, `B2data`, `mo_ener`, `mo_coeff`, `mo_coeff2`, and `corr_ener` variables can also be modified, especially during the first iteration.

## Important Variables/Constants

*   **`G_C`** (`complex(8), allocatable, dimension(:,:)`): The central complex matrix that is constructed and ultimately inverted. It represents a specific block of the Green's function after self-energy corrections.
*   **`G_S`** (`type(energr)`): A derived type variable that stores the Green's function (`G_S%en(i)%gf`) in the AO basis at various energy points. `G_S%en(ener_val)%gf` is computed and then used to extract `G_C`.
*   **`temp_gf`** (`real(8), allocatable, dimension(:,:)`): A temporary matrix holding the Green's function in the MO basis before transformation to the AO basis. Its construction depends on the `singles`/`doubles` flags and `use_b0` options.
*   **`B1data`** (`type(B1)`): Stores one-electron integrals or T1 amplitudes, read from files. Used in constructing `temp_gf`.
*   **`B2data`** (`type(B2)`): Stores two-electron integrals or T2 amplitudes, read from files. Used in constructing `temp_gf`.
*   **`mo_ener`** (`real(8), allocatable, dimension(:)`): Stores molecular orbital energies.
*   **`mo_coeff`** (`real(8), allocatable, dimension(:,:)`): Stores MO coefficients (AO to MO transformation).
*   **`cisd_b0`** (`real(8)`): Stores the B0 term calculated by `Build_B0_CISD` if `use_b0` is true and `b0_type` is "cisd". It's a scalar factor used in `temp_gf` construction.
*   **`numfcore`** (`integer`): Number of frozen core orbitals. Affects loop bounds and active orbital counts.
*   **`numocc`** (`integer`): Total number of occupied orbitals.
*   **`numfvirt`** (`integer`): Number of frozen virtual orbitals. Affects loop bounds and active orbital counts. (Note: `numvirt_act = numvirt - numfvirt` is the number of active virtuals).
*   **`numocc_act`** (`integer`): Number of active occupied orbitals (`numocc - numfcore`).
*   **`numvirt_act`** (`integer`): Number of active virtual orbitals (`numvirt - numfvirt`).
*   **`gamess`** (`logical`): Input flag. If true, the subroutine reads data from files formatted for GAMESS output (e.g., `inputfile // "T2"`, `inputfile // ".mo_dat"`).
*   **`maple`** (`logical`): Input flag. If true, the subroutine reads data from files formatted for Maple output (e.g., `inputfile // ".T2"`, `inputfile // ".scf_dat"`).
*   **`iter`** (`integer`): Iteration counter. `iter.eq.1` triggers initialization routines (reading files, allocating arrays).
*   **`doubles`** (`logical`): Input flag. If true, enables calculations including electron correlation effects from double excitations. Sets `singles` to `.false.`.
*   **`singles`** (`logical`): Local flag, set based on `doubles`. If `doubles` is false, `singles` is true.
*   **`use_b0`** (`logical`): Input flag. If true, specific B0 terms (e.g., `cisd_b0` or `B0_coeff`) are incorporated into the diagonal elements of `temp_gf`.
*   **`b0_type`** (`character(len=40)`): Input string. Determines the nature of the B0 term if `use_b0` is true (e.g., "cisd", "rdm").

## Usage Examples

There are no explicit Fortran code snippets demonstrating how to call this subroutine within the provided code. It is a core computational routine expected to be called by a higher-level driver program managing the overall Green's function or electronic structure calculation.

A hypothetical call might look like:
```fortran
! Assume all parameters and data types are properly declared and initialized
! For example:
! use TypeMod
! use InterfaceMod2 ! For energr type if not in TypeMod
! implicit none
! complex(8), dimension(size_c, size_c) :: G_calculated, Sigma_Left, Sigma_Right
! real(8) :: current_energy, correlation_e
! integer :: n_orb, n_occ, n_virt, iteration, n_atomic
! character(len=40) :: inp_file, b0_specification
! type(B1) :: B1_data_struct
! type(B2) :: B2_data_struct
! real(8), dimension(n_orb) :: mo_energies, b0_coeffs_rdm
! real(8), dimension(n_atomic, n_orb) :: mo_coeffs_mat, mo_coeffs_mat2
! logical :: use_doubles, current_run_flag, use_b0_term, is_gamess, is_maple
! type(energr) :: Greens_func_data
! integer :: num_energy_points, current_energy_idx, n_frozen_core, n_frozen_virtual
! integer :: s_l, s_c, s_lc, s_lcr ! matrix size parameters

! ... (Initialization of all arguments) ...

call Build_G_SD_Invert(G_calculated, Sigma_Left, Sigma_Right, current_energy, &
                       s_l, s_c, s_lc, s_lcr, n_orb, inp_file, n_occ, n_virt, iteration, &
                       B1_data_struct, B2_data_struct, mo_energies, mo_coeffs_mat, mo_coeffs_mat2, &
                       use_doubles, current_run_flag, num_energy_points, current_energy_idx, &
                       Greens_func_data, correlation_e, n_atomic, b0_coeffs_rdm, use_b0_term, &
                       is_gamess, is_maple, n_frozen_core, n_frozen_virtual, b0_specification)

! ... (Use G_calculated) ...
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `InterfaceMod2`: This module is used via `use InterfaceMod2`. It likely provides matrix manipulation routines such as `matmul_dgemm2` (presumably a wrapper for DGEMM for matrix multiplication) and matrix inversion functions (`inv_real` for real matrices, `inv` for complex matrices). It might also define the `energr` derived type if not from `TypeMod`.
    *   `TypeMod`: This module is used via `use TypeMod`. It is expected to define the derived types `B1` and `B2`, which are used for the `B1data` and `B2data` arguments, respectively. These types encapsulate one-electron and two-electron integral/amplitude data. It might also define the `energr` type.
    *   `FunctionMod`: This module is used via `use FunctionMod`. Its specific role is not fully evident from this subroutine alone, but it might provide other mathematical or utility functions used in the broader project or by the other used modules.
*   **External Libraries:**
    *   None are explicitly mentioned as direct external library calls, but `matmul_dgemm2` suggests a possible indirect dependency on a BLAS/LAPACK library for efficient matrix operations.
*   **Interactions with other components:**
    *   **`Build_B0_CISD` Subroutine:** This subroutine is conditionally called if `doubles.eqv..true.`, `use_b0.eqv..true.`, and `trim(b0_type).eq."cisd"`. It calculates the `cisd_b0` term, which is then used in the construction of `temp_gf`. The path to `Build_B0_CISD.f90` is not specified but it's an internal project dependency.
    *   **Input Files:** The subroutine reads data from several files during the first iteration (`iter.eq.1`). The names of these files are constructed from `inputfile` and depend on the `gamess` or `maple` flags:
        *   GAMESS: `inputfile // "T2"` (for T1/T2 amplitudes) and `inputfile // ".mo_dat"` (for MO coefficients and energies).
        *   Maple: `inputfile // ".T2"` (for T1/T2 amplitudes and correlation energy) and `inputfile // ".scf_dat"` (for MO coefficients and energies).
    *   **Calling Environment:** The subroutine is designed to be called by a higher-level routine that manages the overall calculation, providing it with energy points, matrix size definitions, and control flags. It updates the `G_C` matrix and potentially `G_S` which are then used by the caller.
```
