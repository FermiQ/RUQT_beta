# `InterfaceMod.f90` - Documentation

## Overview

The `InterfaceMod.f90` file serves as a central repository for interface declarations in a quantum chemistry or physics simulation program. It defines three distinct modules: `InterfaceMod`, `InterfaceMod2`, and `FunctionMod`. These modules declare the explicit interfaces for a variety of subroutines and functions, ensuring type safety and clear contracts for how these procedures are called from other parts of the program. The declared procedures cover a range of functionalities, including setting flags, reading input data, retrieving results from quantum chemistry packages, performing core computational steps (like building Green's functions or specific energy terms), and common mathematical operations (matrix manipulations, Fermi function calculation).

## Key Components

This file defines three modules, each containing a set of interface blocks for various subroutines and functions.

### `InterfaceMod`
*   **Description:** This module groups interfaces for a broad range of high-level operations. These include setting up calculation parameters, reading data from various quantum chemistry packages (QChem, Molcas, GAMESS, PySCF, libint), defining system partitions and electrode interactions (especially for transport calculations using Wide-Band Limit - WBL), reading main program inputs, and calling core computational routines.
*   **Type:** Module
*   **Key Subroutines/Functions Declared:**
    *   **`Flag_set`**
        *   **Description:** Sets various logical flags (e.g., `cisd_flag`, `rdm_flag`, `hf_flag`, `dft_flag`) and flags indicating the quantum chemistry package used (e.g., `qchem`, `gamess`) based on input strings.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `inputcode`, `functional` (input strings); various logical flags (output).
        *   **Return Value/Outputs:** Modifies the logical flag arguments.
    *   **`Get_HF_QChem`**, **`Get_HF_Molcas`**, **`Get_HF_GAMESS`**, **`Get_HF_PySCF`**
        *   **Description:** These subroutines appear to retrieve Hartree-Fock (HF) calculation results (e.g., two-electron integrals `H_two`, overlap matrix `Smat`) from output files of different quantum chemistry packages (QChem, Molcas, GAMESS, PySCF respectively).
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `inputfile` (input string); `norb`, `numatomic` (input integers); `H_two`, `Smat` (output, allocatable real matrices).
        *   **Return Value/Outputs:** Populates `H_two` and `Smat` matrices.
    *   **`Get_HF_libint`**
        *   **Description:** Retrieves HF data, potentially using the libint library, providing one-electron (`H_one`), two-electron (`TwoIntsCompact`), and overlap (`Smat`) integrals, along with MO coefficients.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `inputfile` (input string); `norb`, `numact` (input integers); `H_one`, `Smat`, `mo_coeff`, `OneInts` (output, allocatable real matrices); `TwoIntsCompact` (output, allocatable real array).
        *   **Return Value/Outputs:** Populates integral and coefficient matrices/arrays.
    *   **`Calculate_Coupling_MoleculeWBL`**
        *   **Description:** Calculates coupling strengths (`Coupling_R`, `Coupling_L`) for a molecule in the Wide-Band Limit, based on a local density of states at the Fermi level.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `localden_fermi` (input real); `Coupling_R`, `Coupling_L` (output reals).
        *   **Return Value/Outputs:** Sets `Coupling_R` and `Coupling_L`.
    *   **`Electrodes_MoleculeWBL`**
        *   **Description:** Computes self-energy matrices (`Sigma_l`, `Sigma_r`) for a molecule connected to electrodes in the Wide-Band Limit, using coupling strengths and overlap information.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `Coupling_R`, `Coupling_L` (input reals); `Smat` (input real matrix); `size_c`, `size_lc`, `size_lcr` (input integers defining matrix dimensions).
        *   **Return Value/Outputs:** Populates `Sigma_L` and `Sigma_R` (output, allocatable complex matrices).
    *   **`PartitionHS_MoleculeWBL`**
        *   **Description:** Partitions the overlap matrix (`Smat`) into blocks corresponding to different parts of the system (e.g., electrode-molecule connections) for WBL calculations.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `Smat` (input real matrix); `size_l`, `size_r`, `size_c`, etc. (input integers for partitioning); `Smat_el`, `Smat_er` (output real matrices).
        *   **Return Value/Outputs:** Populates partitioned overlap matrices.
    *   **`Electrodes_MetalWBL`**
        *   **Description:** Computes self-energy matrices for metallic electrodes in the WBL, using local densities of states and Hamiltonian/overlap matrix blocks.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `Smat_re`, `Smat_le`, `H_Two_le`, `H_Two_re`, etc. (input real matrices for electrode properties); `localden_fermi_l`, `localden_fermi_r` (input reals); `size_c`, `size_l`, `size_r`, etc. (input integers).
        *   **Return Value/Outputs:** Populates `Sigma_L` and `Sigma_R`.
    *   **`PartitionHS_MetalWBL`**
        *   **Description:** Partitions system Hamiltonian (`H_Two`) and overlap (`Smat`) matrices into blocks for different regions (left electrode, right electrode, central region) in metal WBL calculations.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `Smat`, `H_Two` (input real matrices); `size_l`, `size_r`, `size_c`, etc. (input integers); various partitioned matrices (output).
        *   **Return Value/Outputs:** Populates partitioned Hamiltonian and overlap matrices.
    *   **`ReadInput`**
        *   **Description:** Reads a wide range of simulation parameters from an input file. These include orbital counts, energy ranges, voltage settings, convergence criteria, choice of method, and flags.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `inputfile` (input string); numerous parameters defining the calculation setup (mostly outputs).
        *   **Return Value/Outputs:** Populates many simulation control parameters.
    *   **`Build_G_SD_Invert`**
        *   **Description:** Interface to the main computational subroutine that builds and inverts the Green's function matrix, considering singles and doubles contributions. (This subroutine's detailed documentation is in `Build_G_SD_Invert.md`).
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `G_C` (output complex matrix); `Sigma_l`, `Sigma_r` (input complex matrices); `energy` (input real); `B1data`, `B2data` (input `TypeMod` types), and many other control and data parameters.
        *   **Return Value/Outputs:** Outputs the computed `G_C` matrix.

### `InterfaceMod2`
*   **Description:** This module seems to be a smaller, more focused collection of interfaces, specifically for routines that might be called by components defined in `InterfaceMod` or other core parts of the calculation.
*   **Type:** Module
*   **Key Subroutines/Functions Declared:**
    *   **`Build_B0_CISD`**
        *   **Description:** Interface to the subroutine that calculates the B0 term for CISD (Configuration Interaction Singles and Doubles) energy. (This subroutine's detailed documentation is in `Build_B0_CISD.md`).
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `B0_coeff` (output real); `norb`, `numfcore`, `numfvirt`, `numocc`, `numvirt` (input integers for orbital definitions); `B1data`, `B2data` (input `TypeMod` types).
        *   **Return Value/Outputs:** Outputs the calculated `B0_coeff`.

### `FunctionMod`
*   **Description:** This module provides interfaces for a set of utility functions, primarily focused on mathematical operations such as matrix algebra (multiplication, inversion, adjoint) and the Fermi-Dirac distribution function. It also includes functions for index manipulation.
*   **Type:** Module
*   **Key Subroutines/Functions Declared:**
    *   **`matmul_dgemm`**
        *   **Description:** Matrix multiplication for two double-precision real matrices.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `leftmatrix`, `rightmatrix` (input real matrices).
        *   **Return Value/Outputs:** Returns the resulting product matrix (real).
    *   **`matmul_zgemm`**
        *   **Description:** Matrix multiplication for two double-precision complex matrices.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `leftmatrix`, `rightmatrix` (input complex matrices).
        *   **Return Value/Outputs:** Returns the resulting product matrix (complex).
    *   **`matmul_dgemm2`**
        *   **Description:** Matrix multiplication for two double-precision real matrices, storing the result in a provided output matrix.
        *   **Type:** Subroutine
        *   **Key Arguments/Parameters:** `leftmatrix`, `rightmatrix` (input real matrices); `matsol` (output real matrix).
        *   **Return Value/Outputs:** Populates `matsol`.
    *   **`inv`**
        *   **Description:** Inversion of a complex matrix.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `A` (input complex matrix).
        *   **Return Value/Outputs:** Returns the inverse of `A` (complex matrix).
    *   **`inv_real`**
        *   **Description:** Inversion of a real matrix.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `A` (input real matrix).
        *   **Return Value/Outputs:** Returns the inverse of `A` (real matrix).
    *   **`adjoint`**
        *   **Description:** Calculates the adjoint (conjugate transpose) of a complex matrix.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `A` (input complex matrix); `norb` (integer, likely matrix dimension).
        *   **Return Value/Outputs:** Returns the adjoint of `A` (complex matrix).
    *   **`fermi_function`**
        *   **Description:** Calculates the value of the Fermi-Dirac distribution function for a given energy, Fermi energy, and thermal energy (KT).
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `energy`, `fermi_energy`, `KT` (input reals).
        *   **Return Value/Outputs:** Returns the Fermi function value (real).
    *   **`FirstIndex`**
        *   **Description:** Likely calculates a packed 1D index from two orbital indices `i, k`, often used for symmetric or anti-symmetric two-electron integrals.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `i`, `k` (input integers).
        *   **Return Value/Outputs:** Returns an `integer(8)` index.
    *   **`CompositeIndex`**
        *   **Description:** Likely calculates a packed 1D index from two pairs of orbital indices (`ik`, `jl`), often used for storing two-electron integrals in a compact form.
        *   **Type:** Function
        *   **Key Arguments/Parameters:** `ik`, `jl` (input `integer(8)` which are themselves packed indices).
        *   **Return Value/Outputs:** Returns an `integer(8)` composite index.

## Important Variables/Constants

This file primarily contains interface declarations. Thus, it does not define file-level variables or constants that significantly affect its own behavior. The "variables" of importance are the arguments to the subroutines and functions, which dictate how these procedures operate. These are detailed within each component's description above.

## Usage Examples

The modules defined in this file (`InterfaceMod`, `InterfaceMod2`, `FunctionMod`) are intended to be used by other Fortran files within the project. This is achieved by including a `use` statement at the beginning of the Fortran file that needs to call one of the declared procedures.

```fortran
! Example of using InterfaceMod to call ReadInput
program MySimulation
  use InterfaceMod
  implicit none

  character(len=100) :: inputFile_name = "sim_input.dat"
  ! ... (declare other variables needed by ReadInput) ...
  integer :: norbs, n_occ, n_virt
  ! ... (many other parameters)

  ! Call ReadInput using the interface defined in InterfaceMod
  call ReadInput(inputFile_name, norbs, ... , n_occ, n_virt, ... )

  ! ... (rest of the simulation code) ...
end program MySimulation

! Example of using FunctionMod to invert a matrix
module MyCalculations
  use FunctionMod
  use TypeMod ! If matrix types are involved with types from here
  implicit none

contains

  subroutine process_matrix(matrix_in, matrix_out)
    complex(8), dimension(:,:), allocatable, intent(in) :: matrix_in
    complex(8), dimension(:,:), allocatable, intent(out) :: matrix_out
    
    ! Assuming matrix_in is square and invertible
    allocate(matrix_out(size(matrix_in,1), size(matrix_in,2)))

    matrix_out = inv(matrix_in) ! inv is from FunctionMod
    
  end subroutine process_matrix

end module MyCalculations
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `TypeMod`: The interfaces for `Build_G_SD_Invert` (in `InterfaceMod`) and `Build_B0_CISD` (in `InterfaceMod2`) explicitly use derived types (`B1`, `B2`, `energr`) which are expected to be defined in a module named `TypeMod`. This creates a dependency: `TypeMod` must be compiled before or alongside `InterfaceMod.f90`, and linked for programs using these interfaces.
*   **External Libraries:**
    *   While not directly specified in the interfaces, functions like `matmul_dgemm`, `matmul_zgemm`, `inv`, and `inv_real` often wrap calls to optimized linear algebra libraries such as BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra Package). The actual implementation of these functions (not shown in the interface file) would determine the direct dependency.
*   **Interactions with other components:**
    *   **Calling Modules/Programs:** Any Fortran file that `use`s `InterfaceMod`, `InterfaceMod2`, or `FunctionMod` will be able to call the subroutines and functions whose interfaces are declared therein. This is the primary mode of interaction.
    *   **Implementation Files:** This file (`InterfaceMod.f90`) only provides the "shape" of the functions and subroutines. The actual Fortran code that implements these procedures must exist elsewhere in the project and be linked during compilation to create a working executable. For example, the implementation of `ReadInput` would read from a specified file and parse its contents. The implementation of `Build_G_SD_Invert` would perform the complex calculations to build the Green's function.
```
