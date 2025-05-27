# `Build_B0_CISD.f90` - Documentation

## Overview

The `Build_B0_CISD.f90` file contains the `Build_B0_CISD` subroutine, which is responsible for calculating a specific term (referred to as `cisd_b0`) in a Configuration Interaction Singles and Doubles (CISD) energy calculation. This term is a component of the total electronic energy and involves sums over one-electron and two-electron integrals, represented by `B1data` and `B2data` respectively. The calculation is constrained by the number of core, occupied, and virtual orbitals.

## Key Components

### `Build_B0_CISD`
*   **Description:** This subroutine calculates the `cisd_b0` term, which is a scalar value contributing to the CISD energy. It involves summing contributions from one-electron and two-electron terms based on orbital indices.
*   **Type:** Subroutine
*   **Key Arguments/Parameters:**
    *   `cisd_b0` (output, `real(8)`): The calculated B0 term for the CISD energy.
    *   `norb` (input, `integer`): Total number of orbitals. (Note: This variable is passed as an argument but not used in the provided code snippet.)
    *   `numfcore` (input, `integer`): Number of frozen core orbitals. These orbitals are excluded from the correlation calculation.
    *   `numfvirt` (input, `integer`): Number of frozen virtual orbitals. (Note: The loop limits suggest `numfvirt` might represent the count of active virtual orbitals, not frozen ones, as `b` goes up to `numocc+numfvirt` and `a` goes up to `numocc+numfvirt`. This might need clarification based on the broader project context.)
    *   `numocc` (input, `integer`): Number of occupied orbitals (excluding frozen core).
    *   `numvirt` (input, `integer`): Total number of virtual orbitals. (Note: This variable is passed as an argument but not used in the provided code snippet. `numfvirt` seems to define the upper limit for virtual orbital loops.)
    *   `B1data` (input, `type(B1)`): A derived type instance containing one-electron integral data, specifically `B1data%a%o(b,k)`.
    *   `B2data` (input, `type(B2)`): A derived type instance containing two-electron integral data, specifically `B2data%ab%m(a,k)%n(b,j)` and `B2data%aa%m(a,k)%n(b,j)`.
*   **Return Value/Outputs:** The primary output is the value assigned to the `cisd_b0` argument. The subroutine also prints "Building B0" and "B0 Final" along with the calculated `cisd_b0` value to standard output.

## Important Variables/Constants

*   **`sum_B0_1`** (`real(8)`): Local variable used to accumulate the one-electron contribution to `cisd_b0`. It sums the squares of elements from `B1data%a%o`.
*   **`sum_B0_2`** (`real(8)`): Local variable used to accumulate the two-electron contribution to `cisd_b0`. It sums the squares of elements from `B2data%ab%m` and `B2data%aa%m`.
*   **`cisd_b0`** (`real(8)`): Argument and the final calculated value. It is initialized to `1 + sum_B0_1 + sum_B0_2`.
*   **`numfcore`** (`integer`): Parameter defining the starting index for occupied orbitals involved in the calculation (e.g., `k=numfcore+1,numocc`).
*   **`numocc`** (`integer`): Parameter defining the upper limit for occupied orbitals and the starting point for virtual orbital indices (e.g., `b=numocc+1`).
*   **`numfvirt`** (`integer`): Parameter defining the range of virtual orbitals included in the sums (e.g., `b=numocc+1,numocc+numfvirt`).
*   **`B1data`** (`type(B1)`): Input data structure containing one-electron integrals.
*   **`B2data`** (`type(B2)`): Input data structure containing two-electron integrals.
*   **Loop counters** (`integer`): `a, b, k, j` (and others like `i, p, q, z, y, x, r, s, t, ioerror, size_l, size_lc, size_c, size_lcr, iter, energy_values, numatomic` which are declared but not used in this specific subroutine) are used as indices for iterating over orbitals.

## Usage Examples

There are no explicit Fortran code snippets demonstrating how to call this subroutine within the provided code. However, it would typically be called from a higher-level routine responsible for the overall CISD energy calculation, after `B1data` and `B2data` have been populated.

```fortran
! Hypothetical example of calling Build_B0_CISD
! Assumes B1_data_instance and B2_data_instance are populated
! and other parameters are defined.

! use TypeMod ! For B1 and B2 types
! implicit none
! real(8) :: calculated_b0
! integer :: n_orbitals, n_frozen_core, n_frozen_virtual, n_occupied, n_virtual
! type(B1) :: B1_data_instance
! type(B2) :: B2_data_instance

! ! ... (initialize parameters and data instances) ...

! call Build_B0_CISD(calculated_b0, n_orbitals, n_frozen_core, &
!                    n_frozen_virtual, n_occupied, n_virtual, &
!                    B1_data_instance, B2_data_instance)

! print *, "Calculated B0 value: ", calculated_b0
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `TypeMod`: This module is used via the `use TypeMod` statement. It is expected to define the derived types `B1` and `B2`, which are used for the `B1data` and `B2data` arguments, respectively. These types likely encapsulate the structure of one-electron and two-electron integrals.
    *   `FunctionMod`: This module is used via the `use FunctionMod` statement. Its specific role is not evident from the `Build_B0_CISD` subroutine alone, as no functions from this module are explicitly called within this subroutine. It might provide utility functions or data used elsewhere in the project or by components within `TypeMod`.
*   **External Libraries:**
    *   None are explicitly mentioned or used in this subroutine.
*   **Interactions with other components:**
    *   The `Build_B0_CISD` subroutine is likely called by a parent routine that orchestrates the CISD calculation.
    *   It reads data from `B1data` and `B2data` structures, which are presumably populated by other parts of the system responsible for calculating or retrieving these integrals.
    *   The calculated `cisd_b0` value is passed back to the calling routine and contributes to the total CISD energy.
```
