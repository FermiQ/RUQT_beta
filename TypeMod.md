# `TypeMod.f90` - Documentation

## Overview

The `TypeMod.f90` file defines a Fortran module named `TypeMod`. The primary purpose of this module is to declare several custom data types (derived types) that are used throughout the larger quantum chemistry or physics simulation project. These data structures are designed to organize and manage complex data, such as one- and two-electron integrals or amplitudes, Green's functions, and orbital indices, in a structured and efficient manner.

## Key Components

The key components of this module are the derived type definitions.

### `inner`
*   **Description:** A simple type designed to hold a 2D allocatable array of double-precision real numbers. It serves as a building block for the `outer` type, typically representing the innermost part of a nested data structure for multi-dimensional arrays (e.g., the `n(b,j)` part of `B2data%aa%m(a,k)%n(b,j)`).
*   **Structure:**
    *   `real(8), allocatable, dimension(:,:) :: n`

### `outer`
*   **Description:** This type acts as a container for the `inner` type, creating a nested structure. It holds a 2D allocatable array where each element is of type `inner`. This allows for representing four-dimensional data structures, commonly found in quantum chemistry for storing two-electron integrals or related quantities (e.g., the `m(a,k)` part of `B2data%aa%m(a,k)%n(b,j)`).
*   **Structure:**
    *   `type(inner), allocatable, dimension(:,:) :: m`

### `B2`
*   **Description:** This type is designed to store two-electron integral data or related quantities (like T2 amplitudes in coupled-cluster theory or components of a two-particle Green's function). It groups data based on spin combinations: `aa` (alpha-alpha), `ab` (alpha-beta), `bb` (beta-beta), and potentially a `nospin` component for spin-unrestricted or averaged cases. Each component is of type `outer`, allowing for the storage of four-index quantities (e.g., `<ij|kl>` or `t_ij^ab`).
*   **Structure:**
    *   `type(outer) :: aa`
    *   `type(outer) :: ab`
    *   `type(outer) :: bb`
    *   `type(outer) :: nospin`

### `virtorb`
*   **Description:** A type designed to hold a 2D allocatable array of double-precision real numbers, likely representing terms involving virtual orbitals and occupied orbitals (e.g., the `o(b,k)` part of `B1data%a%o(b,k)`). It serves as a building block for the `B1` type.
*   **Structure:**
    *   `real(8), allocatable, dimension(:,:) :: o`

### `B1`
*   **Description:** This type is designed to store one-electron integral data or related quantities (like T1 amplitudes in coupled-cluster theory or components of a one-particle Green's function). It groups data into `a` (alpha spin) and `b` (beta spin) components. Each component is of type `virtorb`, allowing for the storage of two-index quantities (e.g., `f_ia` or `t_i^a`).
*   **Structure:**
    *   `type(virtorb) :: a`
    *   `type(virtorb) :: b`

### `ingreen`
*   **Description:** This type is a container for a Green's function (`gf`), represented as a 2D allocatable array of double-precision real numbers. This likely stores the Green's function matrix at a single energy point.
*   **Structure:**
    *   `real(8), allocatable, dimension(:,:) :: gf`

### `energr`
*   **Description:** This type is designed to store Green's function data at multiple energy points. It contains an allocatable array (`en`), where each element is of type `ingreen`. This allows for storing a series of Green's function matrices, typically one for each energy in a specified range.
*   **Structure:**
    *   `type(ingreen), allocatable, dimension(:) :: en`

## Important Variables/Constants

This module also declares module-level allocatable arrays for storing orbital indices.

*   **`Virt_Index`**: `integer, allocatable, dimension(:)`
    *   **Description:** An allocatable 1D integer array intended to store the indices of virtual orbitals. This can be useful for mapping or iterating over specific sets of virtual orbitals in calculations.
*   **`Occ_Index`**: `integer, allocatable, dimension(:)`
    *   **Description:** An allocatable 1D integer array intended to store the indices of occupied orbitals. This can be useful for mapping or iterating over specific sets of occupied orbitals.

## Usage Examples

These types are used by declaring variables of their kind in other parts of the program.

```fortran
! Example of declaring a variable of a custom type
! from TypeMod

module MyCalculationModule
  use TypeMod ! Make the types available
  implicit none

  ! Declare a variable to hold two-electron terms
  type(B2) :: two_electron_amplitudes

  ! Declare a variable to hold one-electron terms
  type(B1) :: one_electron_amplitudes

  ! Declare a variable to hold Green's functions at various energies
  type(energr) :: green_function_data

  ! Declare variables for storing orbital indices
  integer, dimension(:), allocatable :: my_virtual_indices, my_occupied_indices

contains

  subroutine setup_calculation()
    ! Allocate and populate Virt_Index and Occ_Index if needed
    ! allocate(Virt_Index(num_virtual_orbitals))
    ! allocate(Occ_Index(num_occupied_orbitals))
    ! ... logic to fill these arrays ...

    ! Allocate components of the derived types
    ! For example, for B2data%aa%m(a,k)%n(b,j):
    ! allocate(two_electron_amplitudes%aa%m(num_virtual_a, num_occupied_k))
    ! do k_idx = 1, num_occupied_k
    !   do a_idx = 1, num_virtual_a
    !     allocate(two_electron_amplitudes%aa%m(a_idx,k_idx)%n(num_virtual_b, num_occupied_j))
    !   end do
    ! end do
    ! ... and so on for other components and types ...
  end subroutine setup_calculation

end module MyCalculationModule
```

## Dependencies and Interactions

*   **Dependencies:**
    *   The `TypeMod` module itself does not have explicit dependencies on other custom modules for its definitions. It uses intrinsic Fortran types (`real(8)`, `integer`, `allocatable`, `dimension`).
*   **Interactions:**
    *   This module is fundamental to the project. Its defined types (`B1`, `B2`, `energr`, etc.) are expected to be used extensively by many other Fortran files and modules within the system.
    *   Subroutines and functions in other parts of the project (e.g., in `InterfaceMod.f90`, `RUQT.f90`, `Build_G_SD_Invert.f90`, `Build_B0_CISD.f90`) will declare arguments of these types to pass complex data structures.
    *   The main program and various computational subroutines will instantiate variables of these types to store and manipulate data related to electronic structure, integrals, amplitudes, and Green's functions.
    *   The module-level variables `Virt_Index` and `Occ_Index` would be populated and used by routines that need to perform operations specifically on virtual or occupied orbital subspaces.
```
