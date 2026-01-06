# Changelog

All notable changes to SymbolicDiagonalization.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Special angle simplification**: Trig functions at special angles (π/6, π/5, π/4, π/3, 2π/5, π/2, etc.) automatically simplify to exact algebraic values
  - `cos(π/3) → 1/2`, `cos(π/4) → √2/2`, `cos(2π/5) → (√5-1)/4`
  - Path Laplacian eigenvalues `2 - 2cos(kπ/n)` simplify to golden ratio expressions
  - New functions: `simplify_special_angles(expr)`, lookup tables `COS_SPECIAL_VALUES`, `SIN_SPECIAL_VALUES`
- Aqua.jl tests for package quality assurance
- Comprehensive edge case and error handling tests
- Examples directory with interactive pattern exploration script
- CHANGELOG.md for version history tracking
- New Lie group module organization under `src/patterns/lie/`:
  - `SO2.jl` - SO(2) rotation matrices and Kronecker products
  - `SO3.jl` - SO(3) rotation matrices (Rx, Ry, Rz)
  - `SO4.jl` - SO(4) detection
  - `SU2.jl` - SU(2) unitary matrices and Pauli matrices
  - `SU3.jl` - SU(3) and Gell-Mann matrices
  - `Sp.jl` - Symplectic group Sp(2n)
  - `algebras.jl` - Lie algebra generators (spin-j, angular momentum)
  - `common.jl` - Shared detection utilities

### Changed
- Refactored monolithic `diagonalize.jl` (1856 lines) into focused modules:
  - `src/structure.jl` - Structure detection and matrix properties
  - `src/patterns/graphs.jl` - Hypercube and strongly regular graph patterns
  - `src/patterns/circulant.jl` - Circulant and block circulant matrices
  - `src/patterns/tridiagonal.jl` - Toeplitz tridiagonal and anti-diagonal
  - `src/patterns/permutation.jl` - Permutation matrix eigenvalues
  - `src/patterns/kronecker.jl` - Kronecker product detection
- Updated Julia compatibility to 1.12 (for newer Symbolics.jl features)
- Improved `_with_timeout` implementation using Channel-based synchronization
- Extracted `_validate_matrix_input` helper to reduce code duplication
- Improved docstrings for `_try_factor_perfect_square` explaining intentional deferral

### Fixed
- **Critical**: Cube root of unity values in Cardano's formula (was `-0.5+0.5im`, now `-0.5+0.5*sqrt(3)*im`)
- Inconsistent error handling: now uses `ArgumentError` consistently instead of generic `error()`
- Type instability from untyped vectors (`[]` → `Any[]`)
- Removed duplicate `_is_one` helper in permutation.jl
- Removed dead code (`_std_basis` function)
- Race conditions in `_with_timeout` function

### Breaking Changes
- **Lie group function naming convention overhaul** - All functions now use consistent `GroupName_function` pattern:
  - SO(2): `R2(θ)` → `SO2_rotation(θ)`, `so2_kron` → `SO2_kron`, `so2_kron_eigenvalues` → `SO2_kron_eigenvalues`
  - SO(3): `Rx(θ)` → `SO3_Rx(θ)`, `Ry(θ)` → `SO3_Ry(θ)`, `Rz(θ)` → `SO3_Rz(θ)`
  - SU(2): `Ux(θ)` → `SU2_Ux(θ)`, `Uy(θ)` → `SU2_Uy(θ)`, `Uz(θ)` → `SU2_Uz(θ)`, `su2_kron` → `SU2_kron`, `su2_kron_eigenvalues` → `SU2_kron_eigenvalues`
  - SU(3): `su3_kron` → `SU3_kron`, `su3_kron_eigenvalues` → `SU3_kron_eigenvalues`, `su3_diagonal_trig` → `SU3_diagonal_trig`
  - Pauli matrices: `σx`, `σy`, `σz` → `pauli_x`, `pauli_y`, `pauli_z`
  - Gell-Mann matrices: `λ1`...`λ8` → `gellmann_1`...`gellmann_8`
- **Removed deprecated aliases** - No backward compatibility aliases provided; update code to use new names
- **Removed old Lie group files** - `lie_algebras.jl`, `lie_groups.jl`, `rotations.jl` replaced by `patterns/lie/` directory

## [0.1.0] - 2024-12-11

### Added
- Initial release of SymbolicDiagonalization.jl
- Core symbolic eigenvalue computation for matrices up to 4×4
- Closed-form formulas using Cardano's (cubic) and Ferrari's (quartic) methods
- Pattern detection for special matrix structures:
  - Diagonal and triangular matrices
  - Block-diagonal decomposition
  - Persymmetric matrices (half-size reduction)
  - Circulant and block circulant matrices
  - Toeplitz tridiagonal matrices
  - Permutation matrices
  - Hypercube graph adjacency matrices
  - Strongly regular graph adjacency matrices
- LinearAlgebra.jl interface extensions (`eigen`, `eigvals` for symbolic matrices)
- Comprehensive documentation with mathematical background
- Test suite with 240+ tests

### API
- `symbolic_eigenvalues(A)` - Compute eigenvalues only
- `symbolic_eigenpairs(A)` - Compute eigenvalues and eigenvectors
- `symbolic_diagonalize(A)` - Full diagonalization (P, D, λ)
- `characteristic_polynomial(A)` - Compute characteristic polynomial
- `symbolic_roots(coeffs)` - Solve polynomial equations symbolically

[Unreleased]: https://github.com/username/SymbolicDiagonalization.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/username/SymbolicDiagonalization.jl/releases/tag/v0.1.0
