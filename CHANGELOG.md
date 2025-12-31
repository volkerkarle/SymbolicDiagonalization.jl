# Changelog

All notable changes to SymbolicDiagonalization.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Aqua.jl tests for package quality assurance
- Comprehensive edge case and error handling tests
- Examples directory with interactive pattern exploration script
- CHANGELOG.md for version history tracking

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
