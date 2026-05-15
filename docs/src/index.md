# SymbolicDiagonalization.jl

SymbolicDiagonalization.jl computes closed-form eigenvalues and eigenvectors for symbolic Julia matrices when the problem is algebraically solvable or has exploitable structure.

The package extends the `LinearAlgebra` interface:

```julia
using LinearAlgebra, Symbolics, SymbolicDiagonalization

@variables a b
A = [a b; b a]

eigvals(A)  # [a + b, a - b]
```

## What It Does

General symbolic diagonalization is only possible up to degree 4 in radicals. Larger matrices need structure. SymbolicDiagonalization.jl detects common structures and reduces the eigenproblem before falling back to generic polynomial solving.

Use it when your matrix is:

- Small, dense, and symbolic: up to 4x4 by general formulas
- Block diagonal or reducible by symmetry
- Circulant, symmetric circulant, Toeplitz tridiagonal, permutation-like, Hadamard, DFT, or Kronecker-structured
- Built from supported rotation/unitary constructors such as `SO2_rotation`, `SO3_Rz`, or `SU2_Uz`

## API

```julia
eigvals(A)                         # LinearAlgebra eigenvalues
eigen(A)                           # LinearAlgebra eigen decomposition
symbolic_eigenvalues(A)            # returns (values, characteristic_polynomial, lambda)
symbolic_eigenpairs(A)             # returns (pairs, characteristic_polynomial, lambda)
symbolic_diagonalize(A)            # returns (P, D, pairs), with A = P * D * inv(P)
```

Common constructors:

```julia
SO2_rotation(theta)
SO3_Rx(theta), SO3_Ry(theta), SO3_Rz(theta)
SU2_Ux(theta), SU2_Uy(theta), SU2_Uz(theta)
hadamard_matrix(n)
dft_matrix(n; normalized=false)
cartan_matrix_A(n)
```

## Documentation Map

- [User Guide](user_guide.md): choosing APIs, managing complexity, and validating results
- [Pattern Library](pattern_library.md): supported structures and formulas
- [Mathematical Background](mathematical_background.md): why the reductions work

## Limits

- Generic 5x5 and larger symbolic matrices are not solvable in radicals in general.
- Quartic formulas can be correct but too large to be useful.
- Pattern detection is conservative: a matrix may be mathematically structured but still fall back to the generic path if the structure is not represented in a detectable form.
