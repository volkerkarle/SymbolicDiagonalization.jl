# Pattern Library

SymbolicDiagonalization.jl tries specialized reductions before using the generic characteristic-polynomial solver. This page is a compact reference for those reductions.

## Detection Order

Earlier matches avoid more expensive or less structured algorithms.

| Priority | Pattern | Reduction |
|---:|---|---|
| 1 | Diagonal and triangular | Read diagonal entries |
| 2 | Block diagonal | Solve blocks independently |
| 3 | Lie groups | Use trace, determinant, or angle invariants |
| 4 | Kronecker products | Multiply factor eigenvalues |
| 5 | Hadamard and DFT matrices | Use known spectra |
| 6 | Circulant and symmetric circulant | Use Fourier diagonalization |
| 7 | Tridiagonal Toeplitz and Cartan type A | Use cosine/Chebyshev formulas |
| 8 | Persymmetric | Split into symmetric and antisymmetric sectors |
| 9 | Generic size <= 4 | Solve characteristic polynomial in radicals |

## Structural Patterns

| Pattern | Form | Spectrum |
|---|---|---|
| Diagonal | `diag(d1, ..., dn)` | `d1, ..., dn` |
| Triangular | upper or lower triangular | diagonal entries |
| Block diagonal | `diag(B1, ..., Bk)` | union of spectra of the blocks |
| Persymmetric | `J * A * J == A` | spectra of two smaller projected matrices |

Block and persymmetric reductions are useful because they lower the effective polynomial degree. A 6x6 block matrix made of 2x2 blocks is solved as three quadratic problems, not one sextic.

## Finite-Group Patterns

| Pattern | Constructor or form | Eigenvalue formula |
|---|---|---|
| Circulant | rows cyclically shifted from first row `c` | `lambda_k = sum_j c_j * omega^(j*k)` |
| Symmetric circulant | circulant with palindromic first row | real cosine form of the DFT spectrum |
| Permutation matrix | disjoint cycles | roots of unity for each cycle length |
| Quaternion structure | `Q8_invariant_matrix` or equivalent 2x2 form | scalar part plus/minus quaternion norm |
| Graph Laplacians | `path_laplacian`, `cycle_laplacian` | standard path/cycle spectra |
| Coxeter and Cartan | `cartan_matrix_A`, ..., `cartan_matrix_G2` | closed forms for type A and G2; other types may fall back |
| Hadamard | `hadamard_matrix(n)` | `+/-sqrt(2^n)` for order `2^n` |
| DFT | `dft_matrix(n)` | fourth roots of `n`, with multiplicities depending on `n mod 4` |

Example:

```julia
using LinearAlgebra, Symbolics, SymbolicDiagonalization

@variables a b c d
C = [a b c d;
     d a b c;
     c d a b;
     b c d a]

eigvals(C)
```

## Lie-Group Patterns

| Family | Constructors | Spectrum |
|---|---|---|
| `SO(2)` | `SO2_rotation(theta)` | `exp(+/-im*theta)` |
| `SO(3)` | `SO3_Rx`, `SO3_Ry`, `SO3_Rz` | `1`, `exp(+/-im*theta)` |
| `SO(4)` | block/double rotations | two conjugate rotation pairs |
| `SU(2)` | `SU2_Ux`, `SU2_Uy`, `SU2_Uz` | `exp(+/-im*theta/2)` |
| Symplectic | detected symplectic form | reciprocal-pair constraint, not a complete closed-form solver |

The package also exports Pauli matrices, spin generators, and Gell-Mann matrices. The Gell-Mann helpers are constructors/generators; SU(3) matrices use the generic solver unless another supported structure is detected.

## Kronecker Products

For eigenpairs `(lambda, v)` of `A` and `(mu, w)` of `B`:

```math
(A \otimes B)(v \otimes w) = (\lambda \mu)(v \otimes w)
```

So the spectrum of `A \otimes B` is all products `lambda_i * mu_j`. Nested Kronecker products apply the same rule recursively.

```julia
@variables x y theta
K = kron([x y; y x], SO2_rotation(theta))
eigvals(K)
```

Specialized constructors `SO2_kron` and `SU2_kron` return cleaner trigonometric forms for repeated rotation/unitary products.

## Tridiagonal Patterns

For the symmetric Toeplitz tridiagonal matrix with diagonal `a` and off-diagonal `b`,

```math
\lambda_k = a + 2b\cos\left(\frac{k\pi}{n+1}\right), \quad k = 1, \ldots, n.
```

The type-A Cartan matrix is the special case `a = 2`, `b = -1`.

## Generic Fallback

If no pattern matches, the package forms the characteristic polynomial and solves it in radicals for degrees 1 through 4. For generic degree 5 and above, no radical formula exists in general, so the package raises an error instead of returning misleading expressions.
