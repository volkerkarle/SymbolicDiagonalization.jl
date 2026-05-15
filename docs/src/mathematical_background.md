# Mathematical Background

This page records the main mathematical reductions used by SymbolicDiagonalization.jl. It is intentionally concise; implementation details belong in the source code and examples belong in the [User Guide](user_guide.md).

## The Eigenvalue Barrier

For a matrix `A`, eigenvalues are roots of

```math
p(\lambda) = \det(A - \lambda I).
```

General formulas in radicals exist through degree 4. Abel-Ruffini rules out a general radical formula for degree 5 and above, so large symbolic matrices require additional structure.

The package therefore follows two paths:

- Reduce the effective degree using matrix structure.
- Otherwise solve the characteristic polynomial only up to degree 4.

## Polynomial Solvers

| Degree | Method |
|---:|---|
| 1 | Linear solve |
| 2 | Quadratic formula |
| 3 | Cardano formula |
| 4 | Ferrari/resolvent-cubic reduction |

These formulas are exact but expression size can grow quickly, especially for symbolic cubic and quartic coefficients.

## Structure Reductions

### Blocks

For `A = diag(B1, ..., Bk)`,

```math
\sigma(A) = \sigma(B_1) \cup \cdots \cup \sigma(B_k).
```

This can turn an unsolvable high-degree problem into several low-degree ones.

### Circulants

An `n x n` circulant matrix with first row `(c0, ..., c_{n-1})` is diagonalized by the DFT basis:

```math
\lambda_k = \sum_{j=0}^{n-1} c_j \omega^{jk}, \quad \omega = e^{2\pi i/n}.
```

The eigenvectors depend only on `n`; the entries affect only the Fourier coefficients.

### Kronecker Products

If `Av = lambda*v` and `Bw = mu*w`, then

```math
(A \otimes B)(v \otimes w) = \lambda\mu (v \otimes w).
```

This gives all eigenvalues of `A \otimes B` as pairwise products. A nested product of many 2x2 matrices still requires only repeated quadratic solves.

### Persymmetry

If `J` is the exchange matrix and `J*A*J == A`, the matrix commutes with reflection. The space splits into symmetric and antisymmetric subspaces, giving two smaller eigenproblems.

### Symmetric Toeplitz Tridiagonal Matrices

For diagonal `a` and off-diagonal `b`,

```math
\lambda_k = a + 2b\cos\left(\frac{k\pi}{n+1}\right), \quad k = 1, \ldots, n.
```

This follows from the sine eigenbasis and the Chebyshev recurrence.

## Group-Based Spectra

| Family | Reason the spectrum simplifies |
|---|---|
| `SO(2)` | rotations have characteristic polynomial `lambda^2 - 2cos(theta)lambda + 1` |
| `SO(3)` | one fixed axis gives eigenvalue `1`; the orthogonal plane is an `SO(2)` rotation |
| `SO(4)` | generic rotations split into two independent rotation planes |
| `SU(2)` | determinant one and unitarity force conjugate eigenvalues on the unit circle |
| Hadamard | `H*H' = nI`, so eigenvalues have magnitude `sqrt(n)` |
| DFT | `F^4 = n^2 I` for the unnormalized convention, giving fourth-root structure after scaling conventions are accounted for |

## Complexity Summary

| Pattern | Nominal degree | Effective work |
|---|---:|---|
| Generic `n x n` | `n` | radical formulas only for `n <= 4` |
| Block diagonal | `n` | largest block degree |
| Circulant | `n` | closed Fourier sum |
| Kronecker `A \otimes B` | `mn` | solve factors separately |
| Nested 2x2 Kronecker | `2^k` | `k` quadratic problems |
| Toeplitz tridiagonal | `n` | closed cosine formula |
| Hadamard/DFT | `n` | known finite spectra |

## References

1. P. J. Davis, *Circulant Matrices*, 1979.
2. B. C. Hall, *Lie Groups, Lie Algebras, and Representations*, 2015.
3. I. N. Herstein, *Topics in Algebra*, for the Abel-Ruffini theorem and solvability by radicals.
4. G. H. Golub and C. F. Van Loan, *Matrix Computations*, for structured eigenvalue problems.
