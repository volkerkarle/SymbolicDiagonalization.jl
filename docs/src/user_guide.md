# User Guide

This page focuses on usage decisions rather than listing every supported matrix family. For the full catalog, see the [Pattern Library](pattern_library.md).

## Choose the Right Entry Point

Use the `LinearAlgebra` functions unless you need intermediate symbolic data:

```julia
using LinearAlgebra, Symbolics, SymbolicDiagonalization

@variables a b
A = [a b; b a]

eigvals(A)
```

The lower-level functions expose extra data:

```julia
values, poly, lambda = symbolic_eigenvalues(A)
pairs, poly, lambda = symbolic_eigenpairs(A)
P, D, pairs = symbolic_diagonalize(A)
```

Use `eigvals(A)` for eigenvalues, `eigen(A)` for vectors, and `symbolic_diagonalize(A)` only when you need an explicit `P * D * inv(P)` factorization.

## Prefer Detectable Structure

The package works best when structure is explicit in the matrix representation.

```julia
@variables c0 c1 c2 c3
C = [c0 c1 c2 c3;
     c3 c0 c1 c2;
     c2 c3 c0 c1;
     c1 c2 c3 c0]

eigvals(C)
```

For Lie-group matrices, prefer package constructors over manual reconstruction:

```julia
R = SO2_rotation(theta)        # preferred
R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
```

Both may work, but constructors avoid ambiguity and usually lead to cleaner expressions.

## Manage Expression Growth

Symbolic eigenvalues can become large even for 3x3 or 4x4 matrices. Prefer reducing the problem before asking for a full eigensystem.

Practical order of operations:

1. Try `eigvals(A)` before `eigen(A)`.
2. Look for block, Kronecker, circulant, tridiagonal, or Lie-group structure.
3. Use constructors for supported matrix families.
4. Add a timeout for exploratory symbolic work.

```julia
eigvals(A; timeout=60)
eigvals(A; max_terms=1000)
```

If a result is too large to inspect, use substitution to validate it numerically before simplifying further.

```julia
using Symbolics: substitute
vals = eigvals(A)
substitute(vals[1], Dict(a => 1.0, b => 2.0))
```

## Eigenvectors and Parallelism

Eigenvectors are more expensive than eigenvalues because each eigenspace requires a nullspace computation. When Julia workers are available, `symbolic_eigenpairs` can parallelize independent eigenvector computations.

```julia
using Distributed
addprocs(4)
@everywhere using SymbolicDiagonalization

pairs, poly, lambda = symbolic_eigenpairs(A)
```

This helps most when eigenvalues are already known but eigenvectors are expensive.

## Development Tests

The default test suite skips some slow Lie-group cases. Run them explicitly when changing Lie-group or Kronecker logic:

```bash
SYMBOLICDIAG_SLOW_TESTS=1 julia --project -e 'include("test/test_lie_groups.jl")'
```
