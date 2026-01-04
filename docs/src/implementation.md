# Implementation Details

Technical overview of algorithms and architecture.

## Architecture

```
Input Matrix
    ↓
Structure Detection
    ↓ (pattern found?)
    ├─ Yes → Pattern-Specific Solver
    └─ No  → Characteristic Polynomial → Root Solver (degree ≤ 4)
                                              ↓
                                        Eigenvectors (if requested)
```

## File Organization

```
src/
├── SymbolicDiagonalization.jl  # Module, exports
├── charpoly.jl                 # Characteristic polynomial (Bareiss)
├── roots.jl                    # Root solvers (degrees 1-4)
├── rref.jl                     # RREF, nullspace
├── simplify.jl                 # Trigonometric simplification
├── structure.jl                # Structure detection, utilities
├── diagonalize.jl              # Public API
└── patterns/
    ├── circulant.jl            # Circulant, block circulant
    ├── graphs.jl               # Hypercube, strongly regular
    ├── kronecker.jl            # Kronecker products
    ├── permutation.jl          # Permutation matrices
    ├── tridiagonal.jl          # Toeplitz tridiagonal, special 5×5
    └── lie/                    # Lie group/algebra patterns
        ├── common.jl           # Shared detection utilities
        ├── algebras.jl         # Lie algebra generators, spin-j
        ├── SO2.jl              # SO(2) rotations, Kronecker products
        ├── SO3.jl              # SO(3) rotations (Rx, Ry, Rz)
        ├── SO4.jl              # SO(4) detection
        ├── SU2.jl              # SU(2), Pauli matrices
        ├── SU3.jl              # SU(3), Gell-Mann matrices
        └── Sp.jl               # Symplectic groups Sp(2n)
```

## Key Algorithms

### Characteristic Polynomial

**Algorithm**: Bareiss (fraction-free Gaussian elimination)

**Why**: Standard elimination causes expression explosion; Bareiss keeps expressions polynomial-sized.

**Coefficient extraction**: Taylor series at λ=0 using derivatives (avoids polynomial division).

### Root Solvers

| Degree | Method | Notes |
|--------|--------|-------|
| 1 | Direct | $x = -c_0/c_1$ |
| 2 | Quadratic | Standard formula with discriminant simplification |
| 3 | Cardano | Depress → compute cube roots → back-substitute |
| 4 | Ferrari | Resolve cubic → decompose into quadratics |

**Quartic warning**: Produces ~13.5 MB expressions for fully symbolic 4×4 matrices.

### Structure Detection

Detection order (most specific first):

1. Diagonal / Triangular
2. Block-diagonal (greedy connected components)
3. Circulant / Block circulant
4. Toeplitz tridiagonal
5. Permutation matrix
6. Kronecker product
7. Anti-diagonal / Persymmetric
8. Special 5×5 patterns

All detection uses `_issymzero(x)` which conservatively returns `false` if zero cannot be proven.

### Pattern Solvers

| Pattern | Algorithm | Complexity |
|---------|-----------|------------|
| Circulant | DFT of first row | $O(n^2)$ |
| Toeplitz tridiagonal | Cosine formula | $O(n)$ |
| Kronecker | Solve factors, multiply | $O(T(m) + T(n))$ |
| Permutation | Cycle decomposition → roots of unity | $O(n)$ |
| Block-diagonal | Solve blocks recursively | $\sum O(n_i)$ |

### Eigenvector Computation

Two methods:
- **Adjugate** (n ≤ 3): Columns of adj(A - λI) are in nullspace
- **RREF nullspace** (general): Row reduce A - λI, extract free variables

## Expression Management

### Problem

Symbolic expressions grow exponentially without management.

### Solutions

1. **Aggressive simplification**: Applied after every intermediate step
2. **Complexity limits**: `max_terms` parameter (default 10,000)
3. **Timeout**: `timeout` parameter (default 300s)

### Symbolic Square Root

Custom implementation for `Complex{Num}` since Julia's `sqrt` has boolean checks that fail symbolically:
```
sqrt(a + bi) = sqrt((r+a)/2) + i·sqrt((r-a)/2)  where r = |a + bi|
```

## Performance

### Complexity

| Operation | Complexity |
|-----------|------------|
| Bareiss determinant | $O(n^3)$ |
| Root solving | $O(\text{expr\_size}^{\text{degree}})$ |
| RREF | $O(n^3 \times \text{expr\_size})$ |
| Circulant eigenvalues | $O(n^2)$ |
| Tridiagonal eigenvalues | $O(n)$ |

### Memory (fully symbolic matrices)

| Size | Eigenvalues |
|------|-------------|
| 2×2 | ~1 KB |
| 3×3 | ~100 KB |
| 4×4 | ~50 MB |

### Threading

**Not supported**: Symbolics.jl uses task-local storage that is not thread-safe. Use process-based parallelism (Distributed.jl) if needed.

## Extending

### Adding a Pattern

1. Add detection function `_is_pattern(mat)` returning `nothing` or detection info
2. Add solver `_pattern_eigenvalues(mat, info)`
3. Add to detection cascade in `diagonalize.jl`
4. Add tests in `test/test_patterns.jl`

See [Contributing](contributing.md) for details.
