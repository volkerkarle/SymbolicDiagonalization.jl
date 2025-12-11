# User Guide

Practical guide to using `SymbolicDiagonalization.jl` effectively.

## Table of Contents

- [Getting Started](#getting-started)
- [Basic Workflows](#basic-workflows)
- [Working with Special Patterns](#working-with-special-patterns)
- [Performance Optimization](#performance-optimization)
- [Troubleshooting](#troubleshooting)
- [Advanced Usage](#advanced-usage)

---

## Getting Started

### Installation

```bash
# From package directory
cd SymbolicDiagonalization.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

### Basic Setup

```julia
using Symbolics
using SymbolicDiagonalization
using LinearAlgebra

# Define symbolic variables
@variables a b c d

# Create a matrix
M = [a b; b a]

# Compute eigenvalues
λ = eigvals(M)  # [a+b, a-b]
```

### Quick Reference

```julia
# Eigenvalues only (fast)
λ = eigvals(A)

# Eigenvalues + eigenvectors
E = eigen(A)
E.values    # eigenvalues
E.vectors   # eigenvectors as columns

# With options
λ = eigvals(A, structure=:hermitian, timeout=60)
```

---

## Basic Workflows

### Workflow 1: Analyze Small Matrices (≤ 4×4)

**When to use**: You have a symbolic matrix up to 4×4 and need exact eigenvalues.

**Steps**:
```julia
@variables a b c d

# Create your matrix
M = [a  b  0  0;
     b  a  0  0;
     0  0  c  d;
     0  0  d  c]

# Get eigenvalues (structure auto-detected)
λ = eigvals(M)
# Output: [a+b, a-b, c+d, c-d]

# Get eigenvectors too
E = eigen(M)
```

**Tips**:
- Matrices ≤ 3×3 compute quickly
- 4×4 may be slow for fully symbolic (use structure!)
- Add structure when possible (zeros, blocks)

---

### Workflow 2: Block-Diagonal Matrices

**When to use**: Your matrix has block structure with zeros off the blocks.

**Steps**:
```julia
# Large matrix with block structure
@variables a b c d e f

M = [a  b  0  0  0  0;
     b  a  0  0  0  0;
     0  0  c  d  0  0;
     0  0  d  c  0  0;
     0  0  0  0  e  f;
     0  0  0  0  f  e]

# Automatic block detection
λ = eigvals(M)
# Solves three 2×2 blocks independently
```

**Why this is fast**: 6×6 would normally require degree-6 polynomial (impossible), but block structure reduces to three quadratics.

**How to recognize**: Look for zero-blocks separating independent subsystems.

---

### Workflow 3: Circulant Matrices

**When to use**: Each row is a cyclic shift of the previous row.

**Steps**:
```julia
@variables a b c d

# 4×4 circulant
C = [a b c d;
     d a b c;
     c d a b;
     b c d a]

# DFT-based closed form (works for ANY size!)
λ = eigvals(C)
```

**Pattern recognition**:
```
Row 1: [a b c d]
Row 2: [d a b c]  ← shift right by 1
Row 3: [c d a b]  ← shift right by 2
Row 4: [b c d a]  ← shift right by 3
```

**Works for huge matrices**: Even 100×100 circulant has closed form!

---

### Workflow 4: Tridiagonal Matrices

**When to use**: Matrix is tridiagonal with constant diagonals.

**Steps**:
```julia
@variables a b

# Symmetric Toeplitz tridiagonal
n = 5
T = diagm(0 => fill(a, n), 
          1 => fill(b, n-1),
         -1 => fill(b, n-1))

# Closed-form via cosines
λ = eigvals(T)
# λₖ = a + 2b·cos(kπ/(n+1)) for k=1,...,n
```

**Pattern recognition**:
- All diagonal entries the same: a
- All super-diagonal the same: b
- All sub-diagonal the same: b
- Symmetric

---

### Workflow 5: Kronecker Products

**When to use**: Matrix is a tensor product A ⊗ B.

**Steps**:
```julia
@variables a b c d

# Two small matrices
A = [a 0; 0 b]  # 2×2
B = [c 0; 0 d]  # 2×2

# Kronecker product (4×4)
K = kron(A, B)

# Eigenvalues are products
λ = eigvals(K)
# [ac, ad, bc, bd]
```

**Why this works**: 4×4 is manageable, but what if A and B were larger?

**Example**: 6×6 = (2×2) ⊗ (3×3)
```julia
A = [a b; b a]          # 2×2, eigenvalues {a+b, a-b}
B = [c d 0; d c d; 0 d c]  # 3×3, eigenvalues computed via cubic

K = kron(A, B)  # 6×6 matrix

λ = eigvals(K)
# Computes eigenvalues of A and B separately
# Then forms all products λᵢ(A) · λⱼ(B)
```

Without Kronecker structure, degree-6 polynomial is impossible!

---

## Working with Special Patterns

### Hermitian Matrices

```julia
@variables a::Real b::Real

# Hermitian matrix (complex conjugate symmetric)
H = [a b; conj(b) a]

# Hint that it's Hermitian for optimizations
λ = eigvals(H, structure=:hermitian)
```

**Properties**:
- Always real eigenvalues
- Orthogonal eigenvectors
- Diagonalizable

---

### Permutation Matrices

```julia
# 5×5 permutation: (1→2→3→1), (4↔5)
P = [0 1 0 0 0;   # 1→2
     0 0 1 0 0;   # 2→3
     1 0 0 0 0;   # 3→1
     0 0 0 0 1;   # 4→5
     0 0 0 1 0]   # 5→4

λ = eigvals(P)
# 3-cycle: {1, ω, ω²} where ω = exp(2πi/3)
# 2-cycle: {1, -1}
# Total: {1, 1, 1, -1, ω, ω²}
```

**Key insight**: All eigenvalues are roots of unity (magnitude 1).

---

### Anti-Diagonal Matrices

```julia
@variables a b c

# 5×5 anti-diagonal
A = [0 0 0 0 a;
     0 0 0 b 0;
     0 0 c 0 0;
     0 b 0 0 0;
     a 0 0 0 0]

λ = eigvals(A)
# {c, ±a, ±b}
```

**Pattern**: Eigenvalues come in ±pairs (except center for odd n).

---

## Performance Optimization

### Tip 1: Use Structure Hints

```julia
# Without hint (may be slower)
λ₁ = eigvals(M)

# With hint (faster)
λ₂ = eigvals(M, structure=:hermitian)
```

**When to use**:
- `:hermitian` - Hermitian matrices
- `:symmetric` - Real symmetric
- `:unitary` - Unitary matrices
- `:general` - No assumptions

---

### Tip 2: Eigenvalues Only

```julia
# Slow: computes eigenvectors
E = eigen(M)
λ = E.values

# Fast: eigenvalues only
λ = eigvals(M)
```

**Speedup**: 2-5x for medium matrices, more for large ones.

---

### Tip 3: Partial Substitution

If your matrix has many symbolic variables but some can be numeric:

```julia
@variables a b c d e f g h

# Fully symbolic 4×4 (slow, huge expressions)
M = [a b c d; e f g h; ...]

# Substitute some variables
using Symbolics: substitute
M_partial = substitute(M, Dict(
    c => 0, d => 0,  # Add structure with zeros
    e => 0, g => 0
))

# Now much faster
λ = eigvals(M_partial)
```

---

### Tip 4: Increase Timeout for Complex Matrices

```julia
# Default timeout: 300 seconds
λ = eigvals(M)

# Allow more time
λ = eigvals(M, timeout=600)
```

---

### Tip 5: Simplify Results

```julia
using Symbolics: simplify

λ = eigvals(M)

# Simplify expressions
λ_simple = simplify.(λ)
```

---

## Troubleshooting

### Problem: "Cannot solve degree n ≥ 5 without structure"

**Cause**: Matrix is 5×5+ with no detected structure.

**Solutions**:

1. **Check for block structure**:
```julia
# Add zeros to create blocks
M_blocked = [A zeros(2,2); zeros(3,2) B]
```

2. **Look for special patterns**: Is it circulant? Tridiagonal? Kronecker product?

3. **Use numeric substitution**:
```julia
# Substitute some variables with numbers
M_numeric = substitute(M, Dict(a => 1, b => 2))
```

4. **Use numerical methods**:
```julia
# Convert to Float64 and solve numerically
M_float = Float64.(substitute(M, Dict(...)))
λ_numeric = eigvals(M_float)  # Standard LinearAlgebra
```

---

### Problem: Expression Complexity Error

**Cause**: Expressions exceed `max_terms` limit (default 10,000).

**Solutions**:

1. **Increase limit** (watch memory!):
```julia
λ = eigvals(M, max_terms=50000)
```

2. **Add structure**:
```julia
# Replace off-block elements with zeros
M_sparse = ...
```

3. **Use expand=false**:
```julia
λ = eigvals(M, expand=false)
```

---

### Problem: Computation Timeout

**Cause**: Exceeded time limit (default 300s).

**Solutions**:

1. **Increase timeout**:
```julia
λ = eigvals(M, timeout=900)  # 15 minutes
```

2. **Simplify matrix**:
```julia
# Use fewer variables
# Add structure (zeros, blocks)
# Partial numeric substitution
```

3. **Check matrix size**: 4×4 fully symbolic takes ~minutes. 5×5+ may be impossible.

---

### Problem: Wrong Eigenvalues

**Cause**: Incorrect structure hint or numerical issues.

**Solutions**:

1. **Remove structure hint**:
```julia
λ = eigvals(M)  # Auto-detect
```

2. **Verify numerically**:
```julia
# Substitute numeric values
M_test = substitute(M, Dict(a=>1, b=>2, c=>3))
λ_numeric = eigvals(M_test)

# Check against symbolic
λ_sym = eigvals(M)
λ_sym_test = substitute.(λ_sym, Ref(Dict(a=>1, b=>2, c=>3)))

# Compare
@assert isapprox(sort(λ_numeric), sort(λ_sym_test))
```

---

### Problem: Expression Too Large (Memory Issues)

**Cause**: Quartic formula on fully symbolic 4×4 produces ~13.5 MB per eigenvalue.

**Solutions**:

1. **Don't use fully symbolic 4×4**: Add structure!

2. **Partial substitution**:
```julia
M_partial = substitute(M, Dict(
    # Set some entries to simple values
    M[1,3] => 0,
    M[1,4] => 0,
    M[2,4] => 0
))
```

3. **Use block structure**:
```julia
# Instead of 4×4, use 2×2 blocks
M = [A zeros(2,2); zeros(2,2) B]
```

---

## Advanced Usage

### Custom Eigenvalue Variables

```julia
@variables λ a b

M = [a b; b a]

# Use specific variable name
vals, poly, λ_var = symbolic_eigenvalues(M, var=λ)
```

---

### Characteristic Polynomial Only

```julia
M = [a b; b a]

poly, coeffs, λ = characteristic_polynomial(M)
# poly: λ² - 2a·λ + (a² - b²)
# coeffs: [a² - b², -2a, 1]
```

---

### Eigenpairs with Multiplicity

```julia
M = [a 0 0; 0 a 1; 0 0 a]  # Jordan block

pairs = symbolic_eigenpairs(M)
# [(a, [v₁, v₂, ...])]  # May have multiple eigenvectors
```

---

### Verification

```julia
P, D, pairs = symbolic_diagonalize(M)

# Verify: M = P * D * P⁻¹
using Symbolics: simplify
err = simplify.(M - P * D * inv(P))
@assert all(iszero, err)
```

---

### Mixed Numeric-Symbolic

```julia
# Some entries numeric, some symbolic
@variables a b

M = [a   1.5  0   0;
     1.5 b    0   0;
     0   0    2.0 1.0;
     0   0    1.0 3.0]

λ = eigvals(M)  # Works fine!
```

---

## Practical Examples

### Example 1: Vibrating String

Discrete model of vibrating string with n masses:

```julia
@variables k m ω

n = 5
# Stiffness matrix (tridiagonal)
K = diagm(0 => fill(2k, n),
          1 => fill(-k, n-1),
         -1 => fill(-k, n-1))

# Eigenfrequencies: ω² = λ/m
λ = eigvals(K)
ω_squared = λ ./ m

# ωₖ = √(λₖ/m) = √(k/m) · √(4sin²(kπ/(2(n+1))))
```

---

### Example 2: Quantum Tight-Binding Model

1D tight-binding Hamiltonian:

```julia
@variables t ε

n = 4
H = diagm(0 => fill(ε, n),
          1 => fill(-t, n-1),
         -1 => fill(-t, n-1))

# Energy eigenvalues
E = eigvals(H)
# Eₖ = ε - 2t·cos(kπ/(n+1))
```

---

### Example 3: Graph Laplacian (Cycle)

Laplacian of cycle graph Cₙ:

```julia
@variables d

n = 6
# Cycle graph Laplacian (circulant!)
L = [2d -d  0   0  0  -d;
     -d 2d -d   0  0   0;
      0 -d 2d  -d  0   0;
      0  0 -d  2d -d   0;
      0  0  0  -d 2d  -d;
     -d  0  0   0 -d  2d]

# Eigenvalues via circulant formula
λ = eigvals(L)
# λₖ = 2d(1 - cos(2πk/n)) for k=0,...,n-1
```

---

### Example 4: Coupled Oscillators

Two coupled 2D oscillators:

```julia
@variables k₁ k₂ κ m

# Each oscillator: 2×2 block
A = [k₁ 0; 0 k₁] / m
B = [k₂ 0; 0 k₂] / m

# Coupling
C = [κ 0; 0 κ] / m

# Full system (4×4)
M = [A C; C B]

# If κ=0: decoupled (block-diagonal)
# If κ≠0: analyze coupled system

λ = eigvals(M)
```

---

## Next Steps

- **[API Reference](api_reference.md)** - Complete function documentation
- **[Pattern Library](pattern_library.md)** - All special patterns
- **[Mathematical Background](mathematical_background.md)** - Theory
- **[Implementation](implementation.md)** - Algorithm details
- **[Contributing](contributing.md)** - How to add new patterns

---

## Getting Help

1. **Check documentation**: API reference, pattern library
2. **Verify numerically**: Test with numeric values first
3. **Simplify problem**: Remove variables, add structure
4. **Report issues**: GitHub issues with minimal reproducible example
