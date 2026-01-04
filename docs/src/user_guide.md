# User Guide

Practical guide to using SymbolicDiagonalization.jl effectively.

## Installation

```julia
using Pkg
Pkg.add("SymbolicDiagonalization")
```

Or from the repository:

```bash
cd SymbolicDiagonalization.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Basic Usage

```julia
using Symbolics
using SymbolicDiagonalization

@variables a b c d

# Create a symbolic matrix
M = [a b; b a]

# Eigenvalues only (faster)
eigvals(M)  # [a+b, a-b]

# Eigenvalues + eigenvectors
E = eigen(M)
E.values    # eigenvalues
E.vectors   # eigenvectors as columns
```

## Common Workflows

### Small Matrices (≤ 4×4)

Direct closed-form solution via root formulas:

```julia
@variables a b c d
M = [a b 0 0; b a 0 0; 0 0 c d; 0 0 d c]

eigvals(M)  # [a+b, a-b, c+d, c-d]
```

### Block-Diagonal Matrices

Automatically detected and solved recursively:

```julia
@variables a b c d e f
M = [a b 0 0 0 0;
     b a 0 0 0 0;
     0 0 c d 0 0;
     0 0 d c 0 0;
     0 0 0 0 e f;
     0 0 0 0 f e]

eigvals(M)  # Solves three 2×2 blocks
```

### Circulant Matrices

DFT-based solution for any size:

```julia
@variables a b c d
C = [a b c d; d a b c; c d a b; b c d a]

eigvals(C)  # Works for 100×100 too
```

### Tridiagonal Matrices

Closed-form cosine formula:

```julia
@variables a b
n = 5
T = diagm(0 => fill(a, n), 1 => fill(b, n-1), -1 => fill(b, n-1))

eigvals(T)  # λₖ = a + 2b·cos(kπ/(n+1))
```

### Kronecker Products

Products of factor eigenvalues:

```julia
@variables a b c d
A = [a 0; 0 b]
B = [c 0; 0 d]
K = kron(A, B)  # 4×4

eigvals(K)  # [ac, ad, bc, bd]
```

### SO(2) Rotation Kronecker Products

Kronecker products of SO(2) rotation matrices are automatically detected and produce
clean trigonometric eigenvalues using angle addition formulas:

```julia
@variables θ φ
R_θ = SO2_rotation(θ)
R_φ = SO2_rotation(φ)

K = kron(R_θ, R_φ)
eigvals(K)
# [cos(θ + φ) + im*sin(θ + φ),   # e^{i(θ+φ)}
#  cos(θ - φ) + im*sin(θ - φ),   # e^{i(θ-φ)}
#  cos(θ - φ) - im*sin(θ - φ),   # e^{-i(θ-φ)}
#  cos(θ + φ) - im*sin(θ + φ)]   # e^{-i(θ+φ)}
```

The same-angle case `SO2_rotation(θ) ⊗ SO2_rotation(θ)` automatically simplifies using double-angle formulas:

```julia
eigvals(kron(SO2_rotation(θ), SO2_rotation(θ)))
# [cos(2θ) + im*sin(2θ),   # e^{2iθ}
#  1,                       # θ - θ = 0
#  1,                       # -(θ - θ) = 0
#  cos(2θ) - im*sin(2θ)]   # e^{-2iθ}
```

This works recursively for nested SO(2) Kronecker products.

### SU(2) Kronecker Products

SU(2) is the group of 2×2 unitary matrices with determinant 1. In quantum mechanics, SU(2) describes spin-1/2 systems. The package provides Pauli matrices and SU(2) rotation constructors:

```julia
@variables α β

# Pauli matrices
pauli_x()  # [0 1; 1 0]
pauli_y()  # [0 -im; im 0]
pauli_z()  # [1 0; 0 -1]

# SU(2) rotation matrices (spin-1/2 representation)
SU2_Ux(θ)  # exp(-i θ σx/2)
SU2_Uy(θ)  # exp(-i θ σy/2) - same form as SO(2) but with half-angle
SU2_Uz(θ)  # diagonal: [e^{-iθ/2}, e^{iθ/2}]

# SU(2) ⊗ SU(2) Kronecker product
K = SU2_kron([α, β])  # 4×4 matrix
eigvals(K)
# [cos((α+β)/2) + im*sin((α+β)/2),   # e^{i(α+β)/2}
#  cos((α-β)/2) + im*sin((α-β)/2),   # e^{i(α-β)/2}
#  cos((α-β)/2) - im*sin((α-β)/2),   # e^{-i(α-β)/2}
#  cos((α+β)/2) - im*sin((α+β)/2)]   # e^{-i(α+β)/2}
```

The key difference from SO(2) is the half-angle: SU(2) uses `θ/2` instead of `θ`.

### SU(3) Kronecker Products

SU(3) is the group of 3×3 unitary matrices with determinant 1. It's fundamental in particle physics (color charge in QCD). The package provides Gell-Mann matrices and diagonal SU(3) constructors:

```julia
@variables α₁ α₂ β₁ β₂

# Gell-Mann matrices (8 generators of SU(3))
gellmann_matrices()  # Returns [λ1, λ2, ..., λ8]

# Diagonal SU(3) matrix (Cartan subalgebra)
# U = diag(e^{iα₁}, e^{iα₂}, e^{-i(α₁+α₂)})  # det = 1
U = SU3_diagonal_trig(α₁, α₂)

# SU(3) ⊗ SU(3) Kronecker product (9×9 diagonal matrix)
K = SU3_kron((α₁, α₂), (β₁, β₂))
eigvals(K)
# 9 eigenvalues of the form cos(θ) + im*sin(θ) where θ is:
# α₁+β₁, α₁+β₂, α₁-(β₁+β₂), α₂+β₁, α₂+β₂, α₂-(β₁+β₂),
# -(α₁+α₂)+β₁, -(α₁+α₂)+β₂, -(α₁+α₂)-(β₁+β₂)
```

Direct eigenvalue computation without building the matrix:

```julia
SU3_kron_eigenvalues((α₁, α₂), (β₁, β₂))  # Returns 9 eigenvalues directly
```

### Simplifying Eigenvalue Expressions

The package provides aggressive simplification functions for cleaning up symbolic eigenvalue expressions, particularly those involving trigonometric functions.

#### Key Simplification Functions

```julia
using SymbolicDiagonalization: aggressive_simplify, simplify_eigenvalue, simplify_eigenvalues

@variables θ

# aggressive_simplify: Combines trig, algebraic, and sqrt-trig rules
aggressive_simplify(sqrt(1 - cos(θ)^2))  # → sin(θ)
aggressive_simplify(sqrt(1 - sin(θ)^2))  # → cos(θ)
aggressive_simplify(sin(θ)^2 + cos(θ)^2) # → 1

# simplify_eigenvalue: Specialized for eigenvalue expressions
# Handles complex parts separately for clean output
expr = cos(θ) + im*sqrt(1 - cos(θ)^2)
simplify_eigenvalue(expr)  # → cos(θ) + im*sin(θ)

# simplify_eigenvalues: Apply to a vector of eigenvalues
vals = [cos(θ) + im*sqrt(1 - cos(θ)^2), cos(θ) - im*sqrt(1 - cos(θ)^2)]
simplify_eigenvalues(vals)  # → [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
```

#### When to Use

- **Automatic**: Block-diagonal SO(4), SO(3), and SO(2) matrices automatically use these simplifications
- **Manual**: When working with custom matrices that produce `sqrt(1 - cos²θ)` or similar expressions

#### Simplification Rules Applied

| Rule | Before | After |
|------|--------|-------|
| Pythagorean | `sin²θ + cos²θ` | `1` |
| Sqrt-trig | `sqrt(1 - cos²θ)` | `sin(θ)` |
| Sqrt-trig | `sqrt(1 - sin²θ)` | `cos(θ)` |
| Algebraic | `sqrt(x)^2` | `x` |
| Algebraic | `sqrt(a)*sqrt(b)` | `sqrt(a*b)` |

### Nested Kronecker Products

Arbitrary-depth Kronecker products are handled recursively. The eigenvalues of
$A_1 \otimes A_2 \otimes \cdots \otimes A_n$ are all products:
$$\lambda_{i_1, i_2, \ldots, i_n} = \lambda_{i_1}(A_1) \cdot \lambda_{i_2}(A_2) \cdots \lambda_{i_n}(A_n)$$

```julia
using SymbolicDiagonalization, Symbolics

# Create 5 independent 2×2 symmetric matrices (15 parameters total)
matrices = Matrix{Num}[]
for i in 1:5
    a = Symbolics.variable(Symbol("a$i"))
    b = Symbolics.variable(Symbol("b$i"))
    c = Symbolics.variable(Symbol("c$i"))
    push!(matrices, [a b; b c])
end

# Build 32×32 Kronecker product
K = reduce(kron, matrices)

# Solve symbolically - returns 32 eigenvalues
vals, vecs, pairs = symbolic_eigenvalues(K; timeout=120)
```

This scales efficiently to very large matrices:

| Depth | Matrix Size | Parameters | Time |
|-------|-------------|------------|------|
| 5 | 32×32 | 15 | ~12s |
| 7 | 128×128 | 21 | ~23s |
| 10 | 1024×1024 | 30 | ~33s |

### Full 6-Parameter 3×3 Symmetric Matrices

The package handles the most general 3×3 symmetric matrix with 6 independent parameters:

```julia
@variables a b c d e f
A = [a b c; b d e; c e f]

# Uses diagonal shift optimization internally
vals, vecs, pairs = symbolic_eigenvalues(A; timeout=600)
```

This employs a diagonal shift optimization: shift by $f \cdot I$, solve a 5-variable
problem, then back-substitute. The resulting eigenvalue expressions are large but exact.

## Performance Tips

| Tip | Example |
|-----|---------|
| Use structure hints | `eigvals(M, structure=:hermitian)` |
| Eigenvalues only | `eigvals(M)` vs `eigen(M)` |
| Partial substitution | `substitute(M, Dict(c => 0))` |
| Increase timeout | `eigvals(M, timeout=600)` |
| Disable expansion | `eigvals(M, expand=false)` |

## Troubleshooting

### "Cannot solve degree n ≥ 5 without structure"

Matrix is 5×5+ with no detected structure. Solutions:
1. Add zeros to create block structure
2. Check for circulant, tridiagonal, or Kronecker patterns
3. Substitute numeric values for some variables
4. Use numerical methods via `LinearAlgebra.eigvals` on Float64 matrix

### Expression Complexity Error

Expressions exceed term limit. Solutions:
1. Increase limit: `eigvals(M, max_terms=50000)`
2. Add structure (zeros, blocks)
3. Disable expansion: `eigvals(M, expand=false)`

### Computation Timeout

Exceeded time limit. Solutions:
1. Increase timeout: `eigvals(M, timeout=900)`
2. Simplify matrix (fewer variables, more structure)
3. Use partial numeric substitution

### Memory Issues

Fully symbolic 4×4 quartic produces very large expressions. Solutions:
1. Add structure to avoid quartic formula
2. Use block decomposition
3. Partial numeric substitution

## Verification

```julia
# Verify diagonalization
P, D, pairs = symbolic_diagonalize(M)
using Symbolics: simplify
@assert all(iszero, simplify.(M - P * D * inv(P)))

# Verify numerically
M_test = substitute(M, Dict(a=>1, b=>2))
λ_sym = substitute.(eigvals(M), Ref(Dict(a=>1, b=>2)))
λ_num = LinearAlgebra.eigvals(Float64.(M_test))
@assert isapprox(sort(real.(λ_sym)), sort(λ_num))
```

### Robust Numerical Evaluation

For complex eigenvalue expressions (e.g., 6-parameter matrices), use the internal
`_evaluate_symbolic_expr` function which handles floating-point edge cases:

```julia
using SymbolicDiagonalization, Symbolics, LinearAlgebra

@variables a b c d e f
A = [a b c; b d e; c e f]
vals, _, _ = symbolic_eigenvalues(A; timeout=600)

# Test values
test_vals = Dict(a => 1.0, b => 0.5, c => 0.3, d => 2.0, e => 0.4, f => 3.0)

# Robust evaluation (handles sqrt of tiny negative numbers)
computed = [real(SymbolicDiagonalization._evaluate_symbolic_expr(v, test_vals)) 
            for v in vals]

# Compare to numerical eigenvalues
A_num = Float64[1.0 0.5 0.3; 0.5 2.0 0.4; 0.3 0.4 3.0]
true_eigs = eigvals(Symmetric(A_num))

@assert isapprox(sort(computed), sort(true_eigs), atol=1e-10)
```

This function uses `Complex{Float64}` arithmetic internally to avoid `NaN` from
`sqrt` of tiny negative numbers caused by floating-point precision issues.

## Physical Examples

### Vibrating String

```julia
@variables k m
n = 5
K = diagm(0 => fill(2k, n), 1 => fill(-k, n-1), -1 => fill(-k, n-1))
λ = eigvals(K)  # ωₖ² = λₖ/m
```

### Tight-Binding Model

```julia
@variables t ε
n = 4
H = diagm(0 => fill(ε, n), 1 => fill(-t, n-1), -1 => fill(-t, n-1))
E = eigvals(H)  # Eₖ = ε - 2t·cos(kπ/(n+1))
```

### Cycle Graph Laplacian

```julia
@variables d
L = [2d -d 0 0 0 -d; -d 2d -d 0 0 0; 0 -d 2d -d 0 0;
     0 0 -d 2d -d 0; 0 0 0 -d 2d -d; -d 0 0 0 -d 2d]
λ = eigvals(L)  # λₖ = 2d(1 - cos(2πk/n))
```

## See Also

- [API Reference](api_reference.md) - Complete function signatures
- [Pattern Library](pattern_library.md) - All 16+ patterns with details
- [Mathematical Background](mathematical_background.md) - Theory
