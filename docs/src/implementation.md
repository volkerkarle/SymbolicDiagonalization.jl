# Implementation Details

This document provides technical details about the algorithms and implementation strategies used in SymbolicDiagonalization.jl.

## Table of Contents

1. [Overview](#overview)
2. [Characteristic Polynomial Computation](#characteristic-polynomial-computation)
3. [Root Solvers (Degrees 1-4)](#root-solvers-degrees-1-4)
4. [Structure Detection](#structure-detection)
5. [Special Pattern Solvers](#special-pattern-solvers)
6. [Eigenvector Computation](#eigenvector-computation)
7. [Expression Management](#expression-management)
8. [Performance Considerations](#performance-considerations)

---

## Overview

SymbolicDiagonalization.jl uses a multi-layered approach to compute eigenvalues and eigenvectors:

1. **Structure Detection**: Identify exploitable matrix patterns
2. **Pattern-Specific Solvers**: Use specialized algorithms for detected patterns
3. **Characteristic Polynomial**: Compute det(λI - A) using Bareiss algorithm
4. **Root Finding**: Solve polynomials up to degree 4 using closed-form formulas
5. **Eigenvector Computation**: Use RREF-based nullspace or adjugate method

The implementation prioritizes **symbolic correctness** over numeric efficiency, using fraction-free algorithms and avoiding floating-point operations wherever possible.

---

## Characteristic Polynomial Computation

### Algorithm: Bareiss Determinant

The characteristic polynomial det(λI - A) is computed using the **Bareiss algorithm**, a fraction-free variant of Gaussian elimination.

**File**: `src/charpoly.jl`

### Why Bareiss?

Standard Gaussian elimination requires division, which causes intermediate expression blowup in symbolic computation:
- **Standard GE**: Expressions grow like O(2^n) due to nested divisions
- **Bareiss**: Keeps expressions polynomial-sized by avoiding division until the end

### Algorithm Description

Given matrix M, Bareiss computes det(M) without fractions:

```julia
function _bareiss_det(M)
    n = size(M, 1)
    A = Matrix{eltype(M)}(M)
    prev = one(eltype(M))
    
    for k in 1:n-1
        pivot = A[k, k]
        _issymzero(pivot) && error("zero pivot encountered")
        
        # Update submatrix using division-free formula
        for i in k+1:n, j in k+1:n
            A[i, j] = (A[i, j] * pivot - A[i, k] * A[k, j]) / prev
        end
        
        # Zero out column below pivot
        fill!(view(A, k+1:n, k), zero(eltype(M)))
        prev = pivot
    end
    
    return A[n, n]
end
```

**Key insight**: The formula `(A[i,j] * pivot - A[i,k] * A[k,j]) / prev` is always divisible by `prev` (the previous pivot), so no fractions appear.

### Complexity

- **Time**: O(n³) symbolic operations
- **Expression size**: Polynomial growth (much better than standard methods)
- **Space**: O(n²) for the matrix

### Coefficient Extraction

Instead of polynomial division (which is brittle for symbolic expressions), we use **derivative-based extraction**:

```julia
function _collect_coefficients(poly, λ, degree)
    coeffs = Vector{Any}(undef, degree + 1)
    coeffs[1] = substitute(poly, λ => 0)  # Constant term
    
    deriv = poly
    fact = 1
    for k in 1:degree
        deriv = derivative(deriv, λ)
        fact *= k
        coeffs[k + 1] = substitute(deriv, λ => 0) / fact
    end
    
    return coeffs
end
```

This uses the Taylor series expansion: `p(λ) = Σ (p^(k)(0) / k!) * λ^k`

---

## Root Solvers (Degrees 1-4)

### Overview

SymbolicDiagonalization.jl implements closed-form root solvers for polynomials up to degree 4:

| Degree | Method | Formula Name |
|--------|--------|--------------|
| 1 | Linear | Direct |
| 2 | Quadratic | Quadratic formula |
| 3 | Cubic | Cardano's method |
| 4 | Quartic | Ferrari's method |

**File**: `src/roots.jl`

### Linear Solver (Degree 1)

For `ax + b = 0`:

```julia
_roots_linear(c) = [-c[1] / c[2]]  # c = [b, a]
```

**Complexity**: O(1) symbolic operations

### Quadratic Solver (Degree 2)

For `ax² + bx + c = 0`:

```julia
function _roots_quadratic(c)
    a, b, d = c[3], c[2], c[1]  # d is constant term
    disc = b^2 - 4a*d
    disc = _aggressive_simplify(disc)
    sqrt_disc = sqrt(disc)
    return [(-b - sqrt_disc) / (2a), (-b + sqrt_disc) / (2a)]
end
```

**Key features**:
- Discriminant simplification to reduce expression size
- Symbolic square root handling

**Complexity**: O(n²) where n is the expression size

### Cubic Solver (Degree 3)

For `ax³ + bx² + cx + d = 0`, we use **Cardano's method**:

**Step 1**: Depress the cubic (eliminate x² term)

Substitution `x = y - b/(3a)` transforms to: `y³ + py + q = 0`

where:
- `p = c/a - b²/(3a²)`
- `q = 2b³/(27a³) - bc/(3a²) + d/a`

**Step 2**: Compute discriminant

`Δ = (q/2)² + (p/3)³`

**Step 3**: Compute cube roots

```julia
C = cbrt(-q/2 + √Δ)
D = cbrt(-q/2 - √Δ)
```

**Step 4**: Construct roots

Let `ω = exp(2πi/3) = -1/2 + (√3/2)i` be a primitive cube root of unity.

```julia
roots = [
    C + D,
    ω*C + ω²*D,
    ω²*C + ω*D
]
```

**Step 5**: Shift back

`x = y - b/(3a)`

**Complexity**: O(n³) where n is the expression size

**Implementation notes**:
- Aggressive simplification applied to `p`, `q`, and `Δ`
- Complex arithmetic handled symbolically
- Numeric cube roots used with symbolic coefficients

### Quartic Solver (Degree 4)

For `ax⁴ + bx³ + cx² + dx + e = 0`, we use **Ferrari's method** via a resolvent cubic.

**Step 1**: Depress the quartic (eliminate x³ term)

Substitution `x = y - b/(4a)` transforms to: `y⁴ + py² + qy + r = 0`

where:
- `p = c/a - 3b²/(8a²)`
- `q = d/a + b³/(8a³) - bc/(2a²)`
- `r = e/a - 3b⁴/(256a⁴) + b²c/(16a³) - bd/(4a²)`

**Step 2**: Solve resolvent cubic

The resolvent cubic is: `8z³ - 4pz² - 8rz + (4pr - q²) = 0`

Let α be any root of this cubic (we use the first one).

**Step 3**: Decompose into two quadratics

The depressed quartic factors as:

`(y² + β*y + (α + γ))(y² - β*y + (α - γ)) = 0`

where:
- `β² = 2α - p`
- `γ = -q/(2β)`

**Step 4**: Solve quadratics

Each quadratic gives two roots using the quadratic formula.

**Step 5**: Shift back

`x = y - b/(4a)`

**Complexity**: O(n⁴) where n is the expression size

**Implementation details**:
```julia
function _roots_quartic(c)
    # Normalize and extract coefficients
    a, b, cc, d, e = c[5], c[4], c[3], c[2], c[1]
    
    # Depress to y⁴ + py² + qy + r = 0
    p = cc/a - 3b²/(8a²)
    q = d/a + b³/(8a³) - b*cc/(2a²)
    r = e/a - 3b⁴/(256a⁴) + b²*cc/(16a³) - b*d/(4a²)
    
    # Simplify intermediate expressions
    p = _aggressive_simplify(p)
    q = _aggressive_simplify(q)
    r = _aggressive_simplify(r)
    
    # Resolvent cubic coefficients
    resolvent = [4p*r - q², -8r, -4p, 8]
    alphas = _roots_cubic(resolvent)
    alpha = alphas[1]  # Choose first root
    
    # Compute β and γ
    beta_sq = 2alpha - p
    beta_sq = _aggressive_simplify(beta_sq)
    beta = _symbolic_sqrt(beta_sq)
    gamma = -q/(2beta)
    gamma = _aggressive_simplify(gamma)
    
    # Two quadratics: y² ± β*y + (α ± γ)
    t1 = beta² - 4(alpha + gamma)
    t2 = beta² - 4(alpha - gamma)
    t1 = _aggressive_simplify(t1)
    t2 = _aggressive_simplify(t2)
    
    # Solve both quadratics
    roots_y = [
        (-beta - _symbolic_sqrt(t1)) / 2,
        (-beta + _symbolic_sqrt(t1)) / 2,
        (beta - _symbolic_sqrt(t2)) / 2,
        (beta + _symbolic_sqrt(t2)) / 2
    ]
    
    # Shift back
    return roots_y .- b/(4a)
end
```

**Warning**: Quartic formulas produce **extremely large expressions** for fully symbolic 4×4 matrices (~13.5 MB per eigenvalue).

---

## Structure Detection

### Overview

Before attempting expensive polynomial root finding, we try to detect exploitable matrix structure.

**File**: `src/diagonalize.jl` (functions starting with `_is_*` and `_detect_*`)

### Detection Hierarchy

```
Structure Detection
├── Diagonal (_is_diagonal)
├── Triangular (_is_triangular)
├── Block-Diagonal (_detect_multiple_blocks)
├── Persymmetric (_is_persymmetric)
├── Circulant (_is_circulant)
├── Block Circulant (_is_block_circulant)
├── Toeplitz Tridiagonal (_is_toeplitz_tridiagonal)
├── Anti-Diagonal (_is_antidiagonal)
├── Permutation (_is_permutation_matrix)
├── Kronecker Product (_is_kronecker_product)
└── Special 5×5 (_detect_special_5x5_tridiagonal)
```

### Implementation Details

#### Symbolic Zero Detection

All structure detection relies on checking if expressions are zero:

```julia
function _issymzero(x)
    # Try standard iszero
    try
        v = Base.iszero(x)
        if v isa Bool
            return v
        end
    catch
    end
    
    # Try simplify then check
    try
        sx = Symbolics.simplify(x)
        v = Symbolics.iszero(sx)
        return v === true
    catch
    end
    
    # Try direct Symbolics.iszero
    try
        v = Symbolics.iszero(x)
        return v === true
    catch
    end
    
    # Conservative: not proven zero
    return false
end
```

**Key insight**: We return `false` if we can't prove it's zero (conservative approach).

#### Diagonal Detection

```julia
function _is_diagonal(mat)
    m, n = size(mat)
    for i in 1:m, j in 1:n
        i == j && continue
        !_issymzero(mat[i, j]) && return false
    end
    return true
end
```

**Complexity**: O(n²) symbolic zero checks

#### Block-Diagonal Detection

Finds all block-diagonal structure using a greedy algorithm:

```julia
function _detect_multiple_blocks(mat)
    n = size(mat, 1)
    blocks = Tuple{Int,Int}[]
    current_start = 1
    
    while current_start <= n
        # Find smallest block starting at current_start
        block_end = n  # Default: extend to end
        
        for k in current_start:(n-1)
            # Check if [current_start:k, k+1:n] and [k+1:n, current_start:k] are zero
            upper_right_zero = all(_issymzero, mat[i, j] 
                for i in current_start:k, j in (k+1):n)
            lower_left_zero = all(_issymzero, mat[i, j] 
                for i in (k+1):n, j in current_start:k)
            
            if upper_right_zero && lower_left_zero
                block_end = k
                break
            end
        end
        
        push!(blocks, (current_start, block_end))
        current_start = block_end + 1
    end
    
    # Only return if we found actual blocks (more than just the whole matrix)
    return length(blocks) > 1 ? blocks : nothing
end
```

**Complexity**: O(n⁴) symbolic operations in worst case

#### Circulant Detection

```julia
function _is_circulant(mat)
    n = size(mat, 1)
    first_row = mat[1, :]
    
    # Check each row is a cyclic shift of the first row
    for i in 2:n
        for j in 1:n
            expected_idx = mod1(j - (i - 1), n)
            if !_issymzero(mat[i, j] - first_row[expected_idx])
                return false
            end
        end
    end
    
    return true
end
```

**Complexity**: O(n²) symbolic zero checks

#### Block Circulant Detection

Tries different block sizes k where n is divisible by k:

```julia
function _is_block_circulant(mat)
    m, n = size(mat)
    m == n || return nothing
    
    # Try different block sizes
    for k in 2:div(m, 2)
        m % k == 0 || continue
        n_blocks = div(m, k)
        
        # Extract blocks from first block row
        blocks = [mat[1:k, (i-1)*k+1:i*k] for i in 1:n_blocks]
        
        # Check if each block row is a cyclic shift
        is_block_circulant = true
        for block_row in 2:n_blocks
            for block_col in 1:n_blocks
                shift_idx = mod1(block_col - (block_row - 1), n_blocks)
                expected_block = blocks[shift_idx]
                actual_block = mat[(block_row-1)*k+1:block_row*k, 
                                  (block_col-1)*k+1:block_col*k]
                
                # Check if blocks match
                if !all(_issymzero(actual_block[i, j] - expected_block[i, j]) 
                       for i in 1:k, j in 1:k)
                    is_block_circulant = false
                    break
                end
            end
            is_block_circulant || break
        end
        
        if is_block_circulant
            return (n_blocks, k, blocks)
        end
    end
    
    return nothing
end
```

**Complexity**: O(n⁴) in worst case (trying all block sizes)

#### Kronecker Product Detection

Tries different factorizations N = m × n:

```julia
function _is_kronecker_product(mat)
    N = size(mat, 1)
    
    # Try factorizations N = m * n
    for m in 2:div(N, 2)
        N % m == 0 || continue
        n = div(N, m)
        
        # Extract candidate B from (1,1) block
        B_candidate = mat[1:n, 1:n]
        
        # Try to extract A by looking at block pattern
        A = zeros(eltype(mat), m, m)
        
        for i_block in 1:m, j_block in 1:m
            block = mat[(i_block-1)*n+1:i_block*n, (j_block-1)*n+1:j_block*n]
            
            # Find scalar such that block = scalar * B_candidate
            scalar = find_scalar_multiple(block, B_candidate)
            
            # Verify block = scalar * B_candidate
            if verify_match(block, scalar, B_candidate)
                A[i_block, j_block] = scalar
            else
                @goto next_factorization
            end
        end
        
        return (A, B_candidate, m, n)
        
        @label next_factorization
    end
    
    return nothing
end
```

**Complexity**: O(n⁵) in worst case

#### Permutation Matrix Detection

```julia
function _is_permutation_matrix(A)
    n = size(A, 1)
    
    # Check each row has exactly one 1 and rest zeros
    for i in 1:n
        count_ones = sum(_issymzero(A[i, j] - 1) for j in 1:n)
        count_ones == 1 || return false
        
        # Check rest are zeros
        for j in 1:n
            !_issymzero(A[i, j] - 1) && !_issymzero(A[i, j]) && return false
        end
    end
    
    # Check each column has exactly one 1
    for j in 1:n
        count_ones = sum(_issymzero(A[i, j] - 1) for i in 1:n)
        count_ones == 1 || return false
    end
    
    return true
end
```

**Complexity**: O(n²) symbolic zero checks

---

## Special Pattern Solvers

### Circulant Matrices

**Theory**: A circulant matrix has eigenvalues given by the DFT of its first row.

For circulant matrix C = circ(c₀, c₁, ..., c_{n-1}):

`λⱼ = Σₖ cₖ ωʲᵏ` for j = 0, 1, ..., n-1

where ω = exp(2πi/n).

```julia
function _circulant_eigenvalues(mat)
    n = size(mat, 1)
    first_row = mat[1, :]
    eigenvalues = Vector{Any}(undef, n)
    
    for j in 0:(n-1)
        λ = first_row[1]  # c₀ term
        
        for k in 1:(n-1)
            θ = 2π * j * k / n
            ω_power = cos(θ) + im * sin(θ)
            λ += first_row[k+1] * ω_power
        end
        
        eigenvalues[j+1] = Symbolics.simplify(λ)
    end
    
    return eigenvalues
end
```

**Complexity**: O(n²) symbolic operations (much better than O(n³) for general case)

### Block Circulant Matrices

**Theory**: Block circulant reduces to n separate k×k eigenvalue problems.

For block circulant with blocks [A₀, A₁, ..., A_{n-1}]:

Eigenvalues = ⋃ⱼ eigvals(Dⱼ)

where `Dⱼ = Σₖ ωʲᵏ Aₖ` and ω = exp(2πi/n).

```julia
function _block_circulant_eigenvalues(mat, n_blocks, block_size, blocks)
    all_eigenvalues = []
    
    for j in 0:(n_blocks-1)
        # Compute Dⱼ = Σₖ ωʲᵏ Aₖ
        D_j = zeros(eltype(mat), block_size, block_size)
        
        for k in 0:(n_blocks-1)
            θ = 2π * j * k / n_blocks
            ω_power = cos(θ) + im * sin(θ)
            D_j .+= ω_power .* blocks[k + 1]
        end
        
        # Recursively solve k×k eigenvalue problem
        D_j = Symbolics.simplify.(D_j)
        vals_j, _, _ = symbolic_eigenvalues(D_j)
        append!(all_eigenvalues, vals_j)
    end
    
    return all_eigenvalues
end
```

**Complexity**: O(n × T(k)) where T(k) is time to solve k×k matrix

### Toeplitz Tridiagonal Matrices

**Theory**: Symmetric Toeplitz tridiagonal has closed-form eigenvalues via orthogonal polynomials.

For matrix with diagonal a and off-diagonals b:

`λₖ = a + 2b cos(kπ/(n+1))` for k = 1, 2, ..., n

```julia
function _toeplitz_tridiagonal_eigenvalues(n, a, b, c)
    eigenvalues = Vector{Any}(undef, n)
    
    for k in 1:n
        θ = k * π / (n + 1)
        λₖ = a + 2 * b * cos(θ)
        eigenvalues[k] = Symbolics.simplify(λₖ)
    end
    
    return eigenvalues
end
```

**Complexity**: O(n) symbolic operations (optimal!)

### Kronecker Products

**Theory**: If A has eigenvalues {λᵢ} and B has eigenvalues {μⱼ}, then A ⊗ B has eigenvalues {λᵢ μⱼ}.

```julia
function _kronecker_eigenvalues(A, B, m, n)
    # Compute eigenvalues of A
    λ_A, _, _ = symbolic_eigenvalues(A)
    
    # Compute eigenvalues of B
    λ_B, _, _ = symbolic_eigenvalues(B)
    
    # Compute all products λᵢ * μⱼ
    eigenvalues = []
    for λ in λ_A, μ in λ_B
        push!(eigenvalues, λ * μ)
    end
    
    return eigenvalues
end
```

**Complexity**: O(T(m) + T(n) + mn) where T(k) is time to solve k×k matrix

### Permutation Matrices

**Theory**: Eigenvalues are roots of unity determined by cycle structure.

For a cycle of length k: eigenvalues are `exp(2πij/k)` for j = 0, 1, ..., k-1

```julia
function _compute_permutation_eigenvalues(A)
    cycles = _permutation_to_cycles(A)
    eigenvalues = []
    
    for cycle_length in cycles
        if cycle_length == 1
            push!(eigenvalues, 1)  # Fixed point
        elseif cycle_length == 2
            push!(eigenvalues, 1, -1)  # Transposition
        else
            # k-th roots of unity
            for j in 0:(cycle_length - 1)
                angle = 2π * j / cycle_length
                push!(eigenvalues, exp(im * angle))
            end
        end
    end
    
    return eigenvalues
end
```

**Complexity**: O(n) to find cycles, O(n) to generate eigenvalues

### Anti-Diagonal Matrices

**Theory**: Symmetric anti-diagonal has eigenvalues in ±pairs.

For anti-diagonal entries [a₁, a₂, ..., aₙ]:
- Odd n: eigenvalues are [a_{mid}, ±a₁, ±a₂, ..., ±a_{mid-1}]
- Even n: eigenvalues are [±a₁, ±a₂, ..., ±a_{n/2}]

```julia
function _antidiagonal_eigenvalues(mat)
    n = size(mat, 1)
    antidiag = [mat[i, n + 1 - i] for i in 1:n]
    
    if n % 2 == 1
        # Odd: one middle eigenvalue, rest in pairs
        mid = (n + 1) ÷ 2
        eigenvalues = [antidiag[mid]]
        for i in 1:mid-1
            push!(eigenvalues, antidiag[i], -antidiag[i])
        end
    else
        # Even: all in pairs
        eigenvalues = []
        for i in 1:n÷2
            push!(eigenvalues, antidiag[i], -antidiag[i])
        end
    end
    
    return eigenvalues
end
```

**Complexity**: O(n) symbolic operations

### Persymmetric Splitting

**Theory**: Symmetric persymmetric matrices (Q[i,j] = Q[n+1-j, n+1-i]) can be split into two half-sized blocks.

```julia
function _persymmetric_split(mat)
    n = size(mat, 1)
    n % 2 != 0 && return nothing  # Only for even n
    
    # Build transformation matrix P = [(I+J)/√2, (I-J)/√2]
    # where J is the exchange matrix (anti-identity)
    
    half = div(n, 2)
    P = zeros(eltype(mat), n, n)
    
    # First half: e_i + e_{n+1-i}
    for i in 1:half
        P[i, i] = 1
        P[n+1-i, i] = 1
    end
    
    # Second half: e_i - e_{n+1-i}
    for i in 1:half
        P[i, half+i] = 1
        P[n+1-i, half+i] = -1
    end
    
    # Transform: Q_new = P^T * Q * P
    Q_transformed = transpose(P) * mat * P
    
    # Extract blocks (and divide by 2 for normalization)
    block1 = Q_transformed[1:half, 1:half] ./ 2
    block2 = Q_transformed[half+1:end, half+1:end] ./ 2
    
    return (block1, block2, P)
end
```

**Complexity**: O(n³) for matrix multiplication

---

## Eigenvector Computation

### Two Approaches

SymbolicDiagonalization.jl uses two complementary methods:

1. **Adjugate method** (for small matrices, ≤ 3×3)
2. **RREF-based nullspace** (general case)

**File**: `src/rref.jl`, `src/diagonalize.jl`

### Adjugate Method

For A - λI singular (λ is eigenvalue), the adjugate matrix adj(A - λI) has all columns in the nullspace.

```julia
function _adjugate_vectors(M)
    n = size(M, 1)
    n <= 3 || return []  # Only for small matrices
    
    adj = _adjugate(M)
    
    # Return first non-zero column
    for j in 1:n
        col = Symbolics.simplify.(adj[:, j])
        all(_issymzero, col) && continue
        return [col]  # Return immediately
    end
    
    return []
end

function _adjugate(M)
    n = size(M, 1)
    adj = Matrix{eltype(M)}(undef, n, n)
    
    for i in 1:n, j in 1:n
        minor_det = _minor_det(M, i, j)
        adj[j, i] = (-1)^(i + j) * minor_det
    end
    
    return Symbolics.simplify.(adj)
end
```

**Advantages**:
- Avoids RREF pivoting issues
- More compact expressions for small matrices

**Disadvantages**:
- Only works for n ≤ 3 (determinant computation explodes)
- May return zero vector if unlucky

### RREF-Based Nullspace

General method that works for any size:

```julia
function _nullspace(M)
    # Simplify entries first
    simplified = Symbolics.simplify.(M)
    
    # Compute RREF
    R, pivots = _rref(Symbolics.expand.(simplified))
    
    m, n = size(R)
    free = setdiff(1:n, pivots)  # Free variables
    
    vectors = []
    for f in free
        # Set free variable to 1, solve for pivot variables
        v = fill(zero(eltype(R)), n)
        v[f] = one(eltype(R))
        
        for (row, pivot_col) in enumerate(pivots)
            v[pivot_col] = Symbolics.simplify(-R[row, f])
        end
        
        push!(vectors, Symbolics.simplify.(v))
    end
    
    return vectors
end
```

**Complexity**: O(n³) symbolic operations

### RREF Algorithm

Row reduction to reduced row echelon form:

```julia
function _rref(M)
    A = copy(Matrix(M))
    m, n = size(A)
    pivots = Int[]
    row = 1
    
    for col in 1:n
        # Find pivot
        pivot_row = _find_pivot(A, row, col)
        isnothing(pivot_row) && continue
        
        # Swap rows
        if pivot_row != row
            A[row, :], A[pivot_row, :] = A[pivot_row, :], A[row, :]
        end
        
        # Normalize pivot row
        pivot = Symbolics.simplify(A[row, col])
        A[row, :] .= Symbolics.simplify.(A[row, :] ./ pivot)
        
        # Eliminate column in all other rows
        for r in 1:m
            r == row && continue
            factor = A[r, col]
            _issymzero(factor) && continue
            A[r, :] .= Symbolics.simplify.(A[r, :] .- factor .* A[row, :])
        end
        
        push!(pivots, col)
        row += 1
        row > m && break
    end
    
    return A, pivots
end
```

**Key features**:
- Simplification after each operation to keep expressions manageable
- Symbolic pivot selection (uses `_issymzero`)
- Full row elimination (not just below pivot)

**Complexity**: O(mn² × s) where s is average expression size

---

## Expression Management

### The Expression Explosion Problem

Symbolic computation faces a fundamental challenge: **expressions grow exponentially** without careful management.

Example: A 3×3 symbolic matrix can produce eigenvalues with 1000+ terms.
A 4×4 symbolic matrix can produce eigenvalues with 10,000+ terms (~13.5 MB each!).

### Strategy: Aggressive Simplification

We apply simplification at **every intermediate step** of computation:

```julia
function _aggressive_simplify(expr; max_terms = 10000)
    !_is_symbolic_coeff(expr) && return expr
    
    # Expand first to collect all terms
    expanded = Symbolics.expand(expr)
    
    # Check complexity
    _check_expression_size(expanded, max_terms)
    
    # Simplify
    simplified = Symbolics.simplify(expanded)
    
    # Try to factor perfect squares (TODO: not yet implemented)
    factored = _try_factor_perfect_square(simplified)
    
    return factored
end
```

**Applied in**:
- Quadratic discriminant
- Cubic coefficients p, q, Δ
- Quartic coefficients p, q, r, and intermediate β², γ, t₁, t₂
- RREF operations
- Eigenvector construction

### Expression Size Estimation

We estimate expression complexity by counting operations:

```julia
function _estimate_expr_size(expr)
    !_is_symbolic_coeff(expr) && return 1
    
    try
        unwrapped = Symbolics.unwrap(expr)
        return _count_operations(unwrapped)
    catch
        return 1
    end
end

function _count_operations(x)
    # Base case: leaf node
    if x isa Number || !isdefined(x, :f)
        return 1
    end
    
    # Recursive case: count children
    if isdefined(x, :arguments)
        args = getfield(x, :arguments)
        return 1 + sum(_count_operations, args; init=0)
    end
    
    return 1
end
```

### Complexity Threshold Errors

When expressions exceed limits, we throw helpful errors:

```julia
throw(ExpressionComplexityError(
    """Expression has grown too large (≈$size terms, limit: $max_terms).
    
    Suggestions:
    1. Reduce matrix size (try 2×2 or 3×3 instead of 4×4)
    2. Use fewer symbolic variables
    3. Exploit matrix structure
    4. Use numeric eigenvalues
    5. Increase limit with max_terms parameter (caution!)
    """
))
```

### Timeout Mechanism

Long-running computations can be interrupted:

```julia
function _with_timeout(f, timeout_seconds, degree)
    task = @async f()
    timeout_task = @async (sleep(timeout_seconds); true)
    
    # Wait for either task or timeout
    while !istaskdone(task) && !istaskdone(timeout_task)
        sleep(0.1)
    end
    
    if istaskdone(timeout_task)
        # Timeout occurred
        schedule(task, InterruptException(), error=true)
        throw(ComputationTimeoutError("Computation exceeded $timeout_seconds seconds"))
    end
    
    return fetch(task)
end
```

**Default timeout**: 300 seconds (5 minutes)
**Can be disabled**: Set `timeout = nothing`

### Symbolic Square Root for Complex Numbers

Julia's `sqrt(::Complex)` has boolean checks that fail for symbolic values. We implement it directly:

```julia
function _symbolic_sqrt(x)
    !(x isa Complex) && return sqrt(x)
    
    # For Complex{Num}, implement formula manually:
    # sqrt(a + bi) = sqrt((r+a)/2) + i*sign(b)*sqrt((r-a)/2)
    # where r = sqrt(a² + b²)
    
    a = real(x)
    b = imag(x)
    r = sqrt(a^2 + b^2)
    
    real_part = sqrt((r + a) / 2)
    imag_part = sqrt((r - a) / 2)
    
    return Complex(real_part, imag_part)
end
```

---

## Performance Considerations

### Complexity Summary

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Bareiss determinant | O(n³) | Polynomial expression growth |
| Coefficient extraction | O(n) | Using derivatives |
| Quadratic formula | O(n²) | Expression size grows quadratically |
| Cubic formula | O(n³) | Cardano's method |
| Quartic formula | O(n⁴) | Ferrari's method |
| RREF | O(mn² × s) | s = average expression size |
| Nullspace | O(n³ × s) | After RREF |
| Circulant eigenvalues | O(n²) | DFT-based |
| Tridiagonal eigenvalues | O(n) | Closed-form formula |
| Block decomposition | O(k × T(n/k)) | k blocks of size n/k |

### Bottlenecks

1. **Quartic formula**: Produces massive expressions for fully symbolic 4×4 matrices
2. **RREF**: Simplification at each step can be slow for large expressions
3. **Structure detection**: O(n⁴) for some patterns (block circulant, Kronecker)
4. **Symbolic simplification**: Symbolics.jl simplification is not always fast

### Optimization Strategies

**For users**:
1. Use structured matrices (block-diagonal, circulant, tridiagonal)
2. Substitute numeric values for some variables
3. Request eigenvalues only (skip eigenvectors)
4. Use smaller matrices (2×2, 3×3 much faster than 4×4)
5. Set `expand=false` to skip polynomial expansion

**For developers**:
1. Add more special pattern detectors
2. Improve expression simplification (factor perfect squares, etc.)
3. Cache intermediate results
4. Parallelize independent computations (when thread-safe)
5. Add more structural decompositions (Schur complement, etc.)

### Memory Usage

Approximate memory usage for fully symbolic n×n matrices:

| Size | Eigenvalues | With Eigenvectors |
|------|-------------|-------------------|
| 2×2 | ~1 KB | ~5 KB |
| 3×3 | ~100 KB | ~500 KB |
| 4×4 | ~50 MB | ~200 MB |
| 5×5 | N/A* | N/A* |

*5×5 requires special structure (not general case)

### Parallelization Challenges

**Why we don't parallelize**:

Symbolics.jl uses **task-local storage** for hashconsing (expression deduplication), which is **not thread-safe**.

Attempting to use `Threads.@threads` causes crashes:
```julia
# DON'T DO THIS:
Threads.@threads for v in vals
    vecs = _nullspace(mat .- v .* I)  # CRASHES!
end
```

**Potential solution**: Use process-based parallelism (Distributed.jl) instead of threads.

### Test Coverage

The test suite includes 172 tests covering:
- All root solvers (linear through quartic)
- All structure detectors
- All special pattern solvers
- Edge cases (zero matrices, identity, etc.)
- Error handling (timeouts, complexity errors)

**Execution time**: 37.4 seconds

---

## Future Improvements

### Algorithm Enhancements

1. **Schur decomposition**: For upper triangular form
2. **QR algorithm**: Iterative eigenvalue refinement
3. **Lanczos algorithm**: For large sparse matrices
4. **Power method**: For dominant eigenvalue

### Pattern Detection

1. **Hamiltonian matrices**: J-orthogonal structure
2. **Hankel matrices**: Related to Toeplitz
3. **Cauchy matrices**: Explicit determinant formulas
4. **Vandermonde matrices**: Closed-form determinant

### Expression Optimization

1. **Perfect square factoring**: Detect and factor (a-b)² + c²
2. **Common subexpression elimination**: Deduplicate repeated subexpressions
3. **Gröbner basis reduction**: Polynomial ideal membership
4. **Numerical stability analysis**: Detect ill-conditioned expressions

### User Experience

1. **Progress bars**: For long computations
2. **Incremental results**: Return partial results before timeout
3. **Symbolic assumptions**: Propagate assumptions (real, positive, etc.)
4. **Pretty printing**: Better display of large expressions

---

## References

### Algorithms

1. **Bareiss, E.H.** (1968). "Sylvester's identity and multistep integer-preserving Gaussian elimination." *Mathematics of Computation* 22(103): 565-578.

2. **Cardano, G.** (1545). *Ars Magna* (The Great Art). Closed-form solution for cubic equations.

3. **Ferrari, L.** (1545). Solution of quartic equations (published in Cardano's Ars Magna).

4. **Davis, P.J.** (1979). *Circulant Matrices*. Wiley-Interscience. Theory of circulant eigenvalues.

5. **Trench, W.F.** (1999). "Numerical solution of the eigenvalue problem for Hermitian Toeplitz matrices." *SIAM Journal on Matrix Analysis and Applications* 10(2): 135-146.

### Mathematical Background

6. **Abel, N.H.** (1826). "Beweis der Unmöglichkeit, algebraische Gleichungen von höheren Graden als dem vierten allgemein aufzulösen." Proof that degree ≥5 has no general formula.

7. **Galois, É.** (1832). "Mémoire sur les conditions de résolubilité des équations par radicaux." Galois theory foundation.

8. **Horn, R.A. & Johnson, C.R.** (2013). *Matrix Analysis* (2nd ed.). Cambridge University Press. Comprehensive matrix theory reference.

### Implementation

9. **Symbolics.jl Documentation**. https://symbolics.juliasymbolics.org/

10. **SymbolicUtils.jl**. https://symbolicutils.juliasymbolics.org/

---

## Appendix: File Organization

```
src/
├── SymbolicDiagonalization.jl  # Module definition and exports
├── charpoly.jl                 # Characteristic polynomial (Bareiss)
├── roots.jl                    # Root solvers (degrees 1-4)
├── rref.jl                     # RREF and nullspace computation
└── diagonalize.jl              # Main API and pattern detection
    ├── Public API (eigen, eigvals, symbolic_*)
    ├── Structure Detection (_is_*, _detect_*)
    ├── Special Pattern Solvers (_*_eigenvalues)
    ├── Eigenvector Computation (_adjugate, _nullspace)
    └── Utility Functions
```

**Lines of code**:
- `charpoly.jl`: ~60 lines
- `roots.jl`: ~390 lines
- `rref.jl`: ~75 lines
- `diagonalize.jl`: ~1495 lines
- **Total**: ~2020 lines of implementation code

---

*This implementation documentation was last updated: December 11, 2025*
