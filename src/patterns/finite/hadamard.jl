# ============================================================================
# Sylvester-Hadamard and DFT Matrices
# ============================================================================
#
# These matrices have elegant closed-form eigenvalue formulas.
#
# SYLVESTER-HADAMARD MATRICES
# ---------------------------
# Recursive construction: H₁ = [1], Hₙ = [Hₙ₋₁  Hₙ₋₁; Hₙ₋₁ -Hₙ₋₁]
# Size is always 2ⁿ × 2ⁿ
# 
# Properties:
#   - Symmetric: Hᵀ = H
#   - Orthogonal (up to scaling): H·Hᵀ = 2ⁿ·I
#   - Entries are ±1
#   - Eigenvalues: ±√(2ⁿ) = ±2^(n/2)
#   - Multiplicity: each eigenvalue has multiplicity 2^(n-1)
#
# DFT MATRIX (Fourier Matrix)
# ---------------------------
# F[j,k] = ω^(jk) where ω = e^(2πi/n), indices 0-based
# Or normalized: F[j,k] = ω^(jk)/√n
#
# Properties:
#   - Unitary (when normalized): F·F† = I
#   - Eigenvalues are 4th roots of unity: {1, -1, i, -i}
#   - Multiplicities depend on n mod 4
#
# ============================================================================

# ============================================================================
# Sylvester-Hadamard Matrices
# ============================================================================

"""
    _hadamard_matrix(n::Int)

Construct the Sylvester-Hadamard matrix of order 2ⁿ.

Uses the recursive construction:
- H₀ = [1]
- Hₙ = [Hₙ₋₁  Hₙ₋₁; Hₙ₋₁ -Hₙ₋₁]

Returns a 2ⁿ × 2ⁿ matrix with entries ±1.
"""
function _hadamard_matrix(n::Int)
    n >= 0 || error("Hadamard order must be non-negative")
    
    if n == 0
        return reshape([1], 1, 1)  # Return 1×1 matrix, not a vector
    end
    
    H_prev = _hadamard_matrix(n - 1)
    return [H_prev H_prev; H_prev -H_prev]
end

"""
    _is_hadamard_sylvester(mat)

Check if a matrix is a Sylvester-Hadamard matrix.

Criteria:
1. Size is 2ⁿ × 2ⁿ for some n ≥ 0
2. All entries are ±1
3. H·Hᵀ = 2ⁿ·I (orthogonality up to scaling)

Returns the order n if it's a Sylvester-Hadamard matrix, `nothing` otherwise.
"""
function _is_hadamard_sylvester(mat)
    m, n_cols = size(mat)
    m == n_cols || return nothing
    
    # Check size is power of 2
    m > 0 || return nothing
    n = round(Int, log2(m))
    2^n == m || return nothing
    
    # Check all entries are ±1
    for i in 1:m
        for j in 1:m
            val = mat[i, j]
            if !(_issymzero(val - 1) || _issymzero(val + 1))
                return nothing
            end
        end
    end
    
    # Check orthogonality: H·Hᵀ = 2ⁿ·I
    HHt = mat * mat'
    for i in 1:m
        for j in 1:m
            expected = (i == j) ? m : 0
            if !_issymzero(HHt[i, j] - expected)
                return nothing
            end
        end
    end
    
    return n
end

"""
    _hadamard_eigenvalues(n::Int)

Compute eigenvalues of the Sylvester-Hadamard matrix of order 2ⁿ.

The eigenvalues are ±2^(n/2) = ±√(2ⁿ), each with multiplicity 2^(n-1).

For n=0: eigenvalue is 1 (the 1×1 matrix [1])
For n≥1: eigenvalues are +2^(n/2) and -2^(n/2), each with multiplicity 2^(n-1)
"""
function _hadamard_eigenvalues(n::Int)
    if n == 0
        return [1]
    end
    
    size = 2^n
    multiplicity = 2^(n - 1)
    
    # Eigenvalue magnitude: √(2ⁿ) = 2^(n/2)
    λ = 2^(n // 2)  # Use rational to keep exact when possible
    
    eigenvalues = Vector{Any}(undef, size)
    
    # First half: +√(2ⁿ)
    for i in 1:multiplicity
        eigenvalues[i] = sqrt(Rational(2)^n)
    end
    
    # Second half: -√(2ⁿ)
    for i in 1:multiplicity
        eigenvalues[multiplicity + i] = -sqrt(Rational(2)^n)
    end
    
    return eigenvalues
end

# ============================================================================
# DFT Matrix (Discrete Fourier Transform)
# ============================================================================

"""
    _dft_matrix(n::Int; normalized::Bool=false)

Construct the DFT matrix of size n×n.

F[j,k] = ω^((j-1)(k-1)) where ω = e^(2πi/n), using 1-based indexing.

If normalized=true, returns F/√n which is unitary.
"""
function _dft_matrix(n::Int; normalized::Bool=false)
    n >= 1 || error("DFT size must be positive")
    
    ω = exp(2π * im / n)
    F = [ω^((j-1)*(k-1)) for j in 1:n, k in 1:n]
    
    if normalized
        F = F / sqrt(n)
    end
    
    return F
end

"""
    _is_dft_matrix(mat; normalized::Bool=false)

Check if a matrix is a DFT matrix.

Returns the size n if it's a DFT matrix, `nothing` otherwise.
"""
function _is_dft_matrix(mat; normalized::Bool=false)
    m, n = size(mat)
    m == n || return nothing
    n >= 1 || return nothing
    
    # Skip symbolic matrices - DFT detection is for numeric matrices only
    # (comparing symbolic expressions with numeric tolerance doesn't make sense)
    T = eltype(mat)
    if T <: Num || T <: Complex{<:Num}
        return nothing
    end
    
    # Construct expected DFT matrix
    expected = _dft_matrix(n, normalized=normalized)
    
    tol = 1e-10
    for i in 1:n
        for j in 1:n
            diff = mat[i, j] - expected[i, j]
            if abs(diff) > tol
                return nothing
            end
        end
    end
    
    return n
end

"""
    _dft_eigenvalues(n::Int)

Compute eigenvalues of the DFT matrix Fₙ.

The DFT matrix has eigenvalues that are scaled 4th roots of unity: ±√n, ±i√n.

For ω = e^{+2πi/n} (as used in _dft_matrix), the multiplicities are:
- m(√n)  = ceiling((n+1)/4)
- m(-√n) = floor(n/4)  
- m(i√n) = floor(n/4)
- m(-i√n) = floor((n-1)/4)

These satisfy m(√n) + m(-√n) + m(i√n) + m(-i√n) = n.
"""
function _dft_eigenvalues(n::Int)
    n >= 1 || error("DFT size must be positive")
    
    if n == 1
        return [1.0 + 0.0im]
    end
    
    # The correct multiplicities for DFT eigenvalues with ω = e^{+2πi/n}
    # (positive omega convention, as used in _dft_matrix)
    # Note: For ω = e^{-2πi/n}, m_i and m_negi would be swapped.
    r = n % 4
    
    if r == 0
        m1 = n ÷ 4 + 1
        m_neg1 = n ÷ 4
        m_i = n ÷ 4          # swapped from n ÷ 4 - 1
        m_negi = n ÷ 4 - 1   # swapped from n ÷ 4
    elseif r == 1
        m1 = (n + 3) ÷ 4
        m_neg1 = (n - 1) ÷ 4
        m_i = (n - 1) ÷ 4
        m_negi = (n - 1) ÷ 4
    elseif r == 2
        m1 = (n + 2) ÷ 4
        m_neg1 = (n + 2) ÷ 4
        m_i = (n - 2) ÷ 4
        m_negi = (n - 2) ÷ 4
    else  # r == 3
        m1 = (n + 1) ÷ 4
        m_neg1 = (n + 1) ÷ 4
        m_i = (n + 3) ÷ 4    # swapped from (n - 1) ÷ 4
        m_negi = (n - 1) ÷ 4 # swapped from (n + 3) ÷ 4
    end
    
    # Verify total
    total = m1 + m_neg1 + m_i + m_negi
    @assert total == n "Multiplicity error: got $total, expected $n"
    
    sqrtn = sqrt(n)
    
    eigenvalues = Vector{ComplexF64}(undef, n)
    idx = 1
    
    for _ in 1:m1
        eigenvalues[idx] = sqrtn + 0.0im
        idx += 1
    end
    for _ in 1:m_neg1
        eigenvalues[idx] = -sqrtn + 0.0im
        idx += 1
    end
    for _ in 1:m_i
        eigenvalues[idx] = 0.0 + sqrtn*im
        idx += 1
    end
    for _ in 1:m_negi
        eigenvalues[idx] = 0.0 - sqrtn*im
        idx += 1
    end
    
    return eigenvalues
end

"""
    _dft_eigenvalues_normalized(n::Int)

Compute eigenvalues of the normalized DFT matrix Fₙ/√n.

The normalized DFT is unitary with eigenvalues being exactly 4th roots of unity:
{1, -1, i, -i}.
"""
function _dft_eigenvalues_normalized(n::Int)
    n >= 1 || error("DFT size must be positive")
    
    if n == 1
        return [Complex{Float64}(1)]
    end
    
    # Multiplicities for normalized DFT with ω = e^{+2πi/n} (positive omega convention)
    # Eigenvalues are 4th roots of unity: {1, -1, i, -i}
    r = n % 4
    
    if r == 0
        m1 = n ÷ 4 + 1
        m_neg1 = n ÷ 4
        m_i = n ÷ 4          # swapped for positive omega
        m_negi = n ÷ 4 - 1   # swapped for positive omega
    elseif r == 1
        m1 = (n + 3) ÷ 4
        m_neg1 = (n - 1) ÷ 4
        m_i = (n - 1) ÷ 4
        m_negi = (n - 1) ÷ 4
    elseif r == 2
        m1 = (n + 2) ÷ 4
        m_neg1 = (n + 2) ÷ 4
        m_i = (n - 2) ÷ 4
        m_negi = (n - 2) ÷ 4
    else  # r == 3
        m1 = (n + 1) ÷ 4
        m_neg1 = (n + 1) ÷ 4
        m_i = (n + 3) ÷ 4    # swapped for positive omega
        m_negi = (n - 1) ÷ 4 # swapped for positive omega
    end
    
    eigenvalues = Vector{Complex{Float64}}(undef, n)
    idx = 1
    
    for _ in 1:m1
        eigenvalues[idx] = 1.0 + 0.0im
        idx += 1
    end
    for _ in 1:m_neg1
        eigenvalues[idx] = -1.0 + 0.0im
        idx += 1
    end
    for _ in 1:m_i
        eigenvalues[idx] = 0.0 + 1.0im
        idx += 1
    end
    for _ in 1:m_negi
        eigenvalues[idx] = 0.0 - 1.0im
        idx += 1
    end
    
    return eigenvalues
end

# ============================================================================
# Public API
# ============================================================================

"""
    hadamard_matrix(n::Int)

Construct the Sylvester-Hadamard matrix of order 2ⁿ.

The Hadamard matrix is constructed recursively:
- H₀ = [1]
- Hₙ = [Hₙ₋₁  Hₙ₋₁; Hₙ₋₁ -Hₙ₋₁]

Properties:
- Symmetric with entries ±1
- Satisfies H·Hᵀ = 2ⁿ·I
- Eigenvalues: ±√(2ⁿ), each with multiplicity 2^(n-1)

# Examples
```julia
julia> hadamard_matrix(2)
4×4 Matrix{Int64}:
 1   1   1   1
 1  -1   1  -1
 1   1  -1  -1
 1  -1  -1   1
```
"""
hadamard_matrix(n::Int) = _hadamard_matrix(n)

"""
    dft_matrix(n::Int; normalized::Bool=false)

Construct the Discrete Fourier Transform matrix of size n×n.

F[j,k] = ω^((j-1)(k-1)) where ω = e^(2πi/n).

If `normalized=true`, returns F/√n which is unitary.

Eigenvalues (unnormalized): ±√n, ±i√n
Eigenvalues (normalized): 1, -1, i, -i

# Examples
```julia
julia> dft_matrix(4)
4×4 Matrix{ComplexF64}:
 1.0+0.0im   1.0+0.0im   1.0+0.0im   1.0+0.0im
 1.0+0.0im   0.0+1.0im  -1.0+0.0im  -0.0-1.0im
 1.0+0.0im  -1.0+0.0im   1.0-0.0im  -1.0+0.0im
 1.0+0.0im  -0.0-1.0im  -1.0+0.0im   0.0+1.0im
```
"""
dft_matrix(n::Int; normalized::Bool=false) = _dft_matrix(n, normalized=normalized)
