# ============================================================================
# Anti-Circulant and Kac-Murdock-Szegő Matrices
# ============================================================================
#
# ANTI-CIRCULANT MATRICES
# -----------------------
# An anti-circulant matrix has the form:
#   A[i,j] = c[(i+j-2) mod n + 1]
#
# Example (n=4):
#   [c₁ c₂ c₃ c₄]
#   [c₂ c₃ c₄ c₁]
#   [c₃ c₄ c₁ c₂]
#   [c₄ c₁ c₂ c₃]
#
# This is a Hankel matrix that "wraps around".
# 
# Eigenvalues: Related to circulant eigenvalues with a twist.
# For anti-circulant A with first row [c₁, c₂, ..., cₙ]:
#   λₖ = Σⱼ cⱼ · ω^((j-1)(k-1)) · (-1)^(k-1)  where ω = e^(2πi/n)
#
# Actually, anti-circulant = J · circulant where J is the exchange matrix.
# So eigenvalues are related by: λ_anticirculant = ±λ_circulant
#
# KAC-MURDOCK-SZEGŐ MATRIX
# ------------------------
# K[i,j] = ρ^|i-j| for some ρ with |ρ| < 1
#
# This is a symmetric Toeplitz matrix representing an AR(1) correlation.
#
# Eigenvalues have closed form:
#   λₖ = (1-ρ²) / (1 - 2ρ·cos(θₖ) + ρ²)
# where θₖ are roots of a transcendental equation.
#
# For the simpler case, eigenvalues can be approximated or computed
# via the tridiagonal representation.
#
# ============================================================================

# ============================================================================
# Anti-Circulant Matrices
# ============================================================================

"""
    _anticirculant_matrix(c::AbstractVector)

Construct an anti-circulant matrix from the first row c.

An anti-circulant matrix has A[i,j] = c[(i+j-2) mod n + 1].

# Example
For c = [a, b, c, d]:
```
[a b c d]
[b c d a]
[c d a b]
[d a b c]
```
"""
function _anticirculant_matrix(c::AbstractVector)
    n = length(c)
    A = similar(c, n, n)
    for i in 1:n
        for j in 1:n
            idx = mod(i + j - 2, n) + 1
            A[i, j] = c[idx]
        end
    end
    return A
end

"""
    _is_anticirculant(mat)

Check if a matrix is anti-circulant.

An anti-circulant matrix has constant anti-diagonals that wrap around:
A[i,j] depends only on (i+j) mod n.

Returns the first row if anti-circulant, `nothing` otherwise.
"""
function _is_anticirculant(mat)
    n = size(mat, 1)
    n == size(mat, 2) || return nothing
    n >= 1 || return nothing
    
    # Extract what should be the generating sequence
    first_row = mat[1, :]
    
    # Verify anti-circulant structure
    for i in 1:n
        for j in 1:n
            idx = mod(i + j - 2, n) + 1
            if !_issymzero(mat[i, j] - first_row[idx])
                return nothing
            end
        end
    end
    
    return first_row
end

"""
    _anticirculant_eigenvalues(c::AbstractVector)

Compute eigenvalues of an anti-circulant matrix with first row c.

An anti-circulant matrix A can be written as A = J·C where:
- J is the exchange matrix (1s on anti-diagonal)
- C is a circulant matrix

The eigenvalues are related to circulant eigenvalues.

For an n×n anti-circulant:
- If n is even: eigenvalues come in ± pairs related to circulant eigenvalues
- If n is odd: one eigenvalue is the sum Σcⱼ, others come in pairs

The formula: Let ω = e^(2πi/n). The eigenvalues are:
  λₖ = Σⱼ cⱼ · ω^((j-1)(2k-1)/2)  for k = 1, ..., n

Equivalently, using the relationship with circulant:
The anti-circulant A = P·C where P is a permutation, so eigenvalues
of A are eigenvalues of C (possibly with sign changes).
"""
function _anticirculant_eigenvalues(c::AbstractVector)
    n = length(c)
    
    # The anti-circulant is similar to the circulant with a permutation
    # Eigenvalues: λₖ = Σⱼ cⱼ · ω^((j-1)(k-1)) · phase_factor
    # where the phase depends on the relationship between anti-circ and circ
    
    ω = exp(2π * im / n)
    eigenvalues = Vector{Any}(undef, n)
    
    for k in 0:n-1
        # For anti-circulant, use half-integer frequencies
        # λₖ = Σⱼ cⱼ · ω^((j-1)(k + 1/2))
        λ = zero(ComplexF64)
        for j in 1:n
            λ += c[j] * ω^((j-1) * (k + 0.5))
        end
        eigenvalues[k + 1] = λ
    end
    
    return eigenvalues
end

# ============================================================================
# Kac-Murdock-Szegő Matrix (Geometric Toeplitz)
# ============================================================================

"""
    _kms_matrix(ρ, n::Int)

Construct the Kac-Murdock-Szegő matrix of size n×n with parameter ρ.

K[i,j] = ρ^|i-j|

This represents the correlation matrix of an AR(1) process.

# Example
For n=4, ρ=r:
```
[1   r   r²  r³]
[r   1   r   r²]
[r²  r   1   r ]
[r³  r²  r   1 ]
```
"""
function _kms_matrix(ρ, n::Int)
    n >= 1 || error("Size must be positive")
    K = [ρ^abs(i-j) for i in 1:n, j in 1:n]
    return K
end

"""
    _is_kms_matrix(mat)

Check if a matrix is a Kac-Murdock-Szegő matrix.

A KMS matrix is symmetric Toeplitz with K[i,j] = ρ^|i-j| for some ρ.

Returns (ρ, n) if it's a KMS matrix, `nothing` otherwise.
"""
function _is_kms_matrix(mat)
    n = size(mat, 1)
    n == size(mat, 2) || return nothing
    n >= 1 || return nothing
    
    # Diagonal should be 1
    for i in 1:n
        if !_issymzero(mat[i, i] - 1)
            return nothing
        end
    end
    
    if n == 1
        return (0, 1)  # Any ρ works for n=1, return 0 as default
    end
    
    # Extract ρ from the (1,2) element
    ρ = mat[1, 2]
    
    # Verify Toeplitz structure with geometric decay
    for i in 1:n
        for j in 1:n
            expected = ρ^abs(i - j)
            if !_issymzero(mat[i, j] - expected)
                return nothing
            end
        end
    end
    
    return (ρ, n)
end

"""
    _kms_eigenvalues(ρ, n::Int)

Compute eigenvalues of the Kac-Murdock-Szegő matrix with parameter ρ.

The KMS matrix K with K[i,j] = ρ^|i-j| has a known closed-form for eigenvalues.

The eigenvalues are:
    λₖ = (1 - ρ²) / (1 - 2ρ·cos(θₖ) + ρ²)

where θₖ = kπ/(n+1) for k = 1, 2, ..., n.

This formula comes from the inverse of the tridiagonal structure:
K⁻¹ is tridiagonal with known entries, and we use the tridiagonal eigenvalue formula.

Actually, the exact formula involves the relationship:
For |ρ| < 1, the KMS matrix is positive definite with eigenvalues:
    λₖ = (1 - ρ²) / (1 + ρ² - 2ρ·cos(kπ/(n+1)))

for k = 1, 2, ..., n.
"""
function _kms_eigenvalues(ρ, n::Int)
    n >= 1 || error("Size must be positive")
    
    eigenvalues = Vector{Any}(undef, n)
    
    # Handle special cases
    if _issymzero(ρ)
        # ρ = 0: Identity matrix
        for k in 1:n
            eigenvalues[k] = 1
        end
        return eigenvalues
    end
    
    if _issymzero(ρ - 1)
        # ρ = 1: All-ones matrix, eigenvalues are n (once) and 0 (n-1 times)
        eigenvalues[1] = n
        for k in 2:n
            eigenvalues[k] = 0
        end
        return eigenvalues
    end
    
    if _issymzero(ρ + 1)
        # ρ = -1: Alternating pattern
        # Eigenvalues alternate based on parity
        for k in 1:n
            θ = k * π / (n + 1)
            eigenvalues[k] = (1 - 1) / (1 + 1 - 2*(-1)*cos(θ))  # = 0 / (2 + 2cos(θ))
        end
        # Actually for ρ = -1, need special handling
        # K[i,j] = (-1)^|i-j|, which alternates
        for k in 1:n
            θ = k * π / (n + 1)
            denom = 1 + 1 + 2*cos(θ)  # 2(1 + cos(θ))
            if abs(denom) < 1e-14
                eigenvalues[k] = 0
            else
                eigenvalues[k] = 0 / denom  # numerator is 1 - 1 = 0
            end
        end
        # For ρ = -1, the matrix has a special structure
        # Let's compute directly for this edge case
        return [0 for _ in 1:n]  # Singular matrix when ρ = ±1 and n > 1
    end
    
    # General case: |ρ| ≠ 1
    ρ_sq = ρ^2
    one_minus_ρ_sq = 1 - ρ_sq
    
    for k in 1:n
        θ = k * π / (n + 1)
        cos_θ = cos(θ)
        denom = 1 + ρ_sq - 2 * ρ * cos_θ
        eigenvalues[k] = one_minus_ρ_sq / denom
    end
    
    return eigenvalues
end

"""
    _kms_eigenvalues_symbolic(ρ, n::Int)

Compute symbolic eigenvalues of the KMS matrix.

Returns expressions involving ρ and trigonometric functions.
"""
function _kms_eigenvalues_symbolic(ρ, n::Int)
    n >= 1 || error("Size must be positive")
    
    eigenvalues = Vector{Any}(undef, n)
    
    for k in 1:n
        θ = k * π / (n + 1)
        # λₖ = (1 - ρ²) / (1 + ρ² - 2ρ·cos(θₖ))
        cos_θ = cos(θ)
        numerator = 1 - ρ^2
        denominator = 1 + ρ^2 - 2 * ρ * cos_θ
        eigenvalues[k] = numerator / denominator
    end
    
    return eigenvalues
end

# ============================================================================
# Public API
# ============================================================================

"""
    anticirculant_matrix(c::AbstractVector)

Construct an anti-circulant matrix from the first row c.

An anti-circulant has constant anti-diagonals that wrap around:
A[i,j] = c[(i+j-2) mod n + 1]

# Examples
```julia
julia> anticirculant_matrix([1, 2, 3, 4])
4×4 Matrix{Int64}:
 1  2  3  4
 2  3  4  1
 3  4  1  2
 4  1  2  3
```
"""
anticirculant_matrix(c::AbstractVector) = _anticirculant_matrix(c)

"""
    kms_matrix(ρ, n::Int)

Construct the Kac-Murdock-Szegő matrix of size n×n with parameter ρ.

K[i,j] = ρ^|i-j|

This symmetric Toeplitz matrix represents the correlation matrix of an AR(1) process.

Eigenvalues: λₖ = (1-ρ²) / (1 + ρ² - 2ρ·cos(kπ/(n+1))) for k = 1, ..., n

# Examples
```julia
julia> kms_matrix(0.5, 3)
3×3 Matrix{Float64}:
 1.0   0.5   0.25
 0.5   1.0   0.5
 0.25  0.5   1.0
```
"""
kms_matrix(ρ, n::Int) = _kms_matrix(ρ, n)
