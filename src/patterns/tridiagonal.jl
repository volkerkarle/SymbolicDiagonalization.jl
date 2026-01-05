# ============================================================================
# Structural Pattern: Tridiagonal Matrices
# ============================================================================
#
# Symmetric Toeplitz Tridiagonal:
#   These matrices have known eigenbasis (discrete sine transform). The eigenvalue
#   equation reduces to a three-term recurrence defining Chebyshev polynomials,
#   giving the closed form: λₖ = a + 2b·cos(kπ/(n+1))
#
# Special 5×5 Patterns:
#   Certain near-Toeplitz tridiagonal patterns with exactly one different entry
#   also admit closed-form solutions discovered empirically.
#
# Anti-diagonal Patterns:
#   Matrices with anti-diagonal structure have eigenvalues in ±pairs due to
#   the involutory nature of the flip transformation.
#
# While not strictly group-theoretic, these patterns have solvable structure.
# ============================================================================

"""
    _is_toeplitz_tridiagonal(mat)

Check if a matrix is a symmetric Toeplitz tridiagonal matrix with constant diagonals.
Returns (a, b, c) where:
- a is the main diagonal constant
- b is the subdiagonal constant  
- c is the superdiagonal constant

Returns nothing if the matrix is not of this form or if it's not symmetric.

A symmetric Toeplitz tridiagonal matrix has the form:
```
[a  b  0  0  ...]
[b  a  b  0  ...]
[0  b  a  b  ...]
[...          ...]
[0  0  0  b  a]
```

Note: We only handle the symmetric case (b = c) because the eigenvalue formula
is guaranteed to work for symmetric matrices. For asymmetric tridiagonal matrices,
the formula may not apply if the matrix is defective (non-diagonalizable).
"""
function _is_toeplitz_tridiagonal(mat)
    n = size(mat, 1)
    n == size(mat, 2) || return nothing
    n <= 1 && return nothing  # Need at least 2×2
    
    # Check that matrix is symmetric first
    if !_is_symmetric(mat)
        return nothing
    end
    
    # Check that matrix is tridiagonal (zeros beyond ±1 diagonals)
    for i in 1:n, j in 1:n
        if abs(i - j) > 1 && !_issymzero(mat[i, j])
            return nothing
        end
    end
    
    # Extract diagonal constants
    a = mat[1, 1]
    
    # Check main diagonal is constant
    for i in 2:n
        if !_issymzero(mat[i, i] - a)
            return nothing
        end
    end
    
    # Extract and check subdiagonal (if it exists)
    b = n > 1 ? mat[2, 1] : zero(eltype(mat))
    for i in 3:n
        if !_issymzero(mat[i, i-1] - b)
            return nothing
        end
    end
    
    # For symmetric case, superdiagonal equals subdiagonal
    c = b
    
    return (a, b, c)
end

"""
    _toeplitz_tridiagonal_eigenvalues(n, a, b, c)

Compute eigenvalues of an n×n symmetric Toeplitz tridiagonal matrix with constants (a, b, c).

For a symmetric tridiagonal matrix (b = c), the eigenvalues are:
    λₖ = a + 2b·cos(kπ/(n+1))  for k = 1, 2, ..., n

This formula comes from the theory of orthogonal polynomials and is exact.
"""
function _toeplitz_tridiagonal_eigenvalues(n, a, b, c)
    eigenvalues = Vector{Any}(undef, n)
    
    # For the symmetric case: λₖ = a + 2b·cos(kπ/(n+1))
    # (We verified b = c in the detection function)
    
    for k in 1:n
        θ = k * π / (n + 1)
        λₖ = a + 2 * b * cos(θ)
        eigenvalues[k] = Symbolics.simplify(λₖ)
    end
    
    return eigenvalues
end

"""
Detect specific 5×5 tridiagonal patterns with closed-form eigenvalues.

Pattern 1:                    Pattern 2:
    [a  b  0  0  0]              [a  b  0  0  0]
    [b  a  d  0  0]              [b  a  b  0  0]
    [0  d  a  b  0]     OR       [0  b  a  d  0]
    [0  0  b  a  b]              [0  0  d  a  b]
    [0  0  0  b  a]              [0  0  0  b  a]
    
Both patterns have identical closed-form eigenvalues:
    λ₁ = a - √(2b² + d²)
    λ₂ = a - b
    λ₃ = a
    λ₄ = a + b
    λ₅ = a + √(2b² + d²)
    
Returns (a, b, d) if pattern matches, nothing otherwise.
"""
function _detect_special_5x5_tridiagonal(mat)
    size(mat) == (5, 5) || return nothing
    
    # Check that matrix is tridiagonal (zeros beyond ±1 diagonals)
    for i in 1:5, j in 1:5
        if abs(i - j) > 1 && !_issymzero(mat[i, j])
            return nothing
        end
    end
    
    # Extract diagonal and off-diagonal elements
    a = mat[1, 1]
    
    # Check if all diagonal elements are equal
    for i in 2:5
        if !_issymzero(mat[i, i] - a)
            return nothing
        end
    end
    
    # Check if symmetric
    if !_is_symmetric(mat)
        return nothing
    end
    
    # Extract off-diagonal pattern: should be [b, d, b, b] or [b, b, d, b]
    off1 = mat[1, 2]
    off2 = mat[2, 3]
    off3 = mat[3, 4]
    off4 = mat[4, 5]
    
    # Pattern 1: [b, d, b, b] - positions 1, 3, 4 are the same
    if _issymzero(off1 - off3) && _issymzero(off1 - off4)
        b = off1
        d = off2
        return (a, b, d)
    end
    
    # Pattern 2: [b, b, d, b] - positions 1, 2, 4 are the same
    if _issymzero(off1 - off2) && _issymzero(off1 - off4)
        b = off1
        d = off3
        return (a, b, d)
    end
    
    return nothing
end

"""
    _is_antidiagonal(mat)

Check if a matrix is anti-diagonal (non-zero only on anti-diagonal).
"""
function _is_antidiagonal(mat)
    m, n = size(mat)
    m == n || return false
    # Check that only anti-diagonal entries (i + j = n + 1) are non-zero
    for i in 1:m, j in 1:n
        if i + j == n + 1
            continue  # Anti-diagonal entry, can be anything
        else
            !_issymzero(mat[i, j]) && return false
        end
    end
    return true
end

"""
    _antidiagonal_eigenvalues(mat)

Compute eigenvalues of an anti-diagonal matrix (non-zero only on anti-diagonal).

An anti-diagonal matrix has the form:
```
[0  0  ... 0  a₁]
[0  0  ... a₂ 0 ]
[     ...       ]
[0  aₙ₋₁ ... 0 0]
[aₙ 0  ... 0  0 ]
```

For a symmetric anti-diagonal matrix, the eigenvalues are ±aₖ (with sign pattern
depending on whether n is even or odd). For asymmetric cases, we use the fact that
the matrix is similar to a diagonal matrix via a permutation.
"""
function _antidiagonal_eigenvalues(mat)
    n = size(mat, 1)
    # Extract anti-diagonal elements
    antidiag = [mat[i, n + 1 - i] for i in 1:n]
    
    # For symmetric anti-diagonal: eigenvalues come in ±pairs
    # For general case: we need to be more careful
    # The eigenvalues depend on the parity structure
    
    # For now, handle the symmetric case
    if _is_symmetric(mat)
        eigenvalues = Vector{Any}(undef, n)
        if n % 2 == 1
            # Odd dimension: one zero eigenvalue, rest come in pairs
            mid = (n + 1) ÷ 2
            eigenvalues[1] = antidiag[mid]
            idx = 2
            for i in 1:mid-1
                eigenvalues[idx] = antidiag[i]
                eigenvalues[idx+1] = -antidiag[i]
                idx += 2
            end
        else
            # Even dimension: all come in ±pairs
            idx = 1
            for i in 1:n÷2
                eigenvalues[idx] = antidiag[i]
                eigenvalues[idx+1] = -antidiag[i]
                idx += 2
            end
        end
        return eigenvalues
    else
        # For non-symmetric, the analysis is more complex
        # Return nothing to fall back to general method
        return nothing
    end
end
