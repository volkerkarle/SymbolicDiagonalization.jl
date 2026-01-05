# ============================================================================
# Sp(n) - Symplectic Groups
# ============================================================================
#
# Symplectic groups Sp(2n) consist of 2n x 2n matrices A satisfying A^T J A = J
# where J = [0 I; -I 0] is the standard symplectic form.
#
# Key property: eigenvalues come in reciprocal pairs (lambda, 1/lambda).
# This makes the characteristic polynomial palindromic.
#
# Supported:
# - Sp(2) = SL(2,R): 2x2 symplectic (equivalent to 2x2 real matrices with det=1)
# - Sp(4): 4x4 symplectic matrices
#
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Sp(2) - 2x2 Symplectic (equivalent to SL(2,R))
# ============================================================================

"""
    _is_Sp2(A)

Check if A is a 2x2 symplectic matrix (Sp(2) = SL(2,R)).

A 2x2 matrix is symplectic iff det(A) = 1 (since the symplectic condition
reduces to the determinant condition in 2D).
"""
function _is_Sp2(A)
    size(A) == (2, 2) || return false
    n = _is_symplectic(A)
    return !isnothing(n) && n == 1
end

"""
    _Sp2_eigenvalues(A)

Compute eigenvalues of an Sp(2) matrix symbolically.

For Sp(2) with det = 1, the characteristic polynomial is:
    lambda^2 - tr(A)*lambda + 1 = 0

The eigenvalues are:
    lambda = (tr(A) +/- sqrt(tr(A)^2 - 4)) / 2

Note: eigenvalues satisfy lambda_1 * lambda_2 = 1 (reciprocal pair).
"""
function _Sp2_eigenvalues(A)
    _is_Sp2(A) || return nothing
    
    t = tr(A)
    discriminant = t^2 - 4
    sqrt_disc = sqrt(discriminant)
    
    lambda1 = (t + sqrt_disc) / 2
    lambda2 = (t - sqrt_disc) / 2
    
    return [Symbolics.simplify(lambda1), Symbolics.simplify(lambda2)]
end

# ============================================================================
# Sp(2) Eigenvectors
# ============================================================================

"""
    _Sp2_is_diagonal(A)

Check if A is a diagonal Sp(2) matrix of the form [a 0; 0 1/a].
Returns (a,) if detected, nothing otherwise.
"""
function _Sp2_is_diagonal(A)
    size(A) == (2, 2) || return nothing
    
    # Check off-diagonal elements are zero
    if !_issymzero(A[1, 2]) || !_issymzero(A[2, 1])
        return nothing
    end
    
    # For diagonal Sp(2): a * (1/a) = 1
    # Just verify the determinant condition (which _is_Sp2 already checks)
    return (A[1, 1],)
end

"""
    _Sp2_eigenpairs(A)

Compute eigenvalue-eigenvector pairs for a non-diagonal Sp(2) matrix.

For general Sp(2) matrices, eigenvectors are computed analytically:
- For A = [a b; c d] with eigenvalue λ, eigenvector is [b, λ - a] if b ≠ 0

NOTE: Diagonal Sp(2) [a 0; 0 1/a] returns nothing because:
- Eigenvalues are trivially a and 1/a (on the diagonal)
- Eigenvectors are trivially [1, 0] and [0, 1] (standard basis)

Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing.
"""
function _Sp2_eigenpairs(A)
    _is_Sp2(A) || return nothing
    
    # Check for diagonal structure - if diagonal, return nothing (trivial)
    diag_info = _Sp2_is_diagonal(A)
    if !isnothing(diag_info)
        # Diagonal case is trivial - skip
        return nothing
    end
    
    # For general 2x2 matrix, we can compute eigenvectors analytically
    # For A = [a b; c d] with eigenvalue λ, eigenvector satisfies (A - λI)v = 0
    # If b ≠ 0: v = [b, λ - a] (unnormalized)
    # If c ≠ 0: v = [λ - d, c]
    
    eigenvalues = _Sp2_eigenvalues(A)
    if isnothing(eigenvalues)
        return nothing
    end
    
    λ1, λ2 = eigenvalues
    a, b, c, d = A[1,1], A[1,2], A[2,1], A[2,2]
    
    # Use the formula v = [b, λ - a] if b is non-zero
    # Otherwise use v = [λ - d, c]
    if !_issymzero(b)
        v1 = [b, λ1 - a]
        v2 = [b, λ2 - a]
    elseif !_issymzero(c)
        v1 = [λ1 - d, c]
        v2 = [λ2 - d, c]
    else
        # Both b and c are zero means diagonal - should have been caught above
        return nothing
    end
    
    return [
        (λ1, [v1]),
        (λ2, [v2])
    ]
end

# ============================================================================
# Sp(4) - 4x4 Symplectic
# ============================================================================

"""
    _is_Sp4(A)

Check if A is a 4x4 symplectic matrix (Sp(4)).
"""
function _is_Sp4(A)
    size(A) == (4, 4) || return false
    n = _is_symplectic(A)
    return !isnothing(n) && n == 2
end

"""
    _Sp4_eigenvalues(A)

Compute eigenvalues of an Sp(4) matrix symbolically.

For Sp(4), eigenvalues come in reciprocal pairs: {lambda, 1/lambda, mu, 1/mu}.
The characteristic polynomial is palindromic:
    lambda^4 - a*lambda^3 + b*lambda^2 - a*lambda + 1 = 0

where a = tr(A) and b = (tr(A)^2 - tr(A^2))/2.

This reduces to a quadratic in nu = lambda + 1/lambda:
    nu^2 - a*nu + (b - 2) = 0

Then for each nu, solve lambda^2 - nu*lambda + 1 = 0.
"""
function _Sp4_eigenvalues(A)
    _is_Sp4(A) || return nothing
    
    # Compute trace invariants
    t1 = tr(A)
    t2 = tr(A * A)
    
    # Coefficients of palindromic polynomial
    # c3 = a (coefficient of lambda^3 and lambda^1)
    # c2 = b (coefficient of lambda^2)
    c3 = t1
    c2 = (t1^2 - t2) / 2
    
    # Solve quadratic for nu = lambda + 1/lambda:
    # nu^2 - c3*nu + (c2 - 2) = 0
    a_coeff = c3
    b_coeff = c2 - 2
    
    disc_nu = a_coeff^2 - 4 * b_coeff
    sqrt_disc_nu = sqrt(disc_nu)
    
    nu1 = (a_coeff + sqrt_disc_nu) / 2
    nu2 = (a_coeff - sqrt_disc_nu) / 2
    
    # For each nu, solve lambda^2 - nu*lambda + 1 = 0
    # lambda = (nu +/- sqrt(nu^2 - 4)) / 2
    
    disc1 = nu1^2 - 4
    sqrt_disc1 = sqrt(disc1)
    lambda1 = (nu1 + sqrt_disc1) / 2
    lambda2 = (nu1 - sqrt_disc1) / 2
    
    disc2 = nu2^2 - 4
    sqrt_disc2 = sqrt(disc2)
    lambda3 = (nu2 + sqrt_disc2) / 2
    lambda4 = (nu2 - sqrt_disc2) / 2
    
    return [Symbolics.simplify(lambda1), Symbolics.simplify(lambda2),
            Symbolics.simplify(lambda3), Symbolics.simplify(lambda4)]
end

# ============================================================================
# Sp(4) Eigenvectors
# ============================================================================

"""
    _Sp4_is_diagonal(A)

Check if A is a fully diagonal Sp(4) matrix.

With the standard symplectic form J = [0 I; -I 0], a diagonal Sp(4) matrix
has the form diag(a, b, 1/a, 1/b) where a*1/a = 1 and b*1/b = 1.

Returns (a, b) if detected, nothing otherwise.
"""
function _Sp4_is_diagonal(A)
    size(A) == (4, 4) || return nothing
    
    # Check all off-diagonal elements are zero
    for i in 1:4, j in 1:4
        if i != j && !_issymzero(A[i, j])
            return nothing
        end
    end
    
    # For symplectic with J = [0 I; -I 0]:
    # A^T J A = J requires A[1,1]*A[3,3] = 1 and A[2,2]*A[4,4] = 1
    return (A[1, 1], A[2, 2])
end

"""
    _Sp4_eigenpairs(A)

Compute eigenvalue-eigenvector pairs for a non-diagonal Sp(4) matrix.

NOTE: Diagonal Sp(4) diag(a, b, 1/a, 1/b) returns nothing because:
- Eigenvalues are trivially a, b, 1/a, 1/b (on the diagonal)
- Eigenvectors are trivially standard basis vectors

For general Sp(4), returns nothing (requires nullspace computation).
Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing.
"""
function _Sp4_eigenpairs(A)
    _is_Sp4(A) || return nothing
    
    # Check for fully diagonal structure - if diagonal, return nothing (trivial)
    diag_info = _Sp4_is_diagonal(A)
    if !isnothing(diag_info)
        # Diagonal case is trivial - skip
        return nothing
    end
    
    # For general Sp(4), we need nullspace computation
    return nothing
end

# ============================================================================
# General Sp(2n) detection (for dispatch)
# ============================================================================

"""
    _is_Sp(A)

Check if A is a symplectic matrix of any dimension.
Returns the half-dimension n if A is in Sp(2n), nothing otherwise.
"""
function _is_Sp(A)
    return _is_symplectic(A)
end

"""
    _Sp_eigenvalues(A)

Compute eigenvalues of a symplectic matrix if supported.

Currently supports Sp(2) and Sp(4). Returns nothing for larger dimensions.
"""
function _Sp_eigenvalues(A)
    n = size(A, 1)
    
    if n == 2
        return _Sp2_eigenvalues(A)
    elseif n == 4
        return _Sp4_eigenvalues(A)
    end
    
    # For Sp(2n) with n > 2, we don't have closed-form solutions
    return nothing
end
