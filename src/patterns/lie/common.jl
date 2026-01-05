# ============================================================================
# Compact Lie Groups: Common Detection Utilities
# ============================================================================
#
# This module provides shared infrastructure for detecting matrices belonging
# to compact Lie groups:
#
#   SO(n) - Special orthogonal: AᵀA = I, det(A) = 1
#   SU(n) - Special unitary: A†A = I, det(A) = 1  
#   Sp(2n) - Symplectic: AᵀJA = J where J = [0 I; -I 0]
#
# Compact Lie groups have eigenvalues on the unit circle. Their structure
# allows eigenvalue computation via trace invariants rather than polynomial
# root-finding, bypassing the Abel-Ruffini theorem.
#
# Individual group implementations are in SO2.jl, SO3.jl, SO4.jl, SU2.jl,
# SU3.jl, and Sp.jl.
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Matrix Property Detection
# ============================================================================

"""
    _is_orthogonal(A)

Check if matrix A is orthogonal: A^T * A = I.

Uses trig-aware simplification to handle composed rotations like Euler angles
where entries involve expressions like cos(α)cos(β) + sin(α)sin(β).
"""
function _is_orthogonal(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    ATA = transpose(A) * A
    for i in 1:n, j in 1:n
        target = i == j ? 1 : 0
        diff = ATA[i, j] - target
        if !_issymzero_trig(diff)
            return false
        end
    end
    return true
end

# Note: _is_unitary is defined in structure.jl and used here

"""
    _is_special_orthogonal(A)

Check if matrix A is in SO(n): orthogonal with det = +1.
Returns the dimension n if A is in SO(n), nothing otherwise.
"""
function _is_special_orthogonal(A)
    n = size(A, 1)
    _is_orthogonal(A) || return nothing
    
    d = det(A)
    if _issymzero_trig(d - 1)
        return n
    end
    return nothing
end

"""
    _is_special_unitary(A)

Check if matrix A is in SU(n): unitary with det = 1.
Returns the dimension n if A is in SU(n), nothing otherwise.
"""
function _is_special_unitary(A)
    n = size(A, 1)
    _is_unitary(A) || return nothing
    
    d = det(A)
    if _issymzero(Symbolics.simplify(d - 1))
        return n
    end
    return nothing
end

"""
    _symplectic_J(n)

Construct the standard symplectic form J for Sp(2n):
J = [0  I_n; -I_n  0]
"""
function _symplectic_J(n)
    J = zeros(Int, 2n, 2n)
    for i in 1:n
        J[i, n+i] = 1
        J[n+i, i] = -1
    end
    return J
end

"""
    _is_symplectic(A)

Check if matrix A is symplectic: A^T * J * A = J where J = [0 I; -I 0].
Returns the half-dimension n if A is in Sp(2n), nothing otherwise.
"""
function _is_symplectic(A)
    m = size(A, 1)
    size(A, 2) == m || return nothing
    m % 2 == 0 || return nothing
    
    n = div(m, 2)
    J = _symplectic_J(n)
    
    result = transpose(A) * J * A
    for i in 1:m, j in 1:m
        if !_issymzero(result[i, j] - J[i, j])
            return nothing
        end
    end
    return n
end

"""
    _has_complex_entries(mat)

Check if matrix contains complex (imaginary) components.
Used to disambiguate SU(n) from SO(n).
"""
function _has_complex_entries(mat)
    for entry in mat
        if entry isa Complex
            if !_issymzero(imag(entry))
                return true
            end
        elseif entry isa Num
            # Check if the symbolic expression contains 'im'
            str = string(entry)
            if occursin("im", str)
                return true
            end
        end
    end
    return false
end

"""
    _is_diagonal_matrix(mat)

Check if matrix is diagonal (all off-diagonal elements are zero).
"""
function _is_diagonal_matrix(mat)
    n = size(mat, 1)
    size(mat, 2) == n || return false
    
    for i in 1:n, j in 1:n
        if i != j
            entry = mat[i, j]
            if entry isa Complex
                if !_issymzero(real(entry)) || !_issymzero(imag(entry))
                    return false
                end
            elseif !_issymzero(entry)
                return false
            end
        end
    end
    return true
end

# ============================================================================
# Algebra-related checks
# ============================================================================

"""
    _is_skew_symmetric(A)

Check if A is skew-symmetric: A^T = -A.
Elements of so(n) in the fundamental representation are skew-symmetric.
"""
function _is_skew_symmetric(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    for i in 1:n, j in 1:n
        if !_issymzero(A[i, j] + A[j, i])
            return false
        end
    end
    return true
end

"""
    _is_skew_hermitian(A)

Check if A is skew-Hermitian: A† = -A.
Elements of su(n) in the fundamental representation are skew-Hermitian.
"""
function _is_skew_hermitian(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    for i in 1:n, j in 1:n
        diff = A[i, j] + conj(A[j, i])
        if diff isa Complex
            if !_issymzero(real(diff)) || !_issymzero(imag(diff))
                return false
            end
        elseif !_issymzero(diff)
            return false
        end
    end
    return true
end

"""
    _is_traceless(A)

Check if matrix A is traceless: tr(A) = 0.
"""
function _is_traceless(A)
    t = tr(A)
    if t isa Complex
        return _issymzero(real(t)) && _issymzero(imag(t))
    end
    return _issymzero(t)
end

# ============================================================================
# Master Lie Group Detection and Dispatch
# ============================================================================

"""
    _detect_lie_group(A)

Detect if matrix A belongs to a known Lie group with closed-form eigenvalues.

Returns (group_symbol, params) where group_symbol is one of:
:SO2, :SO3, :SO4, :SU2, :SU3, :Sp2, :Sp4, or nothing
"""
function _detect_lie_group(A)
    n = size(A, 1)
    size(A, 2) == n || return (nothing, nothing)
    
    # 2x2 matrices
    if n == 2
        so2_params = _is_SO2(A)
        if !isnothing(so2_params)
            return (:SO2, so2_params)
        end
        
        if _is_SU2(A)
            return (:SU2, nothing)
        end
        
        if _is_Sp2(A)
            return (:Sp2, nothing)
        end
    end
    
    # 3x3 matrices
    if n == 3
        if _is_SO3(A)
            return (:SO3, nothing)
        end
        
        if _is_SU3(A)
            return (:SU3, nothing)
        end
    end
    
    # 4x4 matrices
    if n == 4
        if _is_SO4(A)
            return (:SO4, nothing)
        end
        
        if _is_Sp4(A)
            return (:Sp4, nothing)
        end
    end
    
    return (nothing, nothing)
end

"""
    _lie_group_eigenvalues(A)

Compute eigenvalues using Lie group structure if detected.

Returns eigenvalues if A belongs to a supported Lie group, nothing otherwise.

Supports: SO(2), SO(3), SO(4), SU(2), SU(3), Sp(2), Sp(4)
"""
function _lie_group_eigenvalues(A)
    group, params = _detect_lie_group(A)
    
    isnothing(group) && return nothing
    
    if group == :SO2
        return _SO2_eigenvalues(A)
    elseif group == :SO3
        return _SO3_eigenvalues(A)
    elseif group == :SO4
        return _SO4_eigenvalues(A)
    elseif group == :SU2
        return _SU2_eigenvalues(A)
    elseif group == :SU3
        return _SU3_eigenvalues(A)
    elseif group == :Sp2
        return _Sp2_eigenvalues(A)
    elseif group == :Sp4
        return _Sp4_eigenvalues(A)
    end
    
    return nothing
end

# ============================================================================
# Master Lie Group Eigenpairs Dispatcher
# ============================================================================

"""
    _lie_group_eigenpairs(A)

Compute eigenvalue-eigenvector pairs using Lie group structure if detected.

Returns a vector of (eigenvalue, [eigenvectors]) tuples if A belongs to a 
supported Lie group with closed-form eigenvectors, nothing otherwise.

Currently supported (non-trivial cases only):
- SO(2): Fixed eigenvectors [1, ±i] for all rotation angles
- SO(3): Axis-aligned rotations (Rx, Ry, Rz) have known eigenvectors
- SO(4): Block-diagonal case [R₁ 0; 0 R₂] has known eigenvectors
- SU(2): Non-diagonal cases (Ux, Uy) have closed-form eigenvectors
- Sp(2): General 2×2 cases have closed-form eigenvectors
- Sp(4): Block-diagonal cases have known eigenvectors

NOT supported (trivial or requires nullspace):
- Diagonal matrices (eigenvalues/vectors are trivially on diagonal/standard basis)
- General rotations that don't match axis-aligned patterns

For complex cases, returns nothing and falls back to nullspace computation.
"""
function _lie_group_eigenpairs(A)
    group, params = _detect_lie_group(A)
    
    isnothing(group) && return nothing
    
    if group == :SO2
        return _SO2_eigenpairs(A)
    elseif group == :SO3
        return _SO3_eigenpairs(A)
    elseif group == :SO4
        return _SO4_eigenpairs(A)
    elseif group == :SU2
        return _SU2_eigenpairs(A)
    elseif group == :Sp2
        return _Sp2_eigenpairs(A)
    elseif group == :Sp4
        return _Sp4_eigenpairs(A)
    end
    # Note: SU3 eigenpairs not supported - diagonal is trivial, general needs nullspace
    
    return nothing
end
