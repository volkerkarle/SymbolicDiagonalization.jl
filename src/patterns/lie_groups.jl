# ============================================================================
# Lie Group Pattern Detection for Symbolic Diagonalization
# Focused on patterns that produce clean symbolic eigenvalues
# ============================================================================

# ============================================================================
# General Detection Functions (Symbolic-Only)
# ============================================================================

"""
    _is_orthogonal(A)

Check if matrix A is orthogonal: A^T * A = I.
Returns true if A is orthogonal, false otherwise.

Only uses symbolic zero checks - no numeric tolerance.
"""
function _is_orthogonal(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    ATA = transpose(A) * A
    for i in 1:n, j in 1:n
        target = i == j ? 1 : 0
        diff = ATA[i, j] - target
        if !_issymzero(diff)
            return false
        end
    end
    return true
end

"""
    _is_special_orthogonal(A)

Check if matrix A is in SO(n): orthogonal with det = +1.
Returns the dimension n if A is in SO(n), nothing otherwise.
"""
function _is_special_orthogonal(A)
    n = size(A, 1)
    _is_orthogonal(A) || return nothing
    
    d = det(A)
    if _issymzero(d - 1)
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
    if _issymzero(d - 1)
        return n
    end
    return nothing
end

"""
    _symplectic_j_matrix(n)

Construct the standard symplectic form J for Sp(2n):
J = [0  I_n; -I_n  0]
"""
function _symplectic_j_matrix(n)
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
    J = _symplectic_j_matrix(n)
    
    result = transpose(A) * J * A
    for i in 1:m, j in 1:m
        if !_issymzero(result[i, j] - J[i, j])
            return nothing
        end
    end
    return n
end

# ============================================================================
# SO(2) - 2D Rotations (Perfect Symbolic Support)
# ============================================================================

"""
    _is_so2(A)

Check if A is a 2x2 rotation matrix (SO(2)).
Returns (cos(theta), sin(theta)) if it is, nothing otherwise.

SO(2) matrices have the form:
[cos(theta)  -sin(theta)]
[sin(theta)   cos(theta)]
"""
function _is_so2(A)
    size(A) == (2, 2) || return nothing
    isnothing(_is_special_orthogonal(A)) && return nothing
    
    c = A[1, 1]  # cos(theta)
    s = A[2, 1]  # sin(theta)
    
    # Verify structure: A[1,2] = -sin(theta), A[2,2] = cos(theta)
    if !_issymzero(A[1, 2] + s) || !_issymzero(A[2, 2] - c)
        return nothing
    end
    
    return (c, s)
end

"""
    _so2_eigenvalues(A)

Compute eigenvalues of an SO(2) matrix (2D rotation).

For rotation by angle theta:
eigenvalues = e^(+/-i*theta) = cos(theta) +/- i*sin(theta)

This is the IDEAL case for symbolic diagonalization - eigenvalues
are expressed directly in terms of the matrix elements.
"""
function _so2_eigenvalues(A)
    cs = _is_so2(A)
    isnothing(cs) && return nothing
    c, s = cs
    
    # eigenvalues: cos(theta) +/- i*sin(theta)
    return [c + im*s, c - im*s]
end

# ============================================================================
# SO(3) - 3D Rotations
# ============================================================================

"""
    _is_so3(A)

Check if A is a 3x3 rotation matrix (SO(3)).
"""
function _is_so3(A)
    size(A) == (3, 3) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so3_eigenvalues(A)

Compute eigenvalues of an SO(3) matrix (3D rotation).

Using Rodrigues' formula:
- One eigenvalue is always 1 (the rotation axis)
- Other two are e^(+/-i*theta) where cos(theta) = (tr(A) - 1) / 2

For symbolic matrices, this produces:
- 1
- (tr(A)-1)/2 + i*sqrt(1 - ((tr(A)-1)/2)^2)
- (tr(A)-1)/2 - i*sqrt(1 - ((tr(A)-1)/2)^2)
"""
function _so3_eigenvalues(A)
    _is_so3(A) || return nothing
    
    tr_A = tr(A)
    cos_theta = (tr_A - 1) / 2
    
    # sin(theta) = sqrt(1 - cos^2(theta))
    # Use aggressive_simplify to get sin(theta) when possible
    sin_theta_sq = 1 - cos_theta^2
    sin_theta = aggressive_simplify(sqrt(sin_theta_sq))
    
    lambda1 = 1
    lambda2 = cos_theta + im * sin_theta
    lambda3 = cos_theta - im * sin_theta
    
    return [lambda1, simplify_eigenvalue(lambda2), simplify_eigenvalue(lambda3)]
end

# ============================================================================
# SO(4) - 4D Rotations
# ============================================================================

"""
    _is_so4(A)

Check if A is a 4x4 rotation matrix (SO(4)).
"""
function _is_so4(A)
    size(A) == (4, 4) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so4_eigenvalues(A)

Compute eigenvalues of an SO(4) matrix (4D rotation).

SO(4) has eigenvalues e^(+/-i*theta1), e^(+/-i*theta2) where the angles
can be extracted from trace invariants:
- tr(A) = 2(cos(theta1) + cos(theta2))
- tr(A^2) = 2(cos(2*theta1) + cos(2*theta2))

For symbolic matrices, this involves solving a quadratic for cos(theta_j).
"""
function _so4_eigenvalues(A)
    _is_so4(A) || return nothing
    
    t1 = tr(A)
    t2 = tr(A * A)
    
    # Let x = cos(theta1), y = cos(theta2)
    # t1 = 2(x + y), t2 = 4(x^2 + y^2) - 4
    # So: x + y = t1/2, x^2 + y^2 = (t2 + 4)/4
    # xy = [(t1/2)^2 - (t2+4)/4] / 2 = (t1^2 - t2 - 4) / 8
    
    sum_cos = t1 / 2
    prod_cos = (t1^2 - t2 - 4) / 8
    
    # x, y are roots of: z^2 - (sum_cos)*z + prod_cos = 0
    discriminant = sum_cos^2 - 4 * prod_cos
    sqrt_disc = sqrt(discriminant)
    
    cos_theta1 = (sum_cos + sqrt_disc) / 2
    cos_theta2 = (sum_cos - sqrt_disc) / 2
    
    # sin(theta) = sqrt(1 - cos^2(theta))
    # Use aggressive_simplify to get sin(theta) when possible
    sin_theta1 = aggressive_simplify(sqrt(1 - cos_theta1^2))
    sin_theta2 = aggressive_simplify(sqrt(1 - cos_theta2^2))
    
    lambda1 = cos_theta1 + im * sin_theta1
    lambda2 = cos_theta1 - im * sin_theta1
    lambda3 = cos_theta2 + im * sin_theta2
    lambda4 = cos_theta2 - im * sin_theta2
    
    return [simplify_eigenvalue(lambda1), simplify_eigenvalue(lambda2), 
            simplify_eigenvalue(lambda3), simplify_eigenvalue(lambda4)]
end

# ============================================================================
# SU(2) - Special Unitary 2x2
# ============================================================================

"""
    _is_su2(A)

Check if A is a 2x2 special unitary matrix (SU(2)).
"""
function _is_su2(A)
    size(A) == (2, 2) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _su2_eigenvalues(A)

Compute eigenvalues of an SU(2) matrix.

For SU(2), eigenvalues are e^(+/-i*theta) where cos(theta) = Re(tr(A))/2.
"""
function _su2_eigenvalues(A)
    _is_su2(A) || return nothing
    
    tr_A = tr(A)
    cos_theta = real(tr_A) / 2
    sin_theta = aggressive_simplify(sqrt(1 - cos_theta^2))
    
    return [simplify_eigenvalue(cos_theta + im * sin_theta), 
            simplify_eigenvalue(cos_theta - im * sin_theta)]
end

# ============================================================================
# SU(3) - Special Unitary 3x3
# ============================================================================

"""
    _is_su3(A)

Check if A is a 3x3 special unitary matrix (SU(3)).
"""
function _is_su3(A)
    size(A) == (3, 3) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _su3_eigenvalues(A)

Compute eigenvalues of an SU(3) matrix using the cubic formula.

For SU(3), eigenvalues are on the unit circle with product = 1.
The characteristic polynomial is: lambda^3 - tr(A)*lambda^2 + e2*lambda - 1 = 0
"""
function _su3_eigenvalues(A)
    _is_su3(A) || return nothing
    
    t1 = tr(A)
    t2 = tr(A * A)
    e2 = (t1^2 - t2) / 2
    
    # Characteristic polynomial: lambda^3 - t1*lambda^2 + e2*lambda - 1 = 0
    coeffs = [-1, e2, -t1, 1]
    return symbolic_roots(coeffs)
end

# ============================================================================
# Sp(2) - Symplectic 2x2 (equivalent to SL(2,R))
# ============================================================================

"""
    _is_sp2(A)

Check if A is a 2x2 symplectic matrix (Sp(2) = SL(2,R)).
"""
function _is_sp2(A)
    size(A) == (2, 2) || return false
    n = _is_symplectic(A)
    return !isnothing(n) && n == 1
end

"""
    _sp2_eigenvalues(A)

Compute eigenvalues of an Sp(2) matrix.

For Sp(2) with det = 1: lambda^2 - tr(A)*lambda + 1 = 0
"""
function _sp2_eigenvalues(A)
    _is_sp2(A) || return nothing
    
    t = tr(A)
    discriminant = t^2 - 4
    sqrt_disc = sqrt(discriminant)
    
    lambda1 = (t + sqrt_disc) / 2
    lambda2 = (t - sqrt_disc) / 2
    
    return [Symbolics.simplify(lambda1), Symbolics.simplify(lambda2)]
end

# ============================================================================
# Sp(4) - Symplectic 4x4
# ============================================================================

"""
    _is_sp4(A)

Check if A is a 4x4 symplectic matrix (Sp(4)).
"""
function _is_sp4(A)
    size(A) == (4, 4) || return false
    n = _is_symplectic(A)
    return !isnothing(n) && n == 2
end

"""
    _sp4_eigenvalues(A)

Compute eigenvalues of an Sp(4) matrix.

Eigenvalues come in reciprocal pairs. The characteristic polynomial
is palindromic: lambda^4 - a*lambda^3 + b*lambda^2 - a*lambda + 1 = 0
This reduces to a quadratic in mu = lambda + 1/lambda.
"""
function _sp4_eigenvalues(A)
    _is_sp4(A) || return nothing
    
    t1 = tr(A)
    t2 = tr(A * A)
    
    c3 = t1
    c2 = (t1^2 - t2) / 2
    
    # Solve mu^2 - c3*mu + (c2 - 2) = 0
    a_coeff = c3
    b_coeff = c2 - 2
    
    disc_mu = a_coeff^2 - 4 * b_coeff
    sqrt_disc_mu = sqrt(disc_mu)
    
    mu1 = (a_coeff + sqrt_disc_mu) / 2
    mu2 = (a_coeff - sqrt_disc_mu) / 2
    
    # For each mu, solve lambda^2 - mu*lambda + 1 = 0
    disc1 = mu1^2 - 4
    sqrt_disc1 = sqrt(disc1)
    lambda1 = (mu1 + sqrt_disc1) / 2
    lambda2 = (mu1 - sqrt_disc1) / 2
    
    disc2 = mu2^2 - 4
    sqrt_disc2 = sqrt(disc2)
    lambda3 = (mu2 + sqrt_disc2) / 2
    lambda4 = (mu2 - sqrt_disc2) / 2
    
    return [Symbolics.simplify(lambda1), Symbolics.simplify(lambda2),
            Symbolics.simplify(lambda3), Symbolics.simplify(lambda4)]
end

# ============================================================================
# Master Detection and Dispatch Function
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
        so2_params = _is_so2(A)
        if !isnothing(so2_params)
            return (:SO2, so2_params)
        end
        
        if _is_su2(A)
            return (:SU2, nothing)
        end
        
        if _is_sp2(A)
            return (:Sp2, nothing)
        end
    end
    
    # 3x3 matrices
    if n == 3
        if _is_so3(A)
            return (:SO3, nothing)
        end
        
        if _is_su3(A)
            return (:SU3, nothing)
        end
    end
    
    # 4x4 matrices
    if n == 4
        if _is_so4(A)
            return (:SO4, nothing)
        end
        
        if _is_sp4(A)
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
        return _so2_eigenvalues(A)
    elseif group == :SO3
        return _so3_eigenvalues(A)
    elseif group == :SO4
        return _so4_eigenvalues(A)
    elseif group == :SU2
        return _su2_eigenvalues(A)
    elseif group == :SU3
        return _su3_eigenvalues(A)
    elseif group == :Sp2
        return _sp2_eigenvalues(A)
    elseif group == :Sp4
        return _sp4_eigenvalues(A)
    end
    
    return nothing
end
