# ============================================================================
# Lie Group Pattern Detection and Eigenvalue Computation
# Covers SO(n), SU(n), Sp(2n), O(p,q) and related matrix Lie groups
# ============================================================================

# ============================================================================
# Helper Functions for Numeric/Symbolic Compatibility
# ============================================================================

"""
    _safe_simplify(x)

Apply Symbolics.simplify only to symbolic expressions.
For numeric values, return as-is.
"""
function _safe_simplify(x)
    if x isa Num || x isa Complex{Num}
        return Symbolics.simplify(x)
    end
    return x
end

"""
    _safe_sqrt_unit_circle(c)

Compute sin(θ) = √(1 - cos²θ) safely for both symbolic and numeric values.
Handles floating-point errors where cos²θ might be slightly > 1.
"""
function _safe_sqrt_unit_circle(c)
    c_sq = c^2
    if c_sq isa Number && !(c_sq isa Num)
        # Numeric: clamp to handle floating point errors
        arg = max(0.0, 1 - real(c_sq))
        return sqrt(Complex(arg))
    else
        # Symbolic: use symbolic sqrt
        return sqrt(1 - c_sq)
    end
end

"""
    _make_numeric_complex(coeffs)

For numeric coefficient arrays, ensure they're Complex to handle
cases where intermediate results might be complex.
"""
function _make_numeric_complex(coeffs)
    if all(c -> c isa Number && !(c isa Num), coeffs)
        return Complex{Float64}.(coeffs)
    end
    return coeffs
end

# ============================================================================
# General Detection Functions
# ============================================================================

"""
    _is_orthogonal(A)

Check if matrix A is orthogonal: A^T * A = I.
Returns true if A is orthogonal, false otherwise.

Handles both symbolic and numeric matrices with appropriate tolerance.
"""
function _is_orthogonal(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    ATA = transpose(A) * A
    for i in 1:n, j in 1:n
        target = i == j ? 1 : 0
        diff = ATA[i, j] - target
        # Handle numeric vs symbolic differently
        if diff isa Num || diff isa Complex{Num}
            if !_issymzero(diff)
                return false
            end
        elseif diff isa Number
            if !isapprox(diff, 0, atol=1e-10)
                return false
            end
        else
            if !_issymzero(diff)
                return false
            end
        end
    end
    return true
end

"""
    _is_special_orthogonal(A)

Check if matrix A is in SO(n): orthogonal with det = +1.
Returns the dimension n if A is in SO(n), nothing otherwise.

Handles both symbolic and numeric matrices.
"""
function _is_special_orthogonal(A)
    n = size(A, 1)
    _is_orthogonal(A) || return nothing
    
    # Check det = 1
    d = det(A)
    if d isa Num || d isa Complex{Num}
        if _issymzero(d - 1)
            return n
        end
    elseif d isa Number
        if isapprox(d, 1, atol=1e-10)
            return n
        end
    else
        if _issymzero(d - 1)
            return n
        end
    end
    return nothing
end

"""
    _is_special_unitary(A)

Check if matrix A is in SU(n): unitary with det = 1.
Returns the dimension n if A is in SU(n), nothing otherwise.

Note: _is_unitary is already defined in structure.jl
"""
function _is_special_unitary(A)
    n = size(A, 1)
    _is_unitary(A) || return nothing
    
    # Check det = 1
    d = det(A)
    # For complex matrices, det should be 1 (real part 1, imag part 0)
    if _issymzero(d - 1)
        return n
    end
    return nothing
end

"""
    _symplectic_j_matrix(n)

Construct the standard symplectic form J for Sp(2n):
J = [0  I_n; -I_n  0]

This is the 2n × 2n matrix representing the symplectic form.
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

Note: The symplectic group Sp(2n) consists of 2n × 2n matrices.
"""
function _is_symplectic(A)
    m = size(A, 1)
    size(A, 2) == m || return nothing
    m % 2 == 0 || return nothing  # Must be even dimension
    
    n = div(m, 2)
    J = _symplectic_j_matrix(n)
    
    # Check A^T * J * A = J
    result = transpose(A) * J * A
    for i in 1:m, j in 1:m
        if !_issymzero(result[i, j] - J[i, j])
            return nothing
        end
    end
    return n
end

"""
    _is_compact_symplectic(A)

Check if matrix A is in the compact symplectic group Sp(2n) ∩ U(2n).
These are unitary symplectic matrices with eigenvalues in conjugate pairs on the unit circle.

Returns the half-dimension n if A is compact symplectic, nothing otherwise.
"""
function _is_compact_symplectic(A)
    n = _is_symplectic(A)
    isnothing(n) && return nothing
    _is_unitary(A) || return nothing
    return n
end

"""
    _indefinite_metric(p, q)

Construct the indefinite metric η = diag(+1,...,+1,-1,...,-1) for O(p,q).
Has p positive and q negative entries on the diagonal.
"""
function _indefinite_metric(p, q)
    n = p + q
    η = zeros(Int, n, n)
    for i in 1:p
        η[i, i] = 1
    end
    for i in 1:q
        η[p+i, p+i] = -1
    end
    return η
end

"""
    _is_indefinite_orthogonal(A, p, q)

Check if matrix A is in O(p,q): preserves the indefinite metric η = diag(+1,...,-1,...).
A^T * η * A = η

Returns true if A is in O(p,q), false otherwise.
"""
function _is_indefinite_orthogonal(A, p, q)
    n = p + q
    size(A) == (n, n) || return false
    
    η = _indefinite_metric(p, q)
    result = transpose(A) * η * A
    
    for i in 1:n, j in 1:n
        diff = result[i, j] - η[i, j]
        # Check symbolic types first (Num and Complex{Num}) since Num <: Real <: Number
        if diff isa Num || diff isa Complex{Num}
            if !_issymzero(diff)
                return false
            end
        elseif diff isa Number
            if !isapprox(diff, 0, atol=1e-10)
                return false
            end
        else
            # Unknown type, try symbolic check
            if !_issymzero(diff)
                return false
            end
        end
    end
    return true
end

"""
    _is_lorentz_group(A)

Check if matrix A is in the Lorentz group SO(1,3) or more generally SO(1,n-1).
These preserve the Minkowski metric η = diag(+1,-1,-1,-1).

Returns (p, q) signature if A is in SO(p,q) with p=1, nothing otherwise.
"""
function _is_lorentz_group(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing
    
    # Check for SO(1, n-1) - Lorentz signature
    if _is_indefinite_orthogonal(A, 1, n-1)
        # Check det = +1 for SO (not just O)
        d = det(A)
        if _issymzero(d - 1)
            return (1, n-1)
        end
    end
    return nothing
end

# ============================================================================
# SO(2) - 2D Rotations
# ============================================================================

"""
    _is_so2(A)

Check if A is a 2×2 rotation matrix (SO(2)).
Returns the rotation angle θ if it is, nothing otherwise.

SO(2) matrices have the form:
[cos(θ)  -sin(θ)]
[sin(θ)   cos(θ)]
"""
function _is_so2(A)
    size(A) == (2, 2) || return nothing
    isnothing(_is_special_orthogonal(A)) && return nothing
    
    # Extract angle: θ = atan(A[2,1], A[1,1])
    # For symbolic matrices, we can't always simplify atan
    # Just verify structure and return the components
    c = A[1, 1]  # cos(θ)
    s = A[2, 1]  # sin(θ)
    
    # Verify structure: A[1,2] = -sin(θ), A[2,2] = cos(θ)
    if !_issymzero(A[1, 2] + s) || !_issymzero(A[2, 2] - c)
        return nothing
    end
    
    return (c, s)  # Return (cos θ, sin θ)
end

"""
    _so2_eigenvalues(A)

Compute eigenvalues of an SO(2) matrix (2D rotation).

For rotation by angle θ:
eigenvalues = e^(±iθ) = cos(θ) ± i·sin(θ)

Returns vector of 2 eigenvalues.
"""
function _so2_eigenvalues(A)
    cs = _is_so2(A)
    isnothing(cs) && return nothing
    c, s = cs
    
    # eigenvalues: cos(θ) ± i*sin(θ)
    return [c + im*s, c - im*s]
end

# ============================================================================
# SO(3) - 3D Rotations
# ============================================================================

"""
    _is_so3(A)

Check if A is a 3×3 rotation matrix (SO(3)).
Returns true if it is, false otherwise.
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
- Other two are e^(±iθ) where cos(θ) = (tr(A) - 1) / 2

Returns vector of 3 eigenvalues: [1, e^(iθ), e^(-iθ)]
"""
function _so3_eigenvalues(A)
    _is_so3(A) || return nothing
    
    # cos(θ) = (tr(A) - 1) / 2
    tr_A = tr(A)
    cos_θ = (tr_A - 1) / 2
    
    # sin(θ) = ±√(1 - cos²(θ))
    # We use the positive root by convention
    sin_θ_squared = 1 - cos_θ^2
    sin_θ = sqrt(sin_θ_squared)
    
    # eigenvalues: 1, cos(θ) + i*sin(θ), cos(θ) - i*sin(θ)
    λ1 = 1
    λ2 = cos_θ + im * sin_θ
    λ3 = cos_θ - im * sin_θ
    
    return [λ1, Symbolics.simplify(λ2), Symbolics.simplify(λ3)]
end

# ============================================================================
# SO(4) - 4D Rotations (via quaternion decomposition)
# ============================================================================

"""
    _is_so4(A)

Check if A is a 4×4 rotation matrix (SO(4)).
Returns true if it is, false otherwise.
"""
function _is_so4(A)
    size(A) == (4, 4) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so4_eigenvalues(A)

Compute eigenvalues of an SO(4) matrix (4D rotation).

SO(4) is locally isomorphic to SO(3) × SO(3) (double cover via quaternions).
Every SO(4) matrix can be written as a product of two commuting rotations
in orthogonal 2-planes, giving eigenvalues:
    e^(±iθ₁), e^(±iθ₂)

The angles can be extracted from:
- tr(A) = 2(cos(θ₁) + cos(θ₂))
- det(A - I) relates to the product of sines

For general SO(4), we use the characteristic polynomial which is biquadratic.
"""
function _so4_eigenvalues(A)
    _is_so4(A) || return nothing
    
    # For SO(4), the characteristic polynomial has special structure
    # λ⁴ - tr(A)λ³ + ... 
    # But eigenvalues come in conjugate pairs on unit circle
    
    # Use the fact that for orthogonal matrices, if λ is eigenvalue, so is 1/λ̄ = λ̄ (on unit circle)
    # So eigenvalues are e^(±iθ₁), e^(±iθ₂)
    
    # Method: Use invariants
    # tr(A) = 2(cos θ₁ + cos θ₂)
    # For SO(4), there's a second invariant from the Pfaffian or from tr(A²)
    # tr(A²) = 2(cos 2θ₁ + cos 2θ₂) = 2(2cos²θ₁ - 1 + 2cos²θ₂ - 1)
    #        = 4(cos²θ₁ + cos²θ₂) - 4
    
    t1 = tr(A)          # = 2(cos θ₁ + cos θ₂)
    t2 = tr(A * A)      # = 4(cos²θ₁ + cos²θ₂) - 4
    
    # Let x = cos θ₁, y = cos θ₂
    # t1 = 2(x + y), t2 = 4(x² + y²) - 4
    # So: x + y = t1/2, x² + y² = (t2 + 4)/4
    # (x + y)² = x² + 2xy + y² → xy = [(t1/2)² - (t2+4)/4] / 2 = (t1² - t2 - 4) / 8
    
    sum_cos = t1 / 2
    prod_cos = (t1^2 - t2 - 4) / 8
    
    # x, y are roots of: z² - (sum_cos)z + prod_cos = 0
    discriminant = sum_cos^2 - 4 * prod_cos
    
    # Handle numeric case where discriminant might be slightly negative
    if discriminant isa Number && !(discriminant isa Num)
        sqrt_disc = sqrt(Complex(discriminant))
    else
        sqrt_disc = sqrt(discriminant)
    end
    
    cos_θ1 = (sum_cos + sqrt_disc) / 2
    cos_θ2 = (sum_cos - sqrt_disc) / 2
    
    # sin θ = √(1 - cos²θ)
    sin_θ1 = _safe_sqrt_unit_circle(cos_θ1)
    sin_θ2 = _safe_sqrt_unit_circle(cos_θ2)
    
    # Eigenvalues: e^(±iθ₁), e^(±iθ₂)
    λ1 = cos_θ1 + im * sin_θ1
    λ2 = cos_θ1 - im * sin_θ1
    λ3 = cos_θ2 + im * sin_θ2
    λ4 = cos_θ2 - im * sin_θ2
    
    return [_safe_simplify(λ1), _safe_simplify(λ2), 
            _safe_simplify(λ3), _safe_simplify(λ4)]
end

# ============================================================================
# SO(5) - 5D Rotations
# ============================================================================

"""
    _is_so5(A)

Check if A is a 5×5 rotation matrix (SO(5)).
Returns true if it is, false otherwise.
"""
function _is_so5(A)
    size(A) == (5, 5) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so5_eigenvalues(A)

Compute eigenvalues of an SO(5) matrix (5D rotation).

For SO(5) (odd dimension), eigenvalues are:
- One eigenvalue is always 1 (the rotation axis in 5D)
- Two conjugate pairs: e^(±iθ₁), e^(±iθ₂)

The angles are extracted from trace invariants:
- tr(A) = 1 + 2(cos θ₁ + cos θ₂)
- tr(A²) = 1 + 2(cos 2θ₁ + cos 2θ₂) = 1 + 4(cos²θ₁ + cos²θ₂) - 4

Returns vector of 5 eigenvalues: [1, e^(±iθ₁), e^(±iθ₂)]
"""
function _so5_eigenvalues(A)
    _is_so5(A) || return nothing
    
    # For SO(5): eigenvalues are 1, e^(±iθ₁), e^(±iθ₂)
    # tr(A) = 1 + 2(cos θ₁ + cos θ₂)
    # tr(A²) = 1 + 2(cos 2θ₁ + cos 2θ₂) = 1 + 4(cos²θ₁ + cos²θ₂) - 4
    
    t1 = tr(A)
    t2 = tr(A * A)
    
    # Let x = cos θ₁, y = cos θ₂
    # t1 = 1 + 2(x + y) → x + y = (t1 - 1) / 2
    # t2 = 1 + 4(x² + y²) - 4 = 4(x² + y²) - 3 → x² + y² = (t2 + 3) / 4
    
    sum_cos = (t1 - 1) / 2
    sum_cos_sq = (t2 + 3) / 4
    
    # xy = [(x+y)² - (x² + y²)] / 2 = [sum_cos² - sum_cos_sq] / 2
    prod_cos = (sum_cos^2 - sum_cos_sq) / 2
    
    # x, y are roots of: z² - (sum_cos)z + prod_cos = 0
    discriminant = sum_cos^2 - 4 * prod_cos
    discriminant = _safe_simplify(discriminant)
    
    # Handle numeric case where discriminant might be slightly negative
    if discriminant isa Number && !(discriminant isa Num)
        sqrt_disc = sqrt(Complex(discriminant))
    else
        sqrt_disc = sqrt(discriminant)
    end
    
    cos_θ1 = (sum_cos + sqrt_disc) / 2
    cos_θ2 = (sum_cos - sqrt_disc) / 2
    
    # sin θ = √(1 - cos²θ)
    sin_θ1 = _safe_sqrt_unit_circle(cos_θ1)
    sin_θ2 = _safe_sqrt_unit_circle(cos_θ2)
    
    # Eigenvalues
    λ1 = one(eltype(A)) + zero(eltype(A))*im
    λ2 = cos_θ1 + im * sin_θ1
    λ3 = cos_θ1 - im * sin_θ1
    λ4 = cos_θ2 + im * sin_θ2
    λ5 = cos_θ2 - im * sin_θ2
    
    return [λ1, _safe_simplify(λ2), _safe_simplify(λ3),
            _safe_simplify(λ4), _safe_simplify(λ5)]
end

# ============================================================================
# SO(6) - 6D Rotations
# ============================================================================

"""
    _is_so6(A)

Check if A is a 6×6 rotation matrix (SO(6)).
Returns true if it is, false otherwise.
"""
function _is_so6(A)
    size(A) == (6, 6) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so6_eigenvalues(A)

Compute eigenvalues of an SO(6) matrix (6D rotation).

For SO(6) (even dimension), eigenvalues are three conjugate pairs:
- e^(±iθ₁), e^(±iθ₂), e^(±iθ₃)

The angles are extracted from trace invariants using Newton's identities.
This requires solving a cubic equation for cos θⱼ.

Returns vector of 6 eigenvalues.
"""
function _so6_eigenvalues(A)
    _is_so6(A) || return nothing
    
    # For SO(6): eigenvalues are e^(±iθ₁), e^(±iθ₂), e^(±iθ₃)
    # tr(A) = 2(cos θ₁ + cos θ₂ + cos θ₃)
    # tr(A²) = 2(cos 2θ₁ + cos 2θ₂ + cos 2θ₃) = 4(cos²θ₁ + cos²θ₂ + cos²θ₃) - 6
    # tr(A³) = 2(cos 3θ₁ + cos 3θ₂ + cos 3θ₃)
    #        = 2(4cos³θⱼ - 3cosθⱼ) = 8(cos³θ₁ + cos³θ₂ + cos³θ₃) - 6(cos θ₁ + cos θ₂ + cos θ₃)
    
    t1 = tr(A)
    t2 = tr(A * A)
    t3 = tr(A * A * A)
    
    # Let p₁ = cos θ₁ + cos θ₂ + cos θ₃ (power sum)
    # Let p₂ = cos²θ₁ + cos²θ₂ + cos²θ₃
    # Let p₃ = cos³θ₁ + cos³θ₂ + cos³θ₃
    
    p1 = t1 / 2
    p2 = (t2 + 6) / 4
    p3 = (t3 + 6 * p1) / 8  # From cos 3θ = 4cos³θ - 3cosθ
    
    # Newton's identities: e₁ = p₁, e₂ = (p₁² - p₂)/2, e₃ = (p₁³ - 3p₁p₂ + 2p₃)/6
    # where e₁, e₂, e₃ are elementary symmetric polynomials
    e1 = p1
    e2 = (p1^2 - p2) / 2
    e3 = (p1^3 - 3*p1*p2 + 2*p3) / 6
    
    # Simplify coefficients (safe for numeric)
    e1 = _safe_simplify(e1)
    e2 = _safe_simplify(e2)
    e3 = _safe_simplify(e3)
    
    # cos θⱼ are roots of: z³ - e₁z² + e₂z - e₃ = 0
    # Coefficients in ascending order: [-e₃, e₂, -e₁, 1]
    coeffs = _make_numeric_complex([-e3, e2, -e1, 1])
    cos_vals = symbolic_roots(coeffs)
    
    # Build eigenvalues from cos values
    eigenvalues = []
    for c in cos_vals
        s = _safe_sqrt_unit_circle(c)
        push!(eigenvalues, _safe_simplify(c + im * s))
        push!(eigenvalues, _safe_simplify(c - im * s))
    end
    
    return eigenvalues
end

# ============================================================================
# SO(7) - 7D Rotations
# ============================================================================

"""
    _is_so7(A)

Check if A is a 7×7 rotation matrix (SO(7)).
Returns true if it is, false otherwise.
"""
function _is_so7(A)
    size(A) == (7, 7) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so7_eigenvalues(A)

Compute eigenvalues of an SO(7) matrix (7D rotation).

For SO(7) (odd dimension), eigenvalues are:
- One eigenvalue is always 1
- Three conjugate pairs: e^(±iθ₁), e^(±iθ₂), e^(±iθ₃)

Uses cubic formula for the three cos θⱼ values.

Returns vector of 7 eigenvalues.
"""
function _so7_eigenvalues(A)
    _is_so7(A) || return nothing
    
    # For SO(7): eigenvalues are 1, e^(±iθ₁), e^(±iθ₂), e^(±iθ₃)
    # tr(A) = 1 + 2(cos θ₁ + cos θ₂ + cos θ₃)
    # tr(A²) = 1 + 4(cos²θ₁ + cos²θ₂ + cos²θ₃) - 6
    # tr(A³) = 1 + 8(cos³θ₁ + cos³θ₂ + cos³θ₃) - 6(cos θ₁ + cos θ₂ + cos θ₃)
    
    t1 = tr(A)
    t2 = tr(A * A)
    t3 = tr(A * A * A)
    
    # Power sums adjusted for the eigenvalue 1
    p1 = (t1 - 1) / 2
    p2 = (t2 - 1 + 6) / 4  # = (t2 + 5) / 4
    p3 = (t3 - 1 + 6 * p1) / 8
    
    # Newton's identities
    e1 = p1
    e2 = (p1^2 - p2) / 2
    e3 = (p1^3 - 3*p1*p2 + 2*p3) / 6
    
    e1 = _safe_simplify(e1)
    e2 = _safe_simplify(e2)
    e3 = _safe_simplify(e3)
    
    # cos θⱼ are roots of: z³ - e₁z² + e₂z - e₃ = 0
    coeffs = _make_numeric_complex([-e3, e2, -e1, 1])
    cos_vals = symbolic_roots(coeffs)
    
    # Build eigenvalues
    eigenvalues = [one(eltype(A)) + zero(eltype(A))*im]  # 1 as first eigenvalue
    for c in cos_vals
        s = _safe_sqrt_unit_circle(c)
        push!(eigenvalues, _safe_simplify(c + im * s))
        push!(eigenvalues, _safe_simplify(c - im * s))
    end
    
    return eigenvalues
end

# ============================================================================
# SO(8) - 8D Rotations
# ============================================================================

"""
    _is_so8(A)

Check if A is an 8×8 rotation matrix (SO(8)).
Returns true if it is, false otherwise.
"""
function _is_so8(A)
    size(A) == (8, 8) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so8_eigenvalues(A)

Compute eigenvalues of an SO(8) matrix (8D rotation).

For SO(8) (even dimension), eigenvalues are four conjugate pairs:
- e^(±iθ₁), e^(±iθ₂), e^(±iθ₃), e^(±iθ₄)

Uses quartic formula for the four cos θⱼ values.

Returns vector of 8 eigenvalues.
"""
function _so8_eigenvalues(A)
    _is_so8(A) || return nothing
    
    # For SO(8): eigenvalues are e^(±iθⱼ) for j=1,2,3,4
    # tr(Aᵏ) = 2 Σⱼ cos(k·θⱼ)
    
    t1 = tr(A)
    t2 = tr(A * A)
    t3 = tr(A * A * A)
    t4 = tr(A * A * A * A)
    
    # Power sums: pₖ = Σⱼ cos^k(θⱼ)
    # From cos(kθ) = Chebyshev polynomial in cos(θ)
    # cos(2θ) = 2cos²θ - 1
    # cos(3θ) = 4cos³θ - 3cosθ  
    # cos(4θ) = 8cos⁴θ - 8cos²θ + 1
    
    p1 = t1 / 2  # Σ cos θⱼ
    p2 = (t2 + 8) / 4  # Σ cos²θⱼ = (t2/2 + 4)/2 = (t2 + 8)/4
    p3 = (t3 + 6 * p1) / 8  # From 2cos(3θ) = 8cos³θ - 6cosθ
    p4 = (t4 - 8 + 16 * p2) / 16  # From 2cos(4θ) = 16cos⁴θ - 16cos²θ + 2
    
    # Newton's identities for elementary symmetric polynomials
    e1 = p1
    e2 = (p1^2 - p2) / 2
    e3 = (p1^3 - 3*p1*p2 + 2*p3) / 6
    e4 = (p1^4 - 6*p1^2*p2 + 3*p2^2 + 8*p1*p3 - 6*p4) / 24
    
    e1 = _safe_simplify(e1)
    e2 = _safe_simplify(e2)
    e3 = _safe_simplify(e3)
    e4 = _safe_simplify(e4)
    
    # cos θⱼ are roots of: z⁴ - e₁z³ + e₂z² - e₃z + e₄ = 0
    coeffs = _make_numeric_complex([e4, -e3, e2, -e1, 1])
    cos_vals = symbolic_roots(coeffs)
    
    # Build eigenvalues
    eigenvalues = []
    for c in cos_vals
        s = _safe_sqrt_unit_circle(c)
        push!(eigenvalues, _safe_simplify(c + im * s))
        push!(eigenvalues, _safe_simplify(c - im * s))
    end
    
    return eigenvalues
end

# ============================================================================
# SO(9) - 9D Rotations
# ============================================================================

"""
    _is_so9(A)

Check if A is a 9×9 rotation matrix (SO(9)).
Returns true if it is, false otherwise.
"""
function _is_so9(A)
    size(A) == (9, 9) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _so9_eigenvalues(A)

Compute eigenvalues of an SO(9) matrix (9D rotation).

For SO(9) (odd dimension), eigenvalues are:
- One eigenvalue is always 1
- Four conjugate pairs: e^(±iθ₁), e^(±iθ₂), e^(±iθ₃), e^(±iθ₄)

Uses quartic formula for the four cos θⱼ values.

Returns vector of 9 eigenvalues.
"""
function _so9_eigenvalues(A)
    _is_so9(A) || return nothing
    
    # For SO(9): eigenvalues are 1, e^(±iθⱼ) for j=1,2,3,4
    
    t1 = tr(A)
    t2 = tr(A * A)
    t3 = tr(A * A * A)
    t4 = tr(A * A * A * A)
    
    # Power sums adjusted for eigenvalue 1
    p1 = (t1 - 1) / 2
    p2 = (t2 - 1 + 8) / 4  # = (t2 + 7) / 4
    p3 = (t3 - 1 + 6 * p1) / 8
    p4 = (t4 - 1 - 8 + 16 * p2) / 16  # = (t4 + 16*p2 - 9) / 16
    
    # Newton's identities
    e1 = p1
    e2 = (p1^2 - p2) / 2
    e3 = (p1^3 - 3*p1*p2 + 2*p3) / 6
    e4 = (p1^4 - 6*p1^2*p2 + 3*p2^2 + 8*p1*p3 - 6*p4) / 24
    
    e1 = _safe_simplify(e1)
    e2 = _safe_simplify(e2)
    e3 = _safe_simplify(e3)
    e4 = _safe_simplify(e4)
    
    # cos θⱼ are roots of: z⁴ - e₁z³ + e₂z² - e₃z + e₄ = 0
    coeffs = _make_numeric_complex([e4, -e3, e2, -e1, 1])
    cos_vals = symbolic_roots(coeffs)
    
    # Build eigenvalues
    eigenvalues = [one(eltype(A)) + zero(eltype(A))*im]  # 1 as first eigenvalue
    for c in cos_vals
        s = _safe_sqrt_unit_circle(c)
        push!(eigenvalues, _safe_simplify(c + im * s))
        push!(eigenvalues, _safe_simplify(c - im * s))
    end
    
    return eigenvalues
end

# ============================================================================
# General SO(n) eigenvalues for n ≤ 9
# ============================================================================

"""
    _general_so_eigenvalues(A)

Compute eigenvalues for any SO(n) matrix with n ≤ 9.

For n ≤ 9, we can solve for the rotation angles symbolically:
- n = 2k (even): k pairs e^(±iθⱼ) → degree k polynomial in cos θⱼ
- n = 2k+1 (odd): eigenvalue 1 + k pairs → same degree k polynomial

The solvability limits:
- n ≤ 3: linear (trivial)
- n ≤ 5: quadratic formula
- n ≤ 7: cubic formula (Cardano)
- n ≤ 9: quartic formula (Ferrari)
- n ≥ 10: degree 5+ polynomial, not solvable in radicals (Abel-Ruffini)

Returns eigenvalues if n ≤ 9 and A is in SO(n), nothing otherwise.
"""
function _general_so_eigenvalues(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    # Dispatch to specific implementations
    if n == 2
        return _so2_eigenvalues(A)
    elseif n == 3
        return _so3_eigenvalues(A)
    elseif n == 4
        return _so4_eigenvalues(A)
    elseif n == 5
        return _so5_eigenvalues(A)
    elseif n == 6
        return _so6_eigenvalues(A)
    elseif n == 7
        return _so7_eigenvalues(A)
    elseif n == 8
        return _so8_eigenvalues(A)
    elseif n == 9
        return _so9_eigenvalues(A)
    end
    
    # n ≥ 10: Cannot solve in radicals
    return nothing
end

# ============================================================================
# SU(2) - Special Unitary 2×2
# ============================================================================

"""
    _is_su2(A)

Check if A is a 2×2 special unitary matrix (SU(2)).
Returns true if it is, false otherwise.

SU(2) matrices have the form:
[α   -β̄]
[β    ᾱ]
where |α|² + |β|² = 1
"""
function _is_su2(A)
    size(A) == (2, 2) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _su2_eigenvalues(A)

Compute eigenvalues of an SU(2) matrix.

For SU(2), eigenvalues are e^(±iθ) where cos(θ) = Re(tr(A))/2.
Since det = 1 and eigenvalues are e^(iθ₁), e^(iθ₂), we have θ₁ + θ₂ = 0.

Returns vector of 2 eigenvalues.
"""
function _su2_eigenvalues(A)
    _is_su2(A) || return nothing
    
    # For SU(2): cos(θ) = Re(tr(A)) / 2
    tr_A = tr(A)
    cos_θ = real(tr_A) / 2
    
    # sin²(θ) = 1 - cos²(θ)
    sin_θ = sqrt(1 - cos_θ^2)
    
    # eigenvalues: e^(±iθ)
    return [cos_θ + im * sin_θ, cos_θ - im * sin_θ]
end

# ============================================================================
# SU(3) - Special Unitary 3×3 (using Cardano's formula)
# ============================================================================

"""
    _is_su3(A)

Check if A is a 3×3 special unitary matrix (SU(3)).
Returns true if it is, false otherwise.
"""
function _is_su3(A)
    size(A) == (3, 3) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _su3_eigenvalues(A)

Compute eigenvalues of an SU(3) matrix.

For SU(3), eigenvalues are e^(iθ₁), e^(iθ₂), e^(iθ₃) where θ₁ + θ₂ + θ₃ = 0 (mod 2π).
The characteristic polynomial is cubic, solvable by Cardano's formula.

Since eigenvalues are on the unit circle and their product is 1 (det = 1),
we can use the trace and tr(A²) to find them.
"""
function _su3_eigenvalues(A)
    _is_su3(A) || return nothing
    
    # For unitary matrices, eigenvalues are on unit circle
    # For SU(3), product of eigenvalues = det(A) = 1
    # Let eigenvalues be e^(iθ₁), e^(iθ₂), e^(iθ₃)
    
    # Use Newton's identities:
    # p₁ = e₁ = tr(A) = λ₁ + λ₂ + λ₃
    # p₂ = tr(A²) = λ₁² + λ₂² + λ₃²
    # e₂ = (p₁² - p₂)/2 = λ₁λ₂ + λ₁λ₃ + λ₂λ₃
    # e₃ = det(A) = λ₁λ₂λ₃ = 1
    
    # The characteristic polynomial is: λ³ - e₁λ² + e₂λ - e₃ = 0
    # i.e., λ³ - tr(A)λ² + e₂λ - 1 = 0
    
    t1 = tr(A)
    t2 = tr(A * A)
    e2 = (t1^2 - t2) / 2
    
    # Coefficients for λ³ + aλ² + bλ + c = 0 (monic form)
    # Here: λ³ - t1·λ² + e2·λ - 1 = 0
    # So a = -t1, b = e2, c = -1
    
    # Use the existing cubic solver via symbolic_roots
    # Coefficients in ascending order: [c, b, a, 1] = [-1, e2, -t1, 1]
    coeffs = [-1, e2, -t1, 1]
    
    # Solve using the package's cubic formula
    vals = symbolic_roots(coeffs)
    
    return vals
end

# ============================================================================
# Sp(2) - Symplectic 2×2 (equivalent to SL(2,R))
# ============================================================================

"""
    _is_sp2(A)

Check if A is a 2×2 symplectic matrix (Sp(2) ≅ SL(2,R)).
Note: Sp(2) means Sp(2·1) = Sp(2), i.e., 2×2 matrices preserving the form.

For 2×2, the symplectic condition A^T J A = J with J = [0 1; -1 0]
is equivalent to det(A) = 1, so Sp(2) ≅ SL(2,R).

Returns true if it is, false otherwise.
"""
function _is_sp2(A)
    size(A) == (2, 2) || return false
    n = _is_symplectic(A)
    return !isnothing(n) && n == 1
end

"""
    _sp2_eigenvalues(A)

Compute eigenvalues of an Sp(2) matrix (2×2 symplectic).

For Sp(2), eigenvalues satisfy λ · (1/λ) = 1 (come in reciprocal pairs).
Since it's 2×2 with det = 1:
    λ² - tr(A)λ + 1 = 0
    λ = (tr(A) ± √(tr(A)² - 4)) / 2

Returns vector of 2 eigenvalues.
"""
function _sp2_eigenvalues(A)
    _is_sp2(A) || return nothing
    
    t = tr(A)
    discriminant = t^2 - 4
    
    # For symbolic expressions, just use sqrt - it will handle complex results
    # For numeric negative values, use explicit imaginary sqrt
    if discriminant isa AbstractFloat && discriminant < 0
        sqrt_disc = im * sqrt(-discriminant)
    elseif discriminant isa Complex && real(discriminant) < 0 && imag(discriminant) ≈ 0
        sqrt_disc = im * sqrt(-real(discriminant))
    else
        # Symbolic or positive numeric: use standard sqrt
        sqrt_disc = sqrt(discriminant)
    end
    
    λ1 = (t + sqrt_disc) / 2
    λ2 = (t - sqrt_disc) / 2
    
    return [Symbolics.simplify(λ1), Symbolics.simplify(λ2)]
end

# ============================================================================
# Sp(4) - Symplectic 4×4
# ============================================================================

"""
    _is_sp4(A)

Check if A is a 4×4 symplectic matrix (Sp(4)).
Returns true if it is, false otherwise.
"""
function _is_sp4(A)
    size(A) == (4, 4) || return false
    n = _is_symplectic(A)
    return !isnothing(n) && n == 2
end

"""
    _sp4_eigenvalues(A)

Compute eigenvalues of an Sp(4) matrix (4×4 symplectic).

For Sp(4), eigenvalues come in reciprocal pairs: λ₁, 1/λ₁, λ₂, 1/λ₂.
The characteristic polynomial has the form:
    λ⁴ - aλ³ + bλ² - aλ + 1 = 0  (palindromic)

This can be reduced to a quadratic in (λ + 1/λ).

Let μ = λ + 1/λ. Then:
    λ² + 1/λ² = μ² - 2
    
Dividing char poly by λ²: λ² - aλ + b - a/λ + 1/λ² = 0
Rearranging: (λ² + 1/λ²) - a(λ + 1/λ) + b = 0
So: μ² - 2 - aμ + b = 0
    μ² - aμ + (b - 2) = 0

Solve for μ, then λ from λ + 1/λ = μ → λ² - μλ + 1 = 0.
"""
function _sp4_eigenvalues(A)
    _is_sp4(A) || return nothing
    
    # Characteristic polynomial coefficients
    # For symplectic matrices, det = 1, so polynomial is palindromic
    t1 = tr(A)        # coefficient of λ³ (with sign)
    t2 = tr(A * A)
    
    # For the characteristic polynomial λ⁴ - c₃λ³ + c₂λ² - c₁λ + c₀
    # c₃ = tr(A), c₀ = det(A) = 1
    # For symplectic: c₁ = c₃ = tr(A) (palindromic)
    # c₂ = (tr(A)² - tr(A²))/2
    
    c3 = t1
    c2 = (t1^2 - t2) / 2
    # c1 = c3 for symplectic
    # c0 = 1
    
    # Solve μ² - c₃μ + (c₂ - 2) = 0
    a_coeff = c3
    b_coeff = c2 - 2
    
    disc_μ = a_coeff^2 - 4 * b_coeff
    sqrt_disc_μ = sqrt(disc_μ)
    
    μ1 = (a_coeff + sqrt_disc_μ) / 2
    μ2 = (a_coeff - sqrt_disc_μ) / 2
    
    # For each μ, solve λ² - μλ + 1 = 0
    # λ = (μ ± √(μ² - 4)) / 2
    
    disc1 = μ1^2 - 4
    sqrt_disc1 = sqrt(disc1)
    λ1 = (μ1 + sqrt_disc1) / 2
    λ2 = (μ1 - sqrt_disc1) / 2  # = 1/λ1
    
    disc2 = μ2^2 - 4
    sqrt_disc2 = sqrt(disc2)
    λ3 = (μ2 + sqrt_disc2) / 2
    λ4 = (μ2 - sqrt_disc2) / 2  # = 1/λ3
    
    return [Symbolics.simplify(λ1), Symbolics.simplify(λ2),
            Symbolics.simplify(λ3), Symbolics.simplify(λ4)]
end

# ============================================================================
# SO(1,1) - Lorentz Boosts in 1+1D
# ============================================================================

"""
    _is_so11(A)

Check if A is a 2×2 Lorentz boost matrix (SO(1,1)).
These have the form:
[cosh(φ)  sinh(φ)]
[sinh(φ)  cosh(φ)]

Returns the (cosh, sinh) pair if it is, nothing otherwise.
"""
function _is_so11(A)
    size(A) == (2, 2) || return nothing
    _is_indefinite_orthogonal(A, 1, 1) || return nothing
    
    # Check det = 1 for SO (not O)
    d = det(A)
    if d isa Num || d isa Complex{Num}
        _issymzero(d - 1) || return nothing
    elseif d isa Number
        isapprox(d, 1, atol=1e-10) || return nothing
    else
        _issymzero(d - 1) || return nothing
    end
    
    # Structure: [c s; s c] where c² - s² = 1
    c = A[1, 1]
    s = A[1, 2]
    
    # Verify structure
    diff1 = A[2, 1] - s
    diff2 = A[2, 2] - c
    
    for diff in [diff1, diff2]
        if diff isa Num || diff isa Complex{Num}
            if !_issymzero(diff)
                return nothing
            end
        elseif diff isa Number
            if !isapprox(diff, 0, atol=1e-10)
                return nothing
            end
        else
            if !_issymzero(diff)
                return nothing
            end
        end
    end
    
    return (c, s)  # (cosh φ, sinh φ)
end

"""
    _so11_eigenvalues(A)

Compute eigenvalues of an SO(1,1) matrix (Lorentz boost in 1+1D).

For boost with rapidity φ (where cosh(φ), sinh(φ) define the matrix):
eigenvalues = e^(±φ)

Since cosh(φ) = (e^φ + e^(-φ))/2 and sinh(φ) = (e^φ - e^(-φ))/2:
e^φ = cosh(φ) + sinh(φ)
e^(-φ) = cosh(φ) - sinh(φ)
"""
function _so11_eigenvalues(A)
    cs = _is_so11(A)
    isnothing(cs) && return nothing
    c, s = cs
    
    # eigenvalues: cosh(φ) ± sinh(φ) = e^(±φ)
    return [c + s, c - s]
end

# ============================================================================
# SO(1,3) - Lorentz Group in 3+1D
# ============================================================================

"""
    _is_so13(A)

Check if A is a 4×4 Lorentz transformation (SO(1,3)).
Returns true if it is, false otherwise.
"""
function _is_so13(A)
    size(A) == (4, 4) || return false
    sig = _is_lorentz_group(A)
    return !isnothing(sig) && sig == (1, 3)
end

"""
    _so13_eigenvalues(A)

Compute eigenvalues of an SO(1,3) matrix (Lorentz transformation).

A general Lorentz transformation can be decomposed into:
- A boost (hyperbolic rotation) in some direction
- A spatial rotation

The eigenvalue structure depends on this decomposition:
- Pure rotations: eigenvalues 1, 1, e^(±iθ) 
- Pure boosts: eigenvalues e^(±φ), 1, 1 (but with different eigenvector structure)
- General: mix of the above

For the general case, we use the quartic formula since the characteristic
polynomial is degree 4.

Note: This is a complex case. For now, we detect the structure and use
the constraint that eigenvalues satisfy |λ| = 1 or come in reciprocal pairs.
"""
function _so13_eigenvalues(A)
    _is_so13(A) || return nothing
    
    # The characteristic polynomial is quartic
    # For Lorentz transformations, eigenvalues have special structure:
    # - They come in pairs: if λ is eigenvalue, so is 1/λ (or λ̄ for complex)
    # - For proper orthochronous Lorentz transformations: one real eigenvalue ≥ 1
    
    # Use the standard quartic solver via characteristic polynomial
    # The structure detection is already done; just solve the polynomial
    
    # Compute characteristic polynomial coefficients
    I4 = Matrix{eltype(A)}(I, 4, 4)
    @variables λ_temp
    char_matrix = A - λ_temp * I4
    
    # Actually, let's just use the existing machinery
    # by falling back to the general eigenvalue computation
    # but we know the eigenvalues have Lorentz structure
    
    # For now, return nothing to indicate we detected the structure
    # but should use the general solver with the knowledge that
    # eigenvalues come in reciprocal pairs
    return nothing  # Fall back to general quartic solver
end

# ============================================================================
# General O(n), U(n), Sp(2n) - Structure Detection with Hints
# ============================================================================

"""
    _orthogonal_eigenvalue_structure(n)

Return a description of eigenvalue structure for O(n) matrices.
Eigenvalues are on the unit circle, come in conjugate pairs, and may include ±1.
"""
function _orthogonal_eigenvalue_structure(n)
    return (
        group = :O,
        dimension = n,
        eigenvalue_constraint = :unit_circle,
        pairing = :conjugate,
        real_eigenvalues = [:plus_one, :minus_one]
    )
end

"""
    _unitary_eigenvalue_structure(n)

Return a description of eigenvalue structure for U(n) matrices.
All eigenvalues are on the unit circle.
"""
function _unitary_eigenvalue_structure(n)
    return (
        group = :U,
        dimension = n,
        eigenvalue_constraint = :unit_circle,
        determinant_constraint = :unit_circle  # det is e^(iθ) for some θ
    )
end

"""
    _special_unitary_eigenvalue_structure(n)

Return a description of eigenvalue structure for SU(n) matrices.
Eigenvalues are on the unit circle and their product is 1.
"""
function _special_unitary_eigenvalue_structure(n)
    return (
        group = :SU,
        dimension = n,
        eigenvalue_constraint = :unit_circle,
        determinant_constraint = :one,  # Product of eigenvalues = 1
        phase_constraint = :sum_zero    # Sum of phases = 0 (mod 2π)
    )
end

"""
    _symplectic_eigenvalue_structure(n)

Return a description of eigenvalue structure for Sp(2n) matrices.
Eigenvalues come in reciprocal pairs: if λ is eigenvalue, so is 1/λ.
"""
function _symplectic_eigenvalue_structure(n)
    return (
        group = :Sp,
        half_dimension = n,
        full_dimension = 2n,
        eigenvalue_constraint = :reciprocal_pairs,
        determinant_constraint = :one
    )
end

# ============================================================================
# Master Detection and Dispatch Function
# ============================================================================

"""
    _detect_lie_group(A)

Detect if matrix A belongs to a known Lie group with closed-form eigenvalues.

Returns a tuple (group_symbol, params) where:
- group_symbol is one of :SO2, :SO3, :SO4, :SO5, :SO6, :SO7, :SO8, :SO9,
  :SU2, :SU3, :Sp2, :Sp4, :SO11, :SO13, 
  :orthogonal, :unitary, :special_unitary, :symplectic, or nothing
- params contains group-specific parameters

Priority is given to specific small groups with closed-form eigenvalue formulas.
For SO(n) with n ≤ 9, we have closed-form eigenvalues using the quartic formula.
"""
function _detect_lie_group(A)
    n = size(A, 1)
    size(A, 2) == n || return (nothing, nothing)
    
    # Try specific groups with closed-form eigenvalues first
    
    # SO(2) - 2D rotations
    if n == 2
        so2_params = _is_so2(A)
        if !isnothing(so2_params)
            return (:SO2, so2_params)
        end
        
        # SO(1,1) - Lorentz boost (check before Sp2 since it's more specific)
        so11_params = _is_so11(A)
        if !isnothing(so11_params)
            return (:SO11, so11_params)
        end
        
        # SU(2)
        if _is_su2(A)
            return (:SU2, nothing)
        end
        
        # Sp(2) ≅ SL(2)
        if _is_sp2(A)
            return (:Sp2, nothing)
        end
    end
    
    # SO(3) - 3D rotations
    if n == 3
        if _is_so3(A)
            return (:SO3, nothing)
        end
        
        # SU(3)
        if _is_su3(A)
            return (:SU3, nothing)
        end
    end
    
    # SO(4) and Sp(4)
    if n == 4
        if _is_so4(A)
            return (:SO4, nothing)
        end
        
        if _is_sp4(A)
            return (:Sp4, nothing)
        end
        
        # SO(1,3) - Lorentz group
        if _is_so13(A)
            return (:SO13, nothing)
        end
    end
    
    # SO(5) - 5D rotations
    if n == 5 && _is_so5(A)
        return (:SO5, nothing)
    end
    
    # SO(6) - 6D rotations
    if n == 6 && _is_so6(A)
        return (:SO6, nothing)
    end
    
    # SO(7) - 7D rotations
    if n == 7 && _is_so7(A)
        return (:SO7, nothing)
    end
    
    # SO(8) - 8D rotations
    if n == 8 && _is_so8(A)
        return (:SO8, nothing)
    end
    
    # SO(9) - 9D rotations (largest with closed-form via quartic)
    if n == 9 && _is_so9(A)
        return (:SO9, nothing)
    end
    
    # General structure detection for larger matrices
    # These provide structure hints but may not have closed-form eigenvalues
    
    so_n = _is_special_orthogonal(A)
    if !isnothing(so_n)
        return (:special_orthogonal, so_n)
    end
    
    if _is_orthogonal(A)
        return (:orthogonal, n)
    end
    
    su_n = _is_special_unitary(A)
    if !isnothing(su_n)
        return (:special_unitary, su_n)
    end
    
    if _is_unitary(A)
        return (:unitary, n)
    end
    
    sp_n = _is_symplectic(A)
    if !isnothing(sp_n)
        return (:symplectic, sp_n)
    end
    
    return (nothing, nothing)
end

"""
    _lie_group_eigenvalues(A)

Attempt to compute eigenvalues using Lie group structure.

Returns the eigenvalues if A belongs to a Lie group with known closed-form formulas,
or nothing if no such structure is detected or the group is too large for closed-form.

Supports:
- SO(n) for n ≤ 9: Uses trace invariants and polynomial root formulas
- SU(2), SU(3): Via trace-based formulas
- Sp(2), Sp(4): Via reciprocal pair structure
- SO(1,1): Hyperbolic eigenvalues
"""
function _lie_group_eigenvalues(A)
    group, params = _detect_lie_group(A)
    
    isnothing(group) && return nothing
    
    # Dispatch to specific eigenvalue functions
    if group == :SO2
        return _so2_eigenvalues(A)
    elseif group == :SO3
        return _so3_eigenvalues(A)
    elseif group == :SO4
        return _so4_eigenvalues(A)
    elseif group == :SO5
        return _so5_eigenvalues(A)
    elseif group == :SO6
        return _so6_eigenvalues(A)
    elseif group == :SO7
        return _so7_eigenvalues(A)
    elseif group == :SO8
        return _so8_eigenvalues(A)
    elseif group == :SO9
        return _so9_eigenvalues(A)
    elseif group == :SU2
        return _su2_eigenvalues(A)
    elseif group == :SU3
        return _su3_eigenvalues(A)
    elseif group == :Sp2
        return _sp2_eigenvalues(A)
    elseif group == :Sp4
        return _sp4_eigenvalues(A)
    elseif group == :SO11
        return _so11_eigenvalues(A)
    elseif group == :SO13
        # SO(1,3) falls back to general solver
        return nothing
    elseif group == :special_orthogonal
        # For n ≤ 9, we should have already matched above
        # For n ≥ 10, no closed form available
        return nothing
    end
    
    # For general orthogonal/unitary/symplectic, we don't have closed forms
    # but we know the eigenvalue structure (unit circle, reciprocal pairs, etc.)
    return nothing
end
