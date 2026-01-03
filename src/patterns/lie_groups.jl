# ============================================================================
# Lie Group Pattern Detection and Eigenvalue Computation
# Covers SO(n), SU(n), Sp(2n), O(p,q) and related matrix Lie groups
# ============================================================================

# ============================================================================
# General Detection Functions
# ============================================================================

"""
    _is_orthogonal(A)

Check if matrix A is orthogonal: A^T * A = I.
Returns true if A is orthogonal, false otherwise.
"""
function _is_orthogonal(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    ATA = transpose(A) * A
    for i in 1:n, j in 1:n
        target = i == j ? 1 : 0
        if !_issymzero(ATA[i, j] - target)
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
    
    # Check det = 1
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
    sqrt_disc = sqrt(discriminant)
    
    cos_θ1 = (sum_cos + sqrt_disc) / 2
    cos_θ2 = (sum_cos - sqrt_disc) / 2
    
    # sin θ = √(1 - cos²θ)
    sin_θ1 = sqrt(1 - cos_θ1^2)
    sin_θ2 = sqrt(1 - cos_θ2^2)
    
    # Eigenvalues: e^(±iθ₁), e^(±iθ₂)
    λ1 = cos_θ1 + im * sin_θ1
    λ2 = cos_θ1 - im * sin_θ1
    λ3 = cos_θ2 + im * sin_θ2
    λ4 = cos_θ2 - im * sin_θ2
    
    return [Symbolics.simplify(λ1), Symbolics.simplify(λ2), 
            Symbolics.simplify(λ3), Symbolics.simplify(λ4)]
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
- group_symbol is one of :SO2, :SO3, :SO4, :SU2, :SU3, :Sp2, :Sp4, :SO11, :SO13, 
  :orthogonal, :unitary, :special_unitary, :symplectic, or nothing
- params contains group-specific parameters

Priority is given to specific small groups with closed-form eigenvalue formulas.
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
    end
    
    # For general orthogonal/unitary/symplectic, we don't have closed forms
    # but we know the eigenvalue structure (unit circle, reciprocal pairs, etc.)
    return nothing
end
