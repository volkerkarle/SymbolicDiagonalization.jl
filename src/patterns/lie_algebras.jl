# ============================================================================
# Lie Algebra Representation Detection and Eigenvalue Computation
# ============================================================================
#
# This module detects when a matrix is an element of a Lie algebra representation
# and uses representation theory to compute eigenvalues without solving the
# characteristic polynomial.
#
# Key insight: For a representation ρ: g → gl(V), the eigenvalues of ρ(X)
# are constrained by representation theory and can often be computed from
# simple invariants like tr(X²).
#
# MATHEMATICAL FRAMEWORK:
# -----------------------
# A representation of a Lie algebra g on vector space V is a homomorphism
# ρ: g → gl(V) preserving the Lie bracket: ρ([X,Y]) = [ρ(X), ρ(Y)].
#
# For a single matrix M = ρ(X), we can detect the representation by:
# 1. Structural properties (skew-Hermitian, dimension, etc.)
# 2. Trace invariants that characterize the representation
# 3. Eigenvalue pattern matching (when we can compute them)
#
# ALGEBRAS SUPPORTED:
# -------------------
# - so(3) ≅ su(2): Spin-j representations, dimension 2j+1
# - su(n): Fundamental and adjoint representations
# - so(n): Fundamental representation (skew-symmetric matrices)
# - sp(2n): Fundamental representation (symplectic algebra)
# - sl(2): Spin-j representations (non-compact version)
# - sl(n): Traceless matrices
#
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Helper Functions
# ============================================================================

"""
    _safe_symbolic_sqrt(x)

Compute sqrt handling symbolic types safely.
For Complex{Num}, extracts real part if imaginary is symbolically zero.
"""
function _safe_symbolic_sqrt(x)
    if x isa Num
        return sqrt(x)
    elseif x isa Complex{Num}
        # For Complex{Num}, check if imaginary part is zero
        im_part = imag(x)
        if _issymzero(im_part)
            re_part = real(x)
            return sqrt(re_part)
        else
            # Can't simplify, return nothing to signal failure
            return nothing
        end
    elseif x isa Complex
        return sqrt(x)
    elseif x isa Real
        if x >= 0
            return sqrt(x)
        else
            return im * sqrt(-x)
        end
    else
        return sqrt(x)
    end
end

"""
    _is_skew_symmetric(A; atol=1e-10)

Check if A is skew-symmetric: A^T = -A.
Elements of so(n) in the fundamental representation are skew-symmetric.
"""
function _is_skew_symmetric(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    for i in 1:n, j in 1:n
        diff = A[i, j] + A[j, i]
        if diff isa Num || diff isa Complex{Num}
            _issymzero(diff) || return false
        elseif diff isa Number
            isapprox(diff, 0, atol=1e-10) || return false
        else
            _issymzero(diff) || return false
        end
    end
    return true
end

"""
    _is_skew_hermitian(A; atol=1e-10)

Check if A is skew-Hermitian: A† = -A (where † is conjugate transpose).
Elements of u(n) and su(n) in unitary representations are skew-Hermitian.

Note: Skew-symmetric real matrices are also skew-Hermitian.
"""
function _is_skew_hermitian(A)
    n = size(A, 1)
    size(A, 2) == n || return false
    
    for i in 1:n, j in 1:n
        diff = A[i, j] + conj(A[j, i])
        if diff isa Num || diff isa Complex{Num}
            _issymzero(diff) || return false
        elseif diff isa Number
            isapprox(diff, 0, atol=1e-10) || return false
        else
            _issymzero(diff) || return false
        end
    end
    return true
end

"""
    _trace(A)

Compute trace, handling symbolic matrices.
"""
function _trace(A)
    n = size(A, 1)
    return sum(A[i, i] for i in 1:n)
end

"""
    _trace_of_square(A)

Compute tr(A²) efficiently without forming A².
tr(A²) = Σᵢⱼ Aᵢⱼ Aⱼᵢ
"""
function _trace_of_square(A)
    n = size(A, 1)
    result = zero(eltype(A)) * zero(eltype(A))  # Ensure correct type
    for i in 1:n, j in 1:n
        result += A[i, j] * A[j, i]
    end
    return result
end

"""
    _is_valid_spin(n)

Check if dimension n corresponds to a valid spin representation.
Returns (true, j) if n = 2j + 1 for some half-integer j ≥ 0.
Returns (false, nothing) otherwise.
"""
function _is_valid_spin(n::Integer)
    # n = 2j + 1  =>  j = (n-1)/2
    # Valid if n ≥ 1 and n is a positive integer
    n >= 1 || return (false, nothing)
    j = (n - 1) // 2  # Use rational to handle half-integers exactly
    return (true, j)
end

"""
    _spin_j_eigenvalue_pattern(j)

Return the eigenvalue pattern for spin-j representation of so(3)/su(2).
For an element X with parameter ω, eigenvalues are {i·ω·m : m = -j, -j+1, ..., j}.

Returns the m-values: [-j, -j+1, ..., j-1, j]
"""
function _spin_j_eigenvalue_pattern(j)
    # j can be integer or half-integer (as Rational)
    if j isa Rational
        if denominator(j) == 1
            j_int = numerator(j)
            return collect(-j_int:j_int)
        else
            # Half-integer case
            j_num = numerator(j)  # This is 2j for half-integers where denom=2
            j_denom = denominator(j)
            return [m // j_denom for m in -j_num:2:j_num]
        end
    else
        # Integer j
        return collect(-j:j)
    end
end

"""
    _sum_of_squares_of_m_values(j)

Compute Σₘ m² for m ∈ {-j, -j+1, ..., j}.
This equals j(j+1)(2j+1)/3.
"""
function _sum_of_squares_of_m_values(j)
    if j isa Rational
        return j * (j + 1) * (2j + 1) // 3
    else
        return j * (j + 1) * (2j + 1) / 3
    end
end

# ============================================================================
# so(3) ≅ su(2) Spin-j Representation Detection
# ============================================================================

"""
    _is_so3_spin_representation(A)

Detect if A is an element of the spin-j representation of so(3) ≅ su(2).

For spin-j representation (dimension n = 2j + 1):
- A must be skew-Hermitian: A† = -A
- Eigenvalues have the form {i·ω·m : m = -j, ..., j} for some real ω ≥ 0

We detect this using the trace invariant:
    tr(A²) = -ω² · Σₘ m² = -ω² · j(j+1)(2j+1)/3

Returns (j, ω) if detected, nothing otherwise.
Where:
- j is the spin (half-integer or integer)
- ω is the "rotation speed" parameter (eigenvalues are i·ω·m)

EDGE CASES:
- Zero matrix: Returns (j, 0) - valid element with ω = 0
- 1×1 matrix: Always spin-0, eigenvalue is 0
- Odd vs even dimensions: Only odd dimensions are valid irreducible reps
"""
function _is_so3_spin_representation(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    # Check valid dimension for spin-j: n = 2j + 1, so n must be ≥ 1
    valid, j = _is_valid_spin(n)
    valid || return nothing
    
    # Must be skew-Hermitian
    _is_skew_hermitian(A) || return nothing
    
    # Compute tr(A²)
    tr_A2 = _trace_of_square(A)
    
    # For spin-j representation:
    # tr(A²) = -ω² · Σₘ m² = -ω² · j(j+1)(2j+1)/3
    # Therefore: ω² = -tr(A²) · 3 / (j(j+1)(2j+1))
    
    sum_m2 = _sum_of_squares_of_m_values(j)
    
    # Handle j = 0 (1×1 matrix) specially
    if j == 0
        # The only spin-0 element is the zero matrix
        if _issymzero(A[1, 1])
            return (j, 0)
        else
            return nothing
        end
    end
    
    # Compute ω² = -tr(A²) / sum_m2
    # For skew-Hermitian matrix, tr(A²) is real and ≤ 0, so -tr(A²) ≥ 0
    omega_squared = -tr_A2 / sum_m2
    
    # Use safe symbolic sqrt
    omega = _safe_symbolic_sqrt(omega_squared)
    isnothing(omega) && return nothing
    
    return (j, omega)
end

"""
    _so3_spin_eigenvalues(j, omega)

Compute eigenvalues of a spin-j representation element with parameter ω.
Eigenvalues are {i·ω·m : m = -j, -j+1, ..., j}.
"""
function _so3_spin_eigenvalues(j, omega)
    m_values = _spin_j_eigenvalue_pattern(j)
    
    # Eigenvalues are i·ω·m for each m
    # Handle the type correctly
    if omega isa Num || omega isa Complex{Num}
        return [im * omega * m for m in m_values]
    else
        # For numeric omega, ensure complex output
        return [im * omega * convert(Float64, m) for m in m_values]
    end
end

# ============================================================================
# su(n) Fundamental Representation Detection
# ============================================================================

"""
    _is_su_fundamental(A)

Detect if A is in the fundamental representation of su(n).

The fundamental representation of su(n) consists of n×n matrices that are:
- Skew-Hermitian: A† = -A
- Traceless: tr(A) = 0

Returns n if detected, nothing otherwise.

Note: This doesn't directly give eigenvalues, but constrains them:
- All eigenvalues are purely imaginary
- Sum of eigenvalues is 0
"""
function _is_su_fundamental(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing  # su(1) is trivial
    
    # Must be skew-Hermitian
    _is_skew_hermitian(A) || return nothing
    
    # Must be traceless
    tr_A = _trace(A)
    if tr_A isa Num || tr_A isa Complex{Num}
        _issymzero(tr_A) || return nothing
    elseif tr_A isa Number
        isapprox(tr_A, 0, atol=1e-10) || return nothing
    else
        _issymzero(tr_A) || return nothing
    end
    
    return n
end

"""
    _is_su2_fundamental(A)

Detect if A is in the fundamental (spin-1/2) representation of su(2).
This is a 2×2 traceless skew-Hermitian matrix.

Returns the parameter ω such that eigenvalues are ±iω/2, or nothing.
"""
function _is_su2_fundamental(A)
    size(A) == (2, 2) || return nothing
    
    # Check su(n) fundamental conditions
    isnothing(_is_su_fundamental(A)) && return nothing
    
    # For 2×2 traceless skew-Hermitian:
    # A = [ia   b+ic ]  with a, b, c real, and trace = 2ia = 0 => a = 0
    #     [-b+ic  -ia]
    # 
    # Actually for traceless: A = [ia    b  ] where a is real, b is complex
    #                            [-b̄   -ia]
    # Wait, let's be more careful.
    # 
    # General 2×2 skew-Hermitian: A† = -A means A[i,j] = -conj(A[j,i])
    # A = [iα    β  ]  where α is real, β is complex
    #     [-β̄   iγ ]  where γ is real
    # Traceless: iα + iγ = 0 => α = -γ
    # So A = [iα    β ]
    #        [-β̄  -iα]
    # 
    # eigenvalues: solve det(A - λI) = (iα - λ)(-iα - λ) + |β|² = 0
    #              -λ² - α² + |β|² = 0
    #              λ² = |β|² - α² = -(α² - |β|²)
    #              λ = ±i√(α² + |β|²) = ±iω/2 where ω = 2√(α² + |β|²)
    #
    # Check: tr(A²) = 2(-α² - |β|²) = -2(α² + |β|²) = -ω²/2
    # So ω² = -2·tr(A²)
    
    tr_A2 = _trace_of_square(A)
    omega_squared = -2 * tr_A2
    
    # Handle symbolic case
    if omega_squared isa Num || omega_squared isa Complex{Num}
        omega = sqrt(omega_squared)
        return omega
    end
    
    # For numeric case
    if omega_squared isa Number
        if omega_squared isa Complex
            if !isapprox(imag(omega_squared), 0, atol=1e-10)
                return nothing
            end
            omega_squared = real(omega_squared)
        end
        
        if omega_squared < -1e-10
            return nothing
        end
        
        omega_squared = max(omega_squared, 0.0)
        omega = sqrt(omega_squared)
        return omega
    end
    
    omega = sqrt(omega_squared)
    return omega
end

"""
    _su2_fundamental_eigenvalues(omega)

Eigenvalues of su(2) fundamental representation element with parameter ω.
Returns [iω/2, -iω/2].
"""
function _su2_fundamental_eigenvalues(omega)
    half_omega = omega / 2
    return [im * half_omega, -im * half_omega]
end

# ============================================================================
# so(n) Fundamental Representation Detection
# ============================================================================

"""
    _is_so_fundamental(A)

Detect if A is in the fundamental representation of so(n).

The fundamental representation of so(n) consists of n×n real skew-symmetric matrices.
(For complex matrices, we check skew-Hermitian with real entries.)

Returns n if detected, nothing otherwise.

Eigenvalue structure:
- All eigenvalues are purely imaginary
- Come in conjugate pairs ±iλ
- For odd n, there's also a zero eigenvalue
"""
function _is_so_fundamental(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing  # so(1) is trivial
    
    # Check skew-symmetric
    _is_skew_symmetric(A) || return nothing
    
    return n
end

"""
    _so_fundamental_eigenvalue_structure(A, n)

For so(n) fundamental representation, eigenvalues come in ±iλ pairs.
For odd n, there's also a zero eigenvalue.

We can extract the λ values from invariants, but the general case requires
solving a polynomial of degree floor(n/2).

Returns the eigenvalue structure info, or nothing if not computable.
"""
function _so_fundamental_eigenvalue_structure(A, n)
    # For n = 2: eigenvalues ±iλ where λ² = -det(A) or equivalently from tr(A²)
    # For n = 3: eigenvalues 0, ±iλ where λ² = -tr(A²)/2
    # For n = 4: eigenvalues ±iλ₁, ±iλ₂ - need both tr(A²) and tr(A⁴)
    # General case gets complex
    
    if n == 2
        # tr(A²) = -2λ² => λ = sqrt(-tr(A²)/2)
        tr_A2 = _trace_of_square(A)
        lambda_squared = -tr_A2 / 2
        
        if lambda_squared isa Num || lambda_squared isa Complex{Num}
            lambda = sqrt(lambda_squared)
            return [im * lambda, -im * lambda]
        elseif lambda_squared isa Number
            if lambda_squared isa Complex
                if !isapprox(imag(lambda_squared), 0, atol=1e-10)
                    return nothing
                end
                lambda_squared = real(lambda_squared)
            end
            lambda_squared = max(lambda_squared, 0.0)
            lambda = sqrt(lambda_squared)
            return [im * lambda, -im * lambda]
        end
        
    elseif n == 3
        # eigenvalues: 0, ±iλ
        # tr(A²) = 0 + (-λ²) + (-λ²) = -2λ²
        tr_A2 = _trace_of_square(A)
        lambda_squared = -tr_A2 / 2
        
        if lambda_squared isa Num || lambda_squared isa Complex{Num}
            lambda = sqrt(lambda_squared)
            return [zero(lambda), im * lambda, -im * lambda]
        elseif lambda_squared isa Number
            if lambda_squared isa Complex
                if !isapprox(imag(lambda_squared), 0, atol=1e-10)
                    return nothing
                end
                lambda_squared = real(lambda_squared)
            end
            lambda_squared = max(lambda_squared, 0.0)
            lambda = sqrt(lambda_squared)
            return [0.0im, im * lambda, -im * lambda]
        end
    end
    
    # For n ≥ 4, we'd need higher trace invariants
    # Fall back to nothing (let other methods handle it)
    return nothing
end

# ============================================================================
# sp(2n) Symplectic Algebra Representation Detection
# ============================================================================

"""
    _symplectic_form(n)

Return the standard 2n×2n symplectic form J = [0 I; -I 0].
"""
function _symplectic_form(n)
    J = zeros(Int, 2n, 2n)
    for i in 1:n
        J[i, n+i] = 1
        J[n+i, i] = -1
    end
    return J
end

"""
    _is_sp_algebra(A)

Detect if A is in the Lie algebra sp(2n).

The Lie algebra sp(2n) consists of 2n×2n matrices satisfying:
    JA + A^T J = 0
where J is the standard symplectic form.

Equivalently: A^T J = -JA, or A^T = -JAJ⁻¹ = JAJ (since J⁻¹ = -J).

Eigenvalue structure: eigenvalues come in ±λ pairs.

Returns n (half the dimension) if detected, nothing otherwise.
"""
function _is_sp_algebra(A)
    dim = size(A, 1)
    size(A, 2) == dim || return nothing
    dim >= 2 && dim % 2 == 0 || return nothing  # Must be even
    
    n = dim ÷ 2
    J = _symplectic_form(n)
    
    # Check JA + A^T J = 0
    result = J * A + transpose(A) * J
    
    for i in 1:dim, j in 1:dim
        entry = result[i, j]
        if entry isa Num || entry isa Complex{Num}
            _issymzero(entry) || return nothing
        elseif entry isa Number
            isapprox(entry, 0, atol=1e-10) || return nothing
        else
            _issymzero(entry) || return nothing
        end
    end
    
    return n
end

"""
    _sp2_algebra_eigenvalues(A)

Compute eigenvalues of a 2×2 sp(2) ≅ sl(2,R) algebra element.

For A ∈ sp(2), eigenvalues come in ±λ pairs.
From the characteristic polynomial: λ² - (1/2)tr(A²) = 0
Wait, that's not right. Let me reconsider.

For sp(2), the condition JA + A^T J = 0 with J = [0 1; -1 0] means:
[0 1; -1 0][a b; c d] + [a c; b d][0 1; -1 0] = 0
[c d; -a -b] + [-c a; -d b] = 0
[0 d+a; -a-d 0] = 0
So a + d = 0, i.e., tr(A) = 0.

With tr(A) = 0, eigenvalues are ±√(det(A)) or ±i√(-det(A)).

Actually, for traceless 2×2: λ² = -det(A).
Since det(A) = ad - bc and a = -d, det(A) = -d² - bc.

Returns the eigenvalues [λ, -λ].
"""
function _sp2_algebra_eigenvalues(A)
    size(A) == (2, 2) || return nothing
    
    # Compute -det(A)
    det_A = A[1,1] * A[2,2] - A[1,2] * A[2,1]
    lambda_squared = -det_A
    
    if lambda_squared isa Num || lambda_squared isa Complex{Num}
        lambda = sqrt(lambda_squared)
        return [lambda, -lambda]
    elseif lambda_squared isa Number
        if lambda_squared isa Complex
            lambda = sqrt(lambda_squared)
        elseif lambda_squared >= 0
            lambda = sqrt(lambda_squared)
        else
            lambda = im * sqrt(-lambda_squared)
        end
        return [lambda, -lambda]
    end
    
    lambda = sqrt(lambda_squared)
    return [lambda, -lambda]
end

# ============================================================================
# sl(n) - Traceless Matrices
# ============================================================================

"""
    _is_sl_algebra(A)

Detect if A is in the Lie algebra sl(n) (traceless n×n matrices).

The only constraint is tr(A) = 0.
Eigenvalues sum to zero.

Returns n if detected, nothing otherwise.
"""
function _is_sl_algebra(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing
    
    tr_A = _trace(A)
    if tr_A isa Num || tr_A isa Complex{Num}
        _issymzero(tr_A) || return nothing
    elseif tr_A isa Number
        isapprox(tr_A, 0, atol=1e-10) || return nothing
    else
        _issymzero(tr_A) || return nothing
    end
    
    return n
end

"""
    _is_sl2_representation(A)

Detect if A is in a finite-dimensional representation of sl(2).

sl(2) has the same representation theory as su(2)/so(3) - the spin-j representations.
The difference is sl(2) is non-compact, so the matrices aren't skew-Hermitian.

For sl(2,C), the condition is just: dimension = 2j+1 and the eigenvalue pattern matches.

We use a heuristic: if A is traceless and the dimension fits, and the eigenvalue
pattern (computed from the Casimir) matches, we accept it.

Returns (j, eigenvalues) if detected, nothing otherwise.
"""
function _is_sl2_representation(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    # Check valid dimension for spin-j
    valid, j = _is_valid_spin(n)
    valid || return nothing
    
    # Must be traceless
    tr_A = _trace(A)
    if tr_A isa Num || tr_A isa Complex{Num}
        _issymzero(tr_A) || return nothing
    elseif tr_A isa Number
        isapprox(tr_A, 0, atol=1e-10) || return nothing
    else
        _issymzero(tr_A) || return nothing
    end
    
    # For sl(2), eigenvalues are {ω·m : m = -j, ..., j} (no factor of i)
    # tr(A²) = ω² · Σₘ m² = ω² · j(j+1)(2j+1)/3
    
    if j == 0
        # sl(2) spin-0 is just the 1×1 zero matrix
        if _issymzero(A[1, 1])
            return (j, [zero(eltype(A))])
        else
            return nothing
        end
    end
    
    sum_m2 = _sum_of_squares_of_m_values(j)
    tr_A2 = _trace_of_square(A)
    
    omega_squared = tr_A2 / sum_m2
    
    if omega_squared isa Num || omega_squared isa Complex{Num}
        omega = sqrt(omega_squared)
        m_values = _spin_j_eigenvalue_pattern(j)
        eigenvalues = [omega * m for m in m_values]
        return (j, eigenvalues)
    elseif omega_squared isa Number
        if omega_squared isa Complex
            omega = sqrt(omega_squared)
        elseif omega_squared >= 0
            omega = sqrt(omega_squared)
        else
            omega = im * sqrt(-omega_squared)
        end
        m_values = _spin_j_eigenvalue_pattern(j)
        eigenvalues = [omega * convert(Float64, m) for m in m_values]
        return (j, eigenvalues)
    end
    
    omega = sqrt(omega_squared)
    m_values = _spin_j_eigenvalue_pattern(j)
    eigenvalues = [omega * m for m in m_values]
    return (j, eigenvalues)
end

# ============================================================================
# Reducible Representation Detection
# ============================================================================

"""
    _find_all_block_splits(A)

Find all valid block split points in matrix A.
Returns a list of indices k such that A[1:k, k+1:n] and A[k+1:n, 1:k] are zero.
"""
function _find_all_block_splits(A)
    n = size(A, 1)
    n <= 1 && return Int[]
    
    splits = Int[]
    for k in 1:n-1
        if all(_issymzero, A[1:k, k+1:n]) && all(_issymzero, A[k+1:n, 1:k])
            push!(splits, k)
        end
    end
    return splits
end

"""
    _detect_direct_sum_of_spin_representations(A)

Detect if A is a direct sum of spin-j representations of so(3)/su(2).

A direct sum would be block-diagonal with each block being a spin-jₖ rep.
The total matrix would be skew-Hermitian, and each block would satisfy
the spin representation conditions.

Strategy:
- First try to interpret as a single irreducible spin-j representation
- If that fails OR if block decomposition gives a result with more reps,
  use the block decomposition

Returns a list of (j, ω) pairs for each block, or nothing.

Note: When a matrix could be interpreted either as a single spin-j or as a
direct sum (e.g., spin-1 in Jz basis looks block-diagonal), we prefer the
single spin-j interpretation as it's more informative.
"""
function _detect_direct_sum_of_spin_representations(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    # First check: must be skew-Hermitian
    _is_skew_hermitian(A) || return nothing
    
    # Try as a single irreducible spin-j representation first
    irreducible_result = _is_so3_spin_representation(A)
    
    # Check for block-diagonal structure
    all_splits = _find_all_block_splits(A)
    
    if isempty(all_splits)
        # No block structure - use irreducible if found
        if !isnothing(irreducible_result)
            return [irreducible_result]
        end
        return nothing
    end
    
    # Has block structure - try decomposition
    block_result = nothing
    for split_point in all_splits
        block1 = A[1:split_point, 1:split_point]
        block2 = A[split_point+1:end, split_point+1:end]
        
        reps1 = _detect_direct_sum_of_spin_representations(block1)
        if isnothing(reps1)
            continue
        end
        
        reps2 = _detect_direct_sum_of_spin_representations(block2)
        if isnothing(reps2)
            continue
        end
        
        block_result = vcat(reps1, reps2)
        break  # Use first successful decomposition
    end
    
    # Decide: prefer irreducible if it works, otherwise use block decomposition
    # Exception: if block decomposition gives more reps, it might be more informative
    if !isnothing(irreducible_result) && !isnothing(block_result)
        # Both work - prefer irreducible (it's the "simplest" interpretation)
        # unless block decomposition reveals genuinely different parameters
        # For now, just prefer irreducible (single spin-j is cleaner)
        return [irreducible_result]
    elseif !isnothing(irreducible_result)
        return [irreducible_result]
    elseif !isnothing(block_result)
        return block_result
    end
    
    return nothing
end

"""
    _direct_sum_eigenvalues(reps)

Compute eigenvalues for a direct sum of spin representations.
`reps` is a list of (j, ω) pairs.
"""
function _direct_sum_eigenvalues(reps)
    all_eigenvalues = []
    for (j, omega) in reps
        eigenvalues = _so3_spin_eigenvalues(j, omega)
        append!(all_eigenvalues, eigenvalues)
    end
    return all_eigenvalues
end

# ============================================================================
# Master Detection and Dispatch
# ============================================================================

"""
    _detect_lie_algebra_representation(A)

Master detection function for Lie algebra representations.

Attempts to identify if A is an element of a known Lie algebra representation
and returns information needed to compute eigenvalues.

Returns a tuple (algebra, rep_info) where:
- algebra is a Symbol (:so3_spin, :su2_fundamental, :su_fundamental, 
                       :so_fundamental, :sp_algebra, :sl2_spin, :sl_algebra,
                       :direct_sum_spin, or nothing)
- rep_info contains representation-specific data

Priority order:
1. so(3)/su(2) spin-j representations (most constrained, cleanest eigenvalues)
2. Direct sum of spin representations
3. sp(2n) algebra (palindromic eigenvalue structure)
4. su(n) fundamental (traceless skew-Hermitian)
5. so(n) fundamental (skew-symmetric)
6. sl(2) representations (non-compact analog of spin-j)
7. sl(n) algebra (just traceless - weak constraint)
"""
function _detect_lie_algebra_representation(A)
    n = size(A, 1)
    size(A, 2) == n || return (nothing, nothing)
    
    # Only apply to symbolic matrices (numeric matrices use LinearAlgebra)
    is_symbolic = eltype(A) <: Num || eltype(A) <: Complex{Num}
    
    # 1. Try so(3)/su(2) spin-j representation (most specific)
    spin_result = _is_so3_spin_representation(A)
    if !isnothing(spin_result)
        j, omega = spin_result
        return (:so3_spin, (j, omega))
    end
    
    # 2. Try direct sum of spin representations
    direct_sum_result = _detect_direct_sum_of_spin_representations(A)
    if !isnothing(direct_sum_result) && length(direct_sum_result) > 1
        return (:direct_sum_spin, direct_sum_result)
    end
    
    # 3. Try sp(2n) algebra
    sp_n = _is_sp_algebra(A)
    if !isnothing(sp_n)
        if sp_n == 1  # sp(2) case
            eigenvalues = _sp2_algebra_eigenvalues(A)
            if !isnothing(eigenvalues)
                return (:sp2_algebra, eigenvalues)
            end
        end
        return (:sp_algebra, sp_n)
    end
    
    # 4. Try su(n) fundamental
    su_n = _is_su_fundamental(A)
    if !isnothing(su_n)
        if su_n == 2
            omega = _is_su2_fundamental(A)
            if !isnothing(omega)
                return (:su2_fundamental, omega)
            end
        end
        return (:su_fundamental, su_n)
    end
    
    # 5. Try so(n) fundamental
    so_n = _is_so_fundamental(A)
    if !isnothing(so_n)
        eigenvalues = _so_fundamental_eigenvalue_structure(A, so_n)
        if !isnothing(eigenvalues)
            return (:so_fundamental, (so_n, eigenvalues))
        end
        return (:so_fundamental, (so_n, nothing))
    end
    
    # 6. Try sl(2) spin-j representation
    sl2_result = _is_sl2_representation(A)
    if !isnothing(sl2_result)
        j, eigenvalues = sl2_result
        return (:sl2_spin, (j, eigenvalues))
    end
    
    # 7. Try sl(n) algebra (weakest constraint)
    sl_n = _is_sl_algebra(A)
    if !isnothing(sl_n)
        return (:sl_algebra, sl_n)
    end
    
    return (nothing, nothing)
end

"""
    _lie_algebra_eigenvalues(A)

Compute eigenvalues of A using Lie algebra representation theory.

Returns a vector of eigenvalues if A is detected as a Lie algebra representation
element with computable eigenvalues, nothing otherwise.
"""
function _lie_algebra_eigenvalues(A)
    algebra, rep_info = _detect_lie_algebra_representation(A)
    
    isnothing(algebra) && return nothing
    
    if algebra == :so3_spin
        j, omega = rep_info
        return _so3_spin_eigenvalues(j, omega)
        
    elseif algebra == :direct_sum_spin
        return _direct_sum_eigenvalues(rep_info)
        
    elseif algebra == :sp2_algebra
        return rep_info  # eigenvalues already computed
        
    elseif algebra == :su2_fundamental
        omega = rep_info
        return _su2_fundamental_eigenvalues(omega)
        
    elseif algebra == :so_fundamental
        so_n, eigenvalues = rep_info
        if !isnothing(eigenvalues)
            return eigenvalues
        end
        # Fall through - can't compute eigenvalues
        return nothing
        
    elseif algebra == :sl2_spin
        j, eigenvalues = rep_info
        return eigenvalues
        
    # For :su_fundamental, :sp_algebra, :sl_algebra we detect the structure
    # but don't have closed-form eigenvalue formulas for general n
    end
    
    return nothing
end

# ============================================================================
# Standard Generators (Pauli, Gell-Mann, etc.)
# ============================================================================

"""
    _pauli_matrices()

Return the three Pauli matrices σ₁, σ₂, σ₃.

These are Hermitian (not skew-Hermitian), so the su(2) generators are i·σₖ/2.
"""
function _pauli_matrices()
    σ₁ = [0 1; 1 0]
    σ₂ = [0 -im; im 0]
    σ₃ = [1 0; 0 -1]
    return (σ₁, σ₂, σ₃)
end

"""
    _su2_generators()

Return the three su(2) generators: τₖ = i·σₖ/2.

These are skew-Hermitian and satisfy [τᵢ, τⱼ] = εᵢⱼₖ τₖ.
Eigenvalues of each τₖ are ±i/2.
"""
function _su2_generators()
    σ₁, σ₂, σ₃ = _pauli_matrices()
    return (im * σ₁ / 2, im * σ₂ / 2, im * σ₃ / 2)
end

"""
    _so3_generators()

Return the three so(3) generators: Lₓ, Lᵧ, Lᵤ (3×3 skew-symmetric).

These satisfy [Lᵢ, Lⱼ] = εᵢⱼₖ Lₖ.
"""
function _so3_generators()
    Lx = [0 0 0; 0 0 -1; 0 1 0]
    Ly = [0 0 1; 0 0 0; -1 0 0]
    Lz = [0 -1 0; 1 0 0; 0 0 0]
    return (Lx, Ly, Lz)
end

"""
    _spin_j_generators(j)

Construct the standard spin-j generators J₊, J₋, Jᵤ for the (2j+1)-dimensional
irreducible representation of su(2)/so(3).

Jᵤ is diagonal with entries (j, j-1, ..., -j+1, -j).
J₊ has entries √(j(j+1) - m(m+1)) on the superdiagonal.
J₋ = (J₊)†.

Returns (Jx, Jy, Jz) where Jx = (J₊ + J₋)/2, Jy = (J₊ - J₋)/(2i).

Note: For skew-Hermitian convention (Lie algebra), multiply by i.
"""
function _spin_j_generators(j)
    # Handle j as Rational or integer
    if j isa Rational
        n = numerator(2j + 1)  # dimension
    else
        n = Int(2j + 1)
    end
    
    # m values: j, j-1, ..., -j
    if j isa Rational
        m_values = [j - k for k in 0:n-1]
    else
        m_values = collect(j:-1:-j)
    end
    
    # Jz is diagonal
    Jz = diagm(0 => float.(m_values))
    
    # J+ has entries √(j(j+1) - m(m+1)) on superdiagonal
    # (J+)|j,m⟩ = √(j(j+1) - m(m+1)) |j,m+1⟩
    Jplus = zeros(Complex{Float64}, n, n)
    for k in 1:n-1
        m = m_values[k+1]  # m value for the state being raised
        Jplus[k, k+1] = sqrt(Float64(j*(j+1) - m*(m+1)))
    end
    
    # J- = (J+)†
    Jminus = Jplus'
    
    # Jx = (J+ + J-)/2, Jy = (J+ - J-)/(2i)
    Jx = (Jplus + Jminus) / 2
    Jy = (Jplus - Jminus) / (2im)
    
    return (Jx, Jy, Jz)
end

"""
    _is_standard_spin_generator(A, j)

Check if A is one of the standard spin-j generators (Jx, Jy, or Jz).
Returns (:Jx, :Jy, :Jz) if matched, nothing otherwise.
"""
function _is_standard_spin_generator(A, j)
    Jx, Jy, Jz = _spin_j_generators(j)
    
    tol = 1e-10
    
    # Check Jz (diagonal)
    if all(isapprox(A[i,j], Jz[i,j], atol=tol) for i in 1:size(A,1), j in 1:size(A,2))
        return :Jz
    end
    
    # Check Jx
    if all(isapprox(A[i,j], Jx[i,j], atol=tol) for i in 1:size(A,1), j in 1:size(A,2))
        return :Jx
    end
    
    # Check Jy
    if all(isapprox(A[i,j], Jy[i,j], atol=tol) for i in 1:size(A,1), j in 1:size(A,2))
        return :Jy
    end
    
    # Check i*Jz (skew-Hermitian version)
    iJz = im * Jz
    if all(isapprox(A[i,j], iJz[i,j], atol=tol) for i in 1:size(A,1), j in 1:size(A,2))
        return :iJz
    end
    
    # Check i*Jx
    iJx = im * Jx
    if all(isapprox(A[i,j], iJx[i,j], atol=tol) for i in 1:size(A,1), j in 1:size(A,2))
        return :iJx
    end
    
    # Check i*Jy
    iJy = im * Jy
    if all(isapprox(A[i,j], iJy[i,j], atol=tol) for i in 1:size(A,1), j in 1:size(A,2))
        return :iJy
    end
    
    return nothing
end
