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
# ALGEBRAS SUPPORTED:
# -------------------
# - so(3) ≅ su(2): Spin-j representations, dimension 2j+1
# - su(n): Fundamental representation (traceless skew-Hermitian)
# - so(n): Fundamental representation (skew-symmetric)
# - sp(2n): Symplectic algebra
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
        im_part = imag(x)
        if _issymzero(im_part)
            re_part = real(x)
            return sqrt(re_part)
        else
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
    _trace_of_square(A)

Compute tr(A²) efficiently without forming A².
tr(A²) = Σᵢⱼ Aᵢⱼ Aⱼᵢ
"""
function _trace_of_square(A)
    n = size(A, 1)
    result = zero(eltype(A)) * zero(eltype(A))
    for i in 1:n, j in 1:n
        result += A[i, j] * A[j, i]
    end
    return result
end

"""
    _is_valid_spin(n)

Check if dimension n corresponds to a valid spin representation.
Returns (true, j) if n = 2j + 1 for some half-integer j ≥ 0.
"""
function _is_valid_spin(n::Integer)
    n >= 1 || return (false, nothing)
    j = (n - 1) // 2
    return (true, j)
end

"""
    _spin_j_eigenvalue_pattern(j)

Return the m-values for spin-j representation: [-j, -j+1, ..., j-1, j]
"""
function _spin_j_eigenvalue_pattern(j)
    if j isa Rational
        if denominator(j) == 1
            j_int = numerator(j)
            return collect(-j_int:j_int)
        else
            j_num = numerator(j)
            j_denom = denominator(j)
            return [m // j_denom for m in -j_num:2:j_num]
        end
    else
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
"""
function _is_so3_spin_representation(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    valid, j = _is_valid_spin(n)
    valid || return nothing
    
    _is_skew_hermitian(A) || return nothing
    
    tr_A2 = _trace_of_square(A)
    sum_m2 = _sum_of_squares_of_m_values(j)
    
    # Handle j = 0 (1×1 matrix) specially
    if j == 0
        if _issymzero(A[1, 1])
            return (j, 0)
        else
            return nothing
        end
    end
    
    # ω² = -tr(A²) / sum_m2
    omega_squared = -tr_A2 / sum_m2
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
    
    if omega isa Num || omega isa Complex{Num}
        return [im * omega * m for m in m_values]
    else
        return [im * omega * convert(Float64, m) for m in m_values]
    end
end

# ============================================================================
# su(n) Fundamental Representation Detection
# ============================================================================

"""
    _is_su_fundamental(A)

Detect if A is in the fundamental representation of su(n).
Returns n if detected (traceless skew-Hermitian), nothing otherwise.
"""
function _is_su_fundamental(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing
    
    _is_skew_hermitian(A) || return nothing
    _is_traceless(A) || return nothing
    
    return n
end

"""
    _is_su2_fundamental(A)

Detect if A is in the fundamental (spin-1/2) representation of su(2).
Returns the parameter ω such that eigenvalues are ±iω/2, or nothing.
"""
function _is_su2_fundamental(A)
    size(A) == (2, 2) || return nothing
    isnothing(_is_su_fundamental(A)) && return nothing
    
    tr_A2 = _trace_of_square(A)
    omega_squared = -2 * tr_A2
    
    if omega_squared isa Num || omega_squared isa Complex{Num}
        return sqrt(omega_squared)
    end
    
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
        return sqrt(max(omega_squared, 0.0))
    end
    
    return sqrt(omega_squared)
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
Returns n if detected (skew-symmetric), nothing otherwise.
"""
function _is_so_fundamental(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing
    
    _is_skew_symmetric(A) || return nothing
    
    return n
end

"""
    _so_fundamental_eigenvalue_structure(A, n)

Compute eigenvalues for so(n) fundamental representation.
Eigenvalues come in ±iλ pairs; for odd n, includes zero.
"""
function _so_fundamental_eigenvalue_structure(A, n)
    if n == 2
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
    
    # For n ≥ 4, need higher trace invariants
    return nothing
end

# ============================================================================
# sp(2n) Symplectic Algebra Detection
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
Returns n (half the dimension) if detected, nothing otherwise.

Condition: JA + A^T J = 0 where J is the standard symplectic form.
"""
function _is_sp_algebra(A)
    dim = size(A, 1)
    size(A, 2) == dim || return nothing
    dim >= 2 && dim % 2 == 0 || return nothing
    
    n = dim ÷ 2
    J = _symplectic_form(n)
    
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
Returns [λ, -λ] where λ² = -det(A).
"""
function _sp2_algebra_eigenvalues(A)
    size(A) == (2, 2) || return nothing
    
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
    
    return [sqrt(lambda_squared), -sqrt(lambda_squared)]
end

# ============================================================================
# sl(n) - Traceless Matrices
# ============================================================================

"""
    _is_sl_algebra(A)

Detect if A is in the Lie algebra sl(n) (traceless n×n matrices).
Returns n if detected, nothing otherwise.
"""
function _is_sl_algebra(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    n >= 2 || return nothing
    
    _is_traceless(A) || return nothing
    
    return n
end

"""
    _is_sl2_representation(A)

Detect if A is in a spin-j representation of sl(2).
Returns (j, eigenvalues) if detected, nothing otherwise.
"""
function _is_sl2_representation(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    valid, j = _is_valid_spin(n)
    valid || return nothing
    
    _is_traceless(A) || return nothing
    
    if j == 0
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
    return (j, [omega * m for m in m_values])
end

# ============================================================================
# Reducible Representation Detection
# ============================================================================

"""
    _find_all_block_splits(A)

Find all valid block split points in matrix A.
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

Detect if A is a direct sum of spin-j representations.
Returns a list of (j, ω) pairs for each block, or nothing.
"""
function _detect_direct_sum_of_spin_representations(A)
    n = size(A, 1)
    size(A, 2) == n || return nothing
    
    _is_skew_hermitian(A) || return nothing
    
    # Try as single irreducible first
    irreducible_result = _is_so3_spin_representation(A)
    
    # Check for block structure
    all_splits = _find_all_block_splits(A)
    
    if isempty(all_splits)
        if !isnothing(irreducible_result)
            return [irreducible_result]
        end
        return nothing
    end
    
    # Try decomposition
    block_result = nothing
    for split_point in all_splits
        block1 = A[1:split_point, 1:split_point]
        block2 = A[split_point+1:end, split_point+1:end]
        
        reps1 = _detect_direct_sum_of_spin_representations(block1)
        isnothing(reps1) && continue
        
        reps2 = _detect_direct_sum_of_spin_representations(block2)
        isnothing(reps2) && continue
        
        block_result = vcat(reps1, reps2)
        break
    end
    
    # Prefer irreducible if both work
    if !isnothing(irreducible_result) && !isnothing(block_result)
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
"""
function _direct_sum_eigenvalues(reps)
    all_eigenvalues = []
    for (j, omega) in reps
        append!(all_eigenvalues, _so3_spin_eigenvalues(j, omega))
    end
    return all_eigenvalues
end

# ============================================================================
# Standard Spin-j Generators
# ============================================================================

"""
    spin_j_generators(j) -> (Jx, Jy, Jz)

Construct the standard spin-j generators for the (2j+1)-dimensional
irreducible representation of su(2)/so(3).

Jz is diagonal with entries (j, j-1, ..., -j).
Jx = (J₊ + J₋)/2, Jy = (J₊ - J₋)/(2i).

For skew-Hermitian convention (Lie algebra), multiply by i.
"""
function spin_j_generators(j)
    if j isa Rational
        n = numerator(2j + 1)
    else
        n = Int(2j + 1)
    end
    
    if j isa Rational
        m_values = [j - k for k in 0:n-1]
    else
        m_values = collect(j:-1:-j)
    end
    
    Jz = diagm(0 => float.(m_values))
    
    Jplus = zeros(Complex{Float64}, n, n)
    for k in 1:n-1
        m = m_values[k+1]
        Jplus[k, k+1] = sqrt(Float64(j*(j+1) - m*(m+1)))
    end
    
    Jminus = Jplus'
    Jx = (Jplus + Jminus) / 2
    Jy = (Jplus - Jminus) / (2im)
    
    return (Jx, Jy, Jz)
end

"""
    so3_generators() -> (Lx, Ly, Lz)

Return the three so(3) generators (3×3 skew-symmetric).
These satisfy [Lᵢ, Lⱼ] = εᵢⱼₖ Lₖ.
"""
function so3_generators()
    Lx = [0 0 0; 0 0 -1; 0 1 0]
    Ly = [0 0 1; 0 0 0; -1 0 0]
    Lz = [0 -1 0; 1 0 0; 0 0 0]
    return (Lx, Ly, Lz)
end

"""
    su2_generators() -> (τ1, τ2, τ3)

Return the three su(2) generators: τₖ = i·σₖ/2.
These are skew-Hermitian and satisfy [τᵢ, τⱼ] = εᵢⱼₖ τₖ.
"""
function su2_generators()
    σ₁, σ₂, σ₃ = pauli_x(), pauli_y(), pauli_z()
    return (im * σ₁ / 2, im * σ₂ / 2, im * σ₃ / 2)
end

# ============================================================================
# Master Detection and Dispatch
# ============================================================================

"""
    _detect_lie_algebra_representation(A)

Master detection function for Lie algebra representations.

Returns (algebra, rep_info) where algebra is one of:
:so3_spin, :su2_fundamental, :su_fundamental, :so_fundamental,
:sp_algebra, :sl2_spin, :sl_algebra, :direct_sum_spin, or nothing
"""
function _detect_lie_algebra_representation(A)
    n = size(A, 1)
    size(A, 2) == n || return (nothing, nothing)
    
    # 1. so(3)/su(2) spin-j representation
    spin_result = _is_so3_spin_representation(A)
    if !isnothing(spin_result)
        j, omega = spin_result
        return (:so3_spin, (j, omega))
    end
    
    # 2. Direct sum of spin representations
    direct_sum_result = _detect_direct_sum_of_spin_representations(A)
    if !isnothing(direct_sum_result) && length(direct_sum_result) > 1
        return (:direct_sum_spin, direct_sum_result)
    end
    
    # 3. sp(2n) algebra
    sp_n = _is_sp_algebra(A)
    if !isnothing(sp_n)
        if sp_n == 1
            eigenvalues = _sp2_algebra_eigenvalues(A)
            if !isnothing(eigenvalues)
                return (:sp2_algebra, eigenvalues)
            end
        end
        return (:sp_algebra, sp_n)
    end
    
    # 4. su(n) fundamental
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
    
    # 5. so(n) fundamental
    so_n = _is_so_fundamental(A)
    if !isnothing(so_n)
        eigenvalues = _so_fundamental_eigenvalue_structure(A, so_n)
        if !isnothing(eigenvalues)
            return (:so_fundamental, (so_n, eigenvalues))
        end
        return (:so_fundamental, (so_n, nothing))
    end
    
    # 6. sl(2) spin-j representation
    sl2_result = _is_sl2_representation(A)
    if !isnothing(sl2_result)
        j, eigenvalues = sl2_result
        return (:sl2_spin, (j, eigenvalues))
    end
    
    # 7. sl(n) algebra
    sl_n = _is_sl_algebra(A)
    if !isnothing(sl_n)
        return (:sl_algebra, sl_n)
    end
    
    return (nothing, nothing)
end

"""
    _lie_algebra_eigenvalues(A)

Compute eigenvalues using Lie algebra representation theory.
Returns eigenvalues if computable, nothing otherwise.
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
        return rep_info
        
    elseif algebra == :su2_fundamental
        omega = rep_info
        return _su2_fundamental_eigenvalues(omega)
        
    elseif algebra == :so_fundamental
        so_n, eigenvalues = rep_info
        return eigenvalues
        
    elseif algebra == :sl2_spin
        j, eigenvalues = rep_info
        return eigenvalues
    end
    
    return nothing
end
