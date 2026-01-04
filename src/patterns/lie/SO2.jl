# ============================================================================
# SO(2) - Special Orthogonal Group in 2D (2D Rotations)
# ============================================================================
#
# This file contains all SO(2) related functionality:
# - Constructors: SO2_rotation(θ)
# - Detection: _is_SO2(A)
# - Eigenvalues: _SO2_eigenvalues(A)
# - Kronecker products: SO2_kron, SO2_kron_eigenvalues, detection
#
# SO(2) is the group of 2×2 rotation matrices with form:
#   [cos(θ)  -sin(θ)]
#   [sin(θ)   cos(θ)]
#
# Eigenvalues are e^{±iθ} = cos(θ) ± i·sin(θ)
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Constructor
# ============================================================================

"""
    SO2_rotation(θ) -> Matrix{Num}

Construct a 2×2 rotation matrix (SO(2) element):

    [cos(θ)  -sin(θ)]
    [sin(θ)   cos(θ)]

Eigenvalues are `cos(θ) ± i·sin(θ) = e^{±iθ}`.

# Example
```julia
@variables θ
R = SO2_rotation(θ)
eigvals(R)  # [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
```
"""
function SO2_rotation(θ)
    c, s = cos(θ), sin(θ)
    return [c -s; s c]
end

# ============================================================================
# Detection
# ============================================================================

"""
    _is_SO2(A)

Check if A is a 2×2 rotation matrix (SO(2)).
Returns (cos(θ), sin(θ)) if it is, nothing otherwise.

SO(2) matrices have the form:
    [cos(θ)  -sin(θ)]
    [sin(θ)   cos(θ)]

Uses _is_special_orthogonal with trig-aware simplification.
"""
function _is_SO2(A)
    size(A) == (2, 2) || return nothing
    isnothing(_is_special_orthogonal(A)) && return nothing
    
    c = A[1, 1]  # cos(θ)
    s = A[2, 1]  # sin(θ)
    
    # Verify structure: A[1,2] = -sin(θ), A[2,2] = cos(θ)
    if !_issymzero(A[1, 2] + s) || !_issymzero(A[2, 2] - c)
        return nothing
    end
    
    return (c, s)
end

# ============================================================================
# Eigenvalues
# ============================================================================

"""
    _SO2_eigenvalues(A)

Compute eigenvalues of an SO(2) matrix (2D rotation).

For rotation by angle θ:
    eigenvalues = e^{±iθ} = cos(θ) ± i·sin(θ)

This is the IDEAL case for symbolic diagonalization - eigenvalues
are expressed directly in terms of the matrix elements.
"""
function _SO2_eigenvalues(A)
    cs = _is_SO2(A)
    isnothing(cs) && return nothing
    c, s = cs
    
    # eigenvalues: cos(θ) ± i·sin(θ)
    return [c + im*s, c - im*s]
end

# ============================================================================
# Kronecker Products - Constructors
# ============================================================================

"""
    SO2_kron_eigenvalues(angles::Vector) -> Vector

Compute eigenvalues of R(θ₁) ⊗ R(θ₂) ⊗ ... ⊗ R(θₖ) in clean trigonometric form.

The eigenvalues are `e^{i(±θ₁±θ₂±...±θₖ)}` for all 2^k sign combinations:
    cos(±θ₁±θ₂±...±θₖ) + i·sin(±θ₁±θ₂±...±θₖ)

# Example
```julia
@variables α β γ
vals = SO2_kron_eigenvalues([α, β, γ])
# 8 eigenvalues: cos(±α±β±γ) + i·sin(±α±β±γ)
```
"""
function SO2_kron_eigenvalues(angles::Vector)
    k = length(angles)
    k >= 1 || error("Need at least one angle")
    
    # Generate all 2^k sign combinations
    eigenvalues = []
    for signs in Iterators.product(fill([-1, 1], k)...)
        # Compute the sum ±θ₁ ± θ₂ ± ... ± θₖ
        angle_sum = sum(s * θ for (s, θ) in zip(signs, angles))
        # Clean eigenvalue: cos(sum) + i·sin(sum)
        push!(eigenvalues, cos(angle_sum) + im * sin(angle_sum))
    end
    
    return eigenvalues
end

"""
    SO2_kron(angles::Vector) -> Matrix

Construct the Kronecker product R(θ₁) ⊗ R(θ₂) ⊗ ... ⊗ R(θₖ).

# Example
```julia
@variables α β
K = SO2_kron([α, β])  # 4×4 matrix
eigvals(K)  # [cos(α+β) + i·sin(α+β), ...]
```
"""
function SO2_kron(angles::Vector)
    length(angles) >= 1 || error("Need at least one angle")
    return reduce(kron, [SO2_rotation(θ) for θ in angles])
end

# ============================================================================
# Kronecker Products - Detection
# ============================================================================

"""
    _detect_SO2_kronecker_product(mat)

Detect if mat is a Kronecker product of SO(2) rotation matrices.

For SO(2) ⊗ SO(2) (4×4 matrix), eigenvalues are e^{i(±θ±φ)}.
For SO(2)^⊗k (2^k × 2^k matrix), eigenvalues are e^{i(±θ₁±θ₂±...±θₖ)}.

Uses clean extraction approach that extracts (cos, sin) pairs directly 
from matrix entries, avoiding sqrt(cos²) expressions.

Returns a vector of eigenvalues if detected, nothing otherwise.
"""
function _detect_SO2_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Only for symbolic matrices - numeric matrices go through generic Kronecker detection
    _is_numeric_matrix(mat) && return nothing
    
    # Check if dimension is a power of 2
    N < 2 && return nothing
    log2_N = log2(N)
    isinteger(log2_N) || return nothing
    k = Int(log2_N)
    
    # For 2×2 matrices, use direct SO(2) detection (handled by _SO2_eigenvalues)
    k < 2 && return nothing
    
    # Try to extract cos/sin pairs directly from the first column structure
    pairs = try
        _extract_SO2_kron_pairs(mat, k)
    catch
        nothing
    end
    
    isnothing(pairs) && return nothing
    
    # Compute clean eigenvalues from the pairs
    return _SO2_kron_eigenvalues_from_pairs(pairs)
end

"""
    _extract_SO2_kron_pairs(mat, k) -> Vector{Tuple}

Extract (cos(θⱼ), sin(θⱼ)) pairs from a 2^k × 2^k SO(2) Kronecker product.
Works by recursively factoring the first column structure.

For mat = R₁ ⊗ M:
- col1[1:half] = c₁ * M[:,1]
- col1[half+1:N] = s₁ * M[:,1]

Extracts c₁, s₁ by dividing block entries by M[1,1] = product of inner cosines.
"""
function _extract_SO2_kron_pairs(mat, k)
    N = 2^k
    col1 = mat[:, 1]
    
    if k == 1
        c, s = col1[1], col1[2]
        # Verify this is actually an SO(2) first column: c² + s² = 1
        if !_issymzero_trig(Symbolics.simplify(c^2 + s^2 - 1))
            return nothing
        end
        return [(c, s)]
    end
    
    half = N ÷ 2
    
    block_11 = mat[1:half, 1:half]
    block_12 = mat[1:half, half+1:N]
    block_21 = mat[half+1:N, 1:half]
    block_22 = mat[half+1:N, half+1:N]
    
    # Verify block structure: block_11 = block_22, block_12 = -block_21
    for i in 1:half, j in 1:half
        if !_issymzero(Symbolics.simplify(block_11[i,j] - block_22[i,j]))
            return nothing
        end
        if !_issymzero(Symbolics.simplify(block_12[i,j] + block_21[i,j]))
            return nothing
        end
    end
    
    # Try pattern matching approach first: extract trig functions from squared sums
    c1_squared = Symbolics.simplify(sum(block_11[i, 1]^2 for i in 1:half))
    s1_squared = Symbolics.simplify(sum(block_21[i, 1]^2 for i in 1:half))
    
    trig_result = _try_extract_trig_from_squared(c1_squared, s1_squared)
    
    if !isnothing(trig_result)
        c1, s1 = trig_result
        # Successfully extracted c₁ = cos(θ₁), s₁ = sin(θ₁)
        M = Symbolics.simplify.(block_11 ./ c1)
        inner_pairs = _extract_SO2_kron_pairs(M, k - 1)
        if !isnothing(inner_pairs)
            return vcat([(c1, s1)], inner_pairs)
        end
    end
    
    # Fallback: recursive structure using M[1,1] = product of inner cosines
    inner_pairs = _extract_SO2_kron_pairs(block_11, k - 1)
    if isnothing(inner_pairs)
        return nothing
    end
    
    # M[1,1] = product of first elements of inner pairs (all the "c" values)
    M_11 = prod(p[1] for p in inner_pairs)
    
    # Extract c₁ and s₁
    c1 = Symbolics.simplify(mat[1, 1] / M_11)
    s1 = Symbolics.simplify(mat[half + 1, 1] / M_11)
    
    # Verify c₁² + s₁² = 1
    check = Symbolics.simplify(c1^2 + s1^2)
    if !_issymzero_trig(Symbolics.simplify(check - 1))
        return nothing
    end
    
    return vcat([(c1, s1)], inner_pairs)
end

"""
    _SO2_kron_eigenvalues_from_pairs(pairs::Vector) -> Vector

Compute eigenvalues from (cos, sin) pairs.
If pairs are directly (cos(θ), sin(θ)), extracts θ and uses the angle-sum formula
for cleaner output: cos(±θ₁±θ₂±...) + i·sin(±θ₁±θ₂±...).
"""
function _SO2_kron_eigenvalues_from_pairs(pairs)
    k = length(pairs)
    
    # Try to extract angles directly from the pairs
    angles = Any[]
    for (c, s) in pairs
        θ = _try_extract_angle_from_cos_sin(c, s)
        if isnothing(θ)
            # Fallback to multiplication approach
            return _SO2_kron_eigenvalues_multiply(pairs)
        end
        push!(angles, θ)
    end
    
    # Use the clean angle-sum formula
    eigenvalues = []
    for signs in Iterators.product(fill([-1, 1], k)...)
        angle_sum = sum(s * θ for (s, θ) in zip(signs, angles))
        push!(eigenvalues, cos(angle_sum) + im * sin(angle_sum))
    end
    
    return eigenvalues
end

"""
    _SO2_kron_eigenvalues_multiply(pairs::Vector) -> Vector

Fallback: compute eigenvalues by multiplying (c ± is) factors.
Used when angle extraction fails.
"""
function _SO2_kron_eigenvalues_multiply(pairs)
    k = length(pairs)
    
    eigenvalues = [1 + 0im]
    
    for (c, s) in pairs
        new_eigenvalues = []
        for λ in eigenvalues
            push!(new_eigenvalues, Symbolics.simplify(λ * (c + im*s)))
            push!(new_eigenvalues, Symbolics.simplify(λ * (c - im*s)))
        end
        eigenvalues = new_eigenvalues
    end
    
    return [trig_simplify(λ) for λ in eigenvalues]
end

"""
    _try_extract_angle_from_cos_sin(c, s) -> Union{Num, Nothing}

Try to extract θ from (cos(θ), sin(θ)) pair by pattern matching.
Returns θ if both c and s are trig functions of the same angle, nothing otherwise.
"""
function _try_extract_angle_from_cos_sin(c, s)
    unwrapped_c = Symbolics.unwrap(c)
    unwrapped_s = Symbolics.unwrap(s)
    
    # Check if c = cos(θ) for some θ
    if !Symbolics.iscall(unwrapped_c) || Symbolics.operation(unwrapped_c) !== cos
        return nothing
    end
    θ_c = Num(Symbolics.arguments(unwrapped_c)[1])
    
    # Check if s = sin(θ) for the same θ
    if !Symbolics.iscall(unwrapped_s) || Symbolics.operation(unwrapped_s) !== sin
        return nothing
    end
    θ_s = Num(Symbolics.arguments(unwrapped_s)[1])
    
    # Verify same angle
    if !_issymzero(Symbolics.simplify(θ_c - θ_s))
        return nothing
    end
    
    return θ_c
end

"""
    _try_extract_trig_from_squared(c_squared, s_squared)

Try to recognize c_squared and s_squared as cos²(θ) and sin²(θ) for some θ.
Returns (cos(θ), sin(θ)) if successful, nothing otherwise.

This avoids the sqrt(cos²(θ)) problem by pattern matching.
"""
function _try_extract_trig_from_squared(c_squared, s_squared)
    unwrapped_c = Symbolics.unwrap(c_squared)
    unwrapped_s = Symbolics.unwrap(s_squared)
    
    # Check if c_squared = f^2 for some f
    try
        if Symbolics.iscall(unwrapped_c) && Symbolics.operation(unwrapped_c) === (^)
            args = Symbolics.arguments(unwrapped_c)
            if length(args) == 2 && isequal(Symbolics.value(args[2]), 2)
                base = args[1]
                # Check if base is cos(something)
                if Symbolics.iscall(base)
                    op = Symbolics.operation(base)
                    if op === cos
                        θ = Symbolics.arguments(base)[1]
                        c = Num(base)  # cos(θ)
                        s = sin(Num(θ))  # sin(θ)
                        
                        # Verify s_squared = sin²(θ)
                        expected_s_squared = Symbolics.simplify(s^2)
                        if _issymzero(Symbolics.simplify(s_squared - expected_s_squared))
                            return (c, s)
                        end
                    end
                end
            end
        end
    catch
        # Pattern matching failed
    end
    
    return nothing
end

# ============================================================================
# SO(2) Kronecker Decomposition (for recursive factorization)
# ============================================================================

"""
    _try_SO2_kronecker_decomposition(mat, k)

Try to decompose a 2^k × 2^k orthogonal matrix as R₁ ⊗ R₂ ⊗ ... ⊗ Rₖ
where each Rᵢ is a 2×2 rotation matrix.

Returns eigenvalues if successful, nothing otherwise.
"""
function _try_SO2_kronecker_decomposition(mat, k)
    result = _try_SO2_kronecker_decomposition_impl(mat, k)
    if isnothing(result)
        return nothing
    end
    
    eigenvalues, _ = result
    return eigenvalues
end

"""
    _try_SO2_kronecker_decomposition_impl(mat, k)

Implementation that returns (eigenvalues, angle_info) where angle_info
contains the (cos, sin) pairs for each rotation factor.
"""
function _try_SO2_kronecker_decomposition_impl(mat, k)
    N = size(mat, 1)
    
    if k == 1
        # Base case: 2×2 matrix, should be SO(2)
        cs = _is_SO2(mat)
        if !isnothing(cs)
            c, s = cs
            eigenvalues = [c + im*s, c - im*s]
            return (eigenvalues, [(c, s)])
        end
        return nothing
    end
    
    # For k > 1: Try to factor as R ⊗ M where R is 2×2 and M is 2^(k-1) × 2^(k-1)
    half_N = N ÷ 2
    
    block_11 = mat[1:half_N, 1:half_N]
    block_12 = mat[1:half_N, half_N+1:N]
    block_21 = mat[half_N+1:N, 1:half_N]
    block_22 = mat[half_N+1:N, half_N+1:N]
    
    # Check: block_11 should equal block_22 (both are c*M)
    for i in 1:half_N, j in 1:half_N
        diff = Symbolics.simplify(block_11[i, j] - block_22[i, j])
        if !_issymzero(diff)
            return nothing
        end
    end
    
    # Check: block_12 should equal -block_21
    for i in 1:half_N, j in 1:half_N
        diff = Symbolics.simplify(block_12[i, j] + block_21[i, j])
        if !_issymzero(diff)
            return nothing
        end
    end
    
    # Extract c² and s² from block norms
    c_squared = Symbolics.simplify(sum(block_11[i, 1]^2 for i in 1:half_N))
    s_squared = Symbolics.simplify(sum(block_21[i, 1]^2 for i in 1:half_N))
    
    # Verify c² + s² = 1
    sum_check = Symbolics.simplify(c_squared + s_squared - 1)
    if !_issymzero_trig(sum_check)
        return nothing
    end
    
    # Special case for k=2: direct extraction from first column
    if k == 2
        a = mat[1, 1]
        b = mat[2, 1]
        c_entry = mat[3, 1]
        d = mat[4, 1]
        
        # Verify cross-ratio: a*d = b*c
        cross_check = Symbolics.simplify(a*d - b*c_entry)
        if !_issymzero(cross_check)
            return nothing
        end
        
        # Build eigenvalues using angle addition formulas
        raw_eigenvalues = [
            (a - d) + im*(b + c_entry),     # e^{i(θ+φ)}
            (a + d) + im*(c_entry - b),     # e^{i(θ-φ)}
            (a + d) + im*(b - c_entry),     # e^{i(φ-θ)}
            (a - d) + im*(-b - c_entry)     # e^{-i(θ+φ)}
        ]
        
        eigenvalues = [trig_simplify(λ) for λ in raw_eigenvalues]
        return (eigenvalues, nothing)
    end
    
    # For k > 2: Try pattern matching first
    trig_result = _try_extract_trig_from_squared(c_squared, s_squared)
    
    if !isnothing(trig_result)
        c, s = trig_result
        M = Symbolics.simplify.(block_11 ./ c)
        result_M = _try_SO2_kronecker_decomposition_impl(M, k - 1)
        if !isnothing(result_M)
            eig_R = [c + im*s, c - im*s]
            eig_M, angle_info_M = result_M
            eigenvalues = [Symbolics.simplify(λ * μ) for λ in eig_R for μ in eig_M]
            angle_info = isnothing(angle_info_M) ? nothing : vcat([(c, s)], angle_info_M)
            return (eigenvalues, angle_info)
        end
    end
    
    # Fallback: use sqrt (less clean)
    c = Symbolics.simplify(sqrt(c_squared))
    s = Symbolics.simplify(sqrt(s_squared))
    eig_R = [c + im*s, c - im*s]
    
    if _issymzero(c)
        M = block_21
    else
        M = Symbolics.simplify.(block_11 ./ c)
    end
    
    result_M = _try_SO2_kronecker_decomposition_impl(M, k - 1)
    if isnothing(result_M)
        return nothing
    end
    
    eig_M, angle_info_M = result_M
    eigenvalues = [Symbolics.simplify(λ * μ) for λ in eig_R for μ in eig_M]
    angle_info = isnothing(angle_info_M) ? nothing : vcat([(c, s)], angle_info_M)
    return (eigenvalues, angle_info)
end
