# ============================================================================
# SO(4) - Special Orthogonal Group in 4D (4D Rotations)
# ============================================================================
#
# This file contains all SO(4) related functionality:
# - Detection: _is_SO4(A)
# - Eigenvalues: _SO4_eigenvalues(A)
#
# SO(4) is the group of 4×4 rotation matrices with eigenvalues:
#   {e^{±iθ₁}, e^{±iθ₂}} where θ₁, θ₂ are the two rotation angles
#
# SO(4) ≅ (SU(2) × SU(2)) / Z₂, so every SO(4) rotation can be expressed
# as a pair of rotations in two planes.
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Detection
# ============================================================================

"""
    _is_SO4(A)

Check if A is a 4×4 rotation matrix (SO(4)).
Returns true if A ∈ SO(4), false otherwise.
"""
function _is_SO4(A)
    size(A) == (4, 4) || return false
    return !isnothing(_is_special_orthogonal(A))
end

# ============================================================================
# Eigenvalues
# ============================================================================

"""
    _SO4_eigenvalues(A)

Compute eigenvalues of an SO(4) matrix (4D rotation).

SO(4) has eigenvalues {e^{±iθ₁}, e^{±iθ₂}} where the angles can be 
extracted from trace invariants:
- tr(A) = 2(cos(θ₁) + cos(θ₂))
- tr(A²) = 2(cos(2θ₁) + cos(2θ₂)) = 4(cos²(θ₁) + cos²(θ₂)) - 4

For symbolic matrices, this involves solving a quadratic for cos(θⱼ).
"""
function _SO4_eigenvalues(A)
    _is_SO4(A) || return nothing
    
    t1 = tr(A)
    t2 = tr(A * A)
    
    # Let x = cos(θ₁), y = cos(θ₂)
    # t1 = 2(x + y), t2 = 4(x² + y²) - 4
    # So: x + y = t1/2, x² + y² = (t2 + 4)/4
    # xy = [(t1/2)² - (t2+4)/4] / 2 = (t1² - t2 - 4) / 8
    
    sum_cos = t1 / 2
    prod_cos = (t1^2 - t2 - 4) / 8
    
    # x, y are roots of: z² - (sum_cos)*z + prod_cos = 0
    discriminant = sum_cos^2 - 4 * prod_cos
    sqrt_disc = sqrt(discriminant)
    
    cos_theta1 = (sum_cos + sqrt_disc) / 2
    cos_theta2 = (sum_cos - sqrt_disc) / 2
    
    # sin(θ) = sqrt(1 - cos²(θ))
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
# Eigenvectors
# ============================================================================

"""
    _SO4_block_structure(A)

Check if A is a block-diagonal SO(4) matrix of the form:
    [R₁  0 ]
    [0   R₂]
where R₁, R₂ are 2×2 rotation matrices.

Returns (R₁, R₂) if detected, nothing otherwise.
"""
function _SO4_block_structure(A)
    size(A) == (4, 4) || return nothing
    
    # Check if off-diagonal blocks are zero
    for i in 1:2, j in 3:4
        if !_issymzero(A[i, j]) || !_issymzero(A[j, i])
            return nothing
        end
    end
    
    # Extract diagonal blocks
    R1 = A[1:2, 1:2]
    R2 = A[3:4, 3:4]
    
    # Check both are SO(2)
    cs1 = _is_SO2(R1)
    cs2 = _is_SO2(R2)
    
    if !isnothing(cs1) && !isnothing(cs2)
        return (R1, R2, cs1, cs2)
    end
    
    return nothing
end

"""
    _SO4_eigenpairs_general(A)

Compute eigenpairs for a general SO(4) matrix using projection operators.

The algorithm:
1. Compute eigenvalues from trace invariants (reduces to quadratic)
2. Build projection operators: P₁ = A² - 2cos(θ₂)A + I projects onto θ₁ eigenspace
3. Use columns of P as basis for invariant 2-planes
4. Project onto 1D eigenspaces using: v = (A - λ̄I) * P[:,1]

This works because:
- P₁ annihilates the θ₂ eigenspace, so its columns span the θ₁ invariant 2-plane
- (A - λ̄I) maps the 2-plane to the 1D eigenspace for λ

Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing on failure.
"""
function _SO4_eigenpairs_general(A)
    # Step 1: Compute rotation angles from trace invariants
    t1 = tr(A)
    t2 = tr(A * A)
    
    # cos(θ₁) + cos(θ₂) = t1/2
    # cos(θ₁) * cos(θ₂) = (t1² - t2 - 4) / 8
    sum_cos = t1 / 2
    prod_cos = (t1^2 - t2 - 4) / 8
    
    # Solve quadratic: z² - sum_cos*z + prod_cos = 0
    discriminant = sum_cos^2 - 4 * prod_cos
    
    # Handle edge case: discriminant should be non-negative for valid SO(4)
    # For symbolic, we proceed and let sqrt handle it
    sqrt_disc = sqrt(discriminant)
    
    cos_θ1 = (sum_cos + sqrt_disc) / 2
    cos_θ2 = (sum_cos - sqrt_disc) / 2
    
    sin_θ1 = sqrt(1 - cos_θ1^2)
    sin_θ2 = sqrt(1 - cos_θ2^2)
    
    # Step 2: Build projection operators onto invariant 2-planes
    # P₁ = (A - λ₂₊I)(A - λ₂₋I) = A² - 2cos(θ₂)A + I projects onto θ₁ eigenspace
    # P₂ = A² - 2cos(θ₁)A + I projects onto θ₂ eigenspace
    I4 = Matrix{eltype(A)}(I, 4, 4)
    P1 = A^2 - 2*cos_θ2*A + I4  # Range = θ₁ eigenspace
    P2 = A^2 - 2*cos_θ1*A + I4  # Range = θ₂ eigenspace
    
    # Step 3: Get basis vectors from columns of P
    # P has rank 2 for generic SO(4), so first column should be nonzero
    u1 = P1[:, 1]  # In θ₁ invariant 2-plane
    u2 = P2[:, 1]  # In θ₂ invariant 2-plane
    
    # Step 4: Project onto 1D eigenspaces
    # For eigenvalue λ = cos(θ) + i*sin(θ), the conjugate is λ̄ = cos(θ) - i*sin(θ)
    # (A - λ̄I) maps the 2-plane to the 1D eigenspace for λ
    λ1_plus = cos_θ1 + im*sin_θ1   # e^{iθ₁}
    λ1_minus = cos_θ1 - im*sin_θ1  # e^{-iθ₁}
    λ2_plus = cos_θ2 + im*sin_θ2   # e^{iθ₂}
    λ2_minus = cos_θ2 - im*sin_θ2  # e^{-iθ₂}
    
    # v for λ₁₊: (A - λ₁₋I) * u1
    v1_plus = (A - λ1_minus*I4) * u1
    # v for λ₁₋: (A - λ₁₊I) * u1
    v1_minus = (A - λ1_plus*I4) * u1
    # v for λ₂₊: (A - λ₂₋I) * u2
    v2_plus = (A - λ2_minus*I4) * u2
    # v for λ₂₋: (A - λ₂₊I) * u2
    v2_minus = (A - λ2_plus*I4) * u2
    
    # Apply simplification to eigenvalues and eigenvectors
    λ1_plus = simplify_eigenvalue(λ1_plus)
    λ1_minus = simplify_eigenvalue(λ1_minus)
    λ2_plus = simplify_eigenvalue(λ2_plus)
    λ2_minus = simplify_eigenvalue(λ2_minus)
    
    v1_plus = aggressive_simplify.(v1_plus)
    v1_minus = aggressive_simplify.(v1_minus)
    v2_plus = aggressive_simplify.(v2_plus)
    v2_minus = aggressive_simplify.(v2_minus)
    
    return [
        (λ1_plus, [v1_plus]),
        (λ1_minus, [v1_minus]),
        (λ2_plus, [v2_plus]),
        (λ2_minus, [v2_minus])
    ]
end

"""
    _SO4_eigenpairs(A)

Compute eigenvalue-eigenvector pairs for an SO(4) matrix.

For block-diagonal SO(4) of the form [R₁ 0; 0 R₂], returns simple closed-form eigenpairs.
For general SO(4), uses the projection operator algorithm.

Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing.
"""
function _SO4_eigenpairs(A)
    size(A) == (4, 4) || return nothing
    
    # Check for block-diagonal structure first (gives simpler expressions)
    block_info = _SO4_block_structure(A)
    
    if !isnothing(block_info)
        R1, R2, (c1, s1), (c2, s2) = block_info
        
        # For block-diagonal, eigenvectors are padded versions of R₁, R₂ eigenvectors
        # R₁ eigenvalues: c1 ± im*s1, eigenvectors [1, ±im]
        # R₂ eigenvalues: c2 ± im*s2, eigenvectors [1, ±im]
        
        λ1_minus = c1 - im*s1  # e^{-iθ₁}
        λ1_plus = c1 + im*s1   # e^{iθ₁}
        λ2_minus = c2 - im*s2  # e^{-iθ₂}
        λ2_plus = c2 + im*s2   # e^{iθ₂}
        
        v1_minus = [1, im, 0, 0]       # For e^{-iθ₁}
        v1_plus = [1, -im, 0, 0]       # For e^{iθ₁}
        v2_minus = [0, 0, 1, im]       # For e^{-iθ₂}
        v2_plus = [0, 0, 1, -im]       # For e^{iθ₂}
        
        return [
            (λ1_minus, [v1_minus]),
            (λ1_plus, [v1_plus]),
            (λ2_minus, [v2_minus]),
            (λ2_plus, [v2_plus])
        ]
    end
    
    # For general SO(4), use the projection operator algorithm
    return _SO4_eigenpairs_general(A)
end
