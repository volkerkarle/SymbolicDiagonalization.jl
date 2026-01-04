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
