# ============================================================================
# SU(2) - Special Unitary Group in 2D
# ============================================================================
#
# This file contains all SU(2) related functionality:
# - Pauli matrices: pauli_x, pauli_y, pauli_z
# - Constructors: SU2_Ux(θ), SU2_Uy(θ), SU2_Uz(θ)
# - Detection: _is_SU2(A)
# - Eigenvalues: _SU2_eigenvalues(A)
# - Kronecker products: SU2_kron, SU2_kron_eigenvalues, detection
#
# SU(2) is the group of 2×2 special unitary matrices.
# Eigenvalues are e^{±iθ/2} = cos(θ/2) ± i·sin(θ/2)
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Pauli Matrices
# ============================================================================

"""
    pauli_x() -> Matrix{Int}

Pauli X matrix (σₓ):
    [0  1]
    [1  0]
"""
pauli_x() = [0 1; 1 0]

"""
    pauli_y() -> Matrix{Complex{Int}}

Pauli Y matrix (σᵧ):
    [0  -i]
    [i   0]
"""
pauli_y() = [0 -im; im 0]

"""
    pauli_z() -> Matrix{Int}

Pauli Z matrix (σᵤ):
    [1   0]
    [0  -1]
"""
pauli_z() = [1 0; 0 -1]

# ============================================================================
# Constructors
# ============================================================================

"""
    SU2_Ux(θ) -> Matrix{Num}

Construct SU(2) rotation matrix around x-axis (spin-1/2):

    U_x(θ) = exp(-i θ σx/2) = cos(θ/2)I - i sin(θ/2)σx

    [cos(θ/2)      -i·sin(θ/2)]
    [-i·sin(θ/2)    cos(θ/2)  ]

Eigenvalues are `e^{±iθ/2}`.
"""
function SU2_Ux(θ)
    c, s = cos(θ/2), sin(θ/2)
    return [c -im*s; -im*s c]
end

"""
    SU2_Uy(θ) -> Matrix{Num}

Construct SU(2) rotation matrix around y-axis (spin-1/2):

    U_y(θ) = exp(-i θ σy/2) = cos(θ/2)I - i sin(θ/2)σy

    [cos(θ/2)   -sin(θ/2)]
    [sin(θ/2)    cos(θ/2)]

Eigenvalues are `e^{±iθ/2}`.

Note: This has the same form as SO(2), but with half-angles.
"""
function SU2_Uy(θ)
    c, s = cos(θ/2), sin(θ/2)
    return [c -s; s c]
end

"""
    SU2_Uz(θ) -> Matrix{Num}

Construct SU(2) rotation matrix around z-axis (spin-1/2):

    U_z(θ) = exp(-i θ σz/2) = cos(θ/2)I - i sin(θ/2)σz

    [e^{-iθ/2}     0    ]
    [   0       e^{iθ/2}]

Equivalently:
    [cos(θ/2) - i·sin(θ/2)          0            ]
    [       0                cos(θ/2) + i·sin(θ/2)]

Eigenvalues are `e^{±iθ/2}`.
"""
function SU2_Uz(θ)
    c, s = cos(θ/2), sin(θ/2)
    return [c - im*s 0; 0 c + im*s]
end

# ============================================================================
# Detection
# ============================================================================

"""
    _is_SU2(A)

Check if A is a 2×2 special unitary matrix (SU(2)).
Returns true if A ∈ SU(2), false otherwise.
"""
function _is_SU2(A)
    size(A) == (2, 2) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _is_SU2_trig(A)

Check if A is an SU(2) matrix using trig-aware simplification.
"""
function _is_SU2_trig(A)
    size(A) == (2, 2) || return false
    
    # Check unitarity: A * A^H = I with trig simplification
    AAH = A * adjoint(A)
    for i in 1:2, j in 1:2
        target = i == j ? 1 : 0
        diff = AAH[i, j] - target
        if !_issymzero_trig(real(diff)) || !_issymzero_trig(imag(diff))
            return false
        end
    end
    
    # Check det = 1 with trig simplification
    d = det(A)
    if !_issymzero_trig(real(d) - 1) || !_issymzero_trig(imag(d))
        return false
    end
    
    return true
end

# ============================================================================
# Eigenvalues
# ============================================================================

"""
    _SU2_eigenvalues(A)

Compute eigenvalues of an SU(2) matrix.

For SU(2), eigenvalues are e^{±iθ/2} where cos(θ/2) = Re(tr(A))/2.
"""
function _SU2_eigenvalues(A)
    _is_SU2(A) || return nothing
    
    tr_A = tr(A)
    cos_theta = real(tr_A) / 2
    sin_theta = aggressive_simplify(sqrt(1 - cos_theta^2))
    
    return [simplify_eigenvalue(cos_theta + im * sin_theta), 
            simplify_eigenvalue(cos_theta - im * sin_theta)]
end

"""
    _SU2_trace(A)

Extract the trace of an SU(2) matrix, returning (cos(θ/2), sin(θ/2)) 
where the eigenvalues are e^{±iθ/2}.

For SU(2), tr(U) = 2·cos(θ/2) where θ is the rotation angle.
"""
function _SU2_trace(A)
    tr_A = tr(A)
    cos_half = real(tr_A) / 2
    sin_half = aggressive_simplify(sqrt(1 - cos_half^2))
    return (cos_half, sin_half)
end

# ============================================================================
# Eigenvectors
# ============================================================================

"""
    _SU2_rotation_type(A)

Detect if A is a rotation around a principal axis (x, y, or z).

Returns:
- :Uz if diagonal (rotation around z-axis): [e^{-iθ/2} 0; 0 e^{iθ/2}]
- :Uy if real (rotation around y-axis): [c -s; s c] (like SO(2) but half-angle)
- :Ux if anti-diagonal imaginary: [c -is; -is c]
- nothing otherwise (general SU(2) rotation)

For axis-aligned rotations, we have closed-form eigenvectors.
"""
function _SU2_rotation_type(A)
    size(A) == (2, 2) || return nothing
    
    # Check for Uz: diagonal matrix [e^{-iθ/2} 0; 0 e^{iθ/2}]
    if _issymzero(A[1, 2]) && _issymzero(A[2, 1])
        return :Uz
    end
    
    # Check for Uy: real matrix [c -s; s c] 
    # (same structure as SO(2), but with half-angle)
    if _issymzero(imag(A[1, 1])) && _issymzero(imag(A[1, 2])) &&
       _issymzero(imag(A[2, 1])) && _issymzero(imag(A[2, 2]))
        # Check SO(2) structure: A[1,1] = A[2,2], A[1,2] = -A[2,1]
        if _issymzero(A[1, 1] - A[2, 2]) && _issymzero(A[1, 2] + A[2, 1])
            return :Uy
        end
    end
    
    # Check for Ux: [c -is; -is c] structure
    # A[1,1] = A[2,2] = c (real)
    # A[1,2] = A[2,1] = -is (purely imaginary)
    if _issymzero(A[1, 1] - A[2, 2]) && _issymzero(A[1, 2] - A[2, 1])
        if _issymzero(imag(A[1, 1])) && _issymzero(real(A[1, 2]))
            return :Ux
        end
    end
    
    return nothing
end

"""
    _SU2_axis_eigenvectors(axis_type)

Return the eigenvectors for an axis-aligned SU(2) rotation.

For eigenvalues {e^{iθ/2}, e^{-iθ/2}}:

Uz (diagonal):
  - λ=e^{-iθ/2}: [1, 0] (standard basis)
  - λ=e^{iθ/2}: [0, 1] (standard basis)

Ux (rotation around x):
  - λ=e^{iθ/2}: [1, 1]/√2
  - λ=e^{-iθ/2}: [1, -1]/√2

Uy (rotation around y):
  - λ=e^{iθ/2}: [1, i]/√2
  - λ=e^{-iθ/2}: [1, -i]/√2

Note: Returns unnormalized eigenvectors.
"""
function _SU2_axis_eigenvectors(axis_type)
    if axis_type == :Uz
        # Diagonal case: standard basis
        return [[1, 0], [0, 1]]
    elseif axis_type == :Ux
        # Ux eigenvectors
        return [[1, 1], [1, -1]]
    elseif axis_type == :Uy
        # Uy eigenvectors (same as SO(2))
        return [[1, im], [1, -im]]
    end
    return nothing
end

"""
    _SU2_eigenpairs(A)

Compute eigenvalue-eigenvector pairs for an SU(2) matrix.

For axis-aligned rotations (Ux, Uy, Uz), returns closed-form eigenpairs.
For general rotations, returns nothing (requires nullspace computation).

Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing.

# Example
```julia
@variables θ
U = SU2_Uz(θ)
pairs = _SU2_eigenpairs(U)
# For Uz diagonal: eigenvectors are [1,0] and [0,1]
```
"""
function _SU2_eigenpairs(A)
    size(A) == (2, 2) || return nothing
    
    # First check axis type - this is fast and doesn't require full SU(2) verification
    axis_type = _SU2_rotation_type(A)
    
    # For general SU(2) rotations, we need nullspace computation
    isnothing(axis_type) && return nothing
    
    # Verify it's actually SU(2) using trig-aware check
    # (The standard _is_SU2 may fail on trig expressions)
    _is_SU2(A) || _is_SU2_trig(A) || return nothing
    
    # Get eigenvectors for this axis type
    vecs = _SU2_axis_eigenvectors(axis_type)
    isnothing(vecs) && return nothing
    
    # For Uz (diagonal), eigenvalues are directly on the diagonal
    if axis_type == :Uz
        λ1 = A[1, 1]  # e^{-iθ/2}
        λ2 = A[2, 2]  # e^{iθ/2}
        return [(λ1, [vecs[1]]), (λ2, [vecs[2]])]
    end
    
    # For Ux and Uy, compute eigenvalues from trace
    # Use trace-based formula which is more robust
    tr_A = tr(A)
    cos_theta = real(tr_A) / 2
    sin_theta = aggressive_simplify(sqrt(1 - cos_theta^2))
    
    λ1 = simplify_eigenvalue(cos_theta + im * sin_theta)
    λ2 = simplify_eigenvalue(cos_theta - im * sin_theta)
    
    # Return as vector of (eigenvalue, [eigenvector]) tuples
    return [(λ1, [vecs[1]]), (λ2, [vecs[2]])]
end

# ============================================================================
# Kronecker Products - Constructors
# ============================================================================

"""
    SU2_kron_eigenvalues(angles::Vector) -> Vector

Compute eigenvalues of U(θ₁) ⊗ U(θ₂) ⊗ ... ⊗ U(θₖ) in clean form.

For SU(2) matrices, eigenvalues are `e^{±iθⱼ/2}` (half-angle compared to SO(3)).
The Kronecker product eigenvalues are products of individual eigenvalues:
    `e^{i(±θ₁±θ₂±...±θₖ)/2}` for all 2^k sign combinations.

This gives:
    cos((±θ₁±θ₂±...±θₖ)/2) + i·sin((±θ₁±θ₂±...±θₖ)/2)
"""
function SU2_kron_eigenvalues(angles::Vector)
    k = length(angles)
    k >= 1 || error("Need at least one angle")
    
    eigenvalues = []
    for signs in Iterators.product(fill([-1, 1], k)...)
        angle_sum = sum(s * θ for (s, θ) in zip(signs, angles))
        # SU(2) eigenvalues use half-angles
        half_angle = angle_sum / 2
        push!(eigenvalues, cos(half_angle) + im * sin(half_angle))
    end
    
    return eigenvalues
end

"""
    SU2_kron(angles::Vector; axis=:z) -> Matrix

Construct the Kronecker product U_axis(θ₁) ⊗ U_axis(θ₂) ⊗ ... ⊗ U_axis(θₖ).

The `axis` parameter can be `:x`, `:y`, or `:z` (default).
"""
function SU2_kron(angles::Vector; axis=:z)
    length(angles) >= 1 || error("Need at least one angle")
    U = axis == :x ? SU2_Ux : axis == :y ? SU2_Uy : SU2_Uz
    return reduce(kron, [U(θ) for θ in angles])
end

# ============================================================================
# Kronecker Products - Detection
# ============================================================================

"""
    _detect_SU2_kronecker_product(mat)

Detect if a matrix is a Kronecker product of SU(2) matrices and compute eigenvalues.

Returns eigenvalues if detected as SU(2)^⊗k, nothing otherwise.
"""
function _detect_SU2_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Must be power of 2
    N >= 4 || return nothing  # Need at least 4×4 for SU(2)⊗SU(2)
    k = Int(round(log2(N)))
    2^k == N || return nothing
    
    # SU(2)⊗SU(2) must have complex entries
    if !_has_complex_entries(mat)
        return nothing
    end
    
    if k == 2  # 4×4 = SU(2) ⊗ SU(2)
        return _detect_SU2_kron_2fold(mat)
    elseif k == 3  # 8×8 = SU(2) ⊗ SU(2) ⊗ SU(2)
        return _detect_SU2_kron_3fold(mat)
    end
    
    return nothing
end

"""
    _detect_SU2_kron_2fold(mat)

Detect if 4×4 matrix is SU(2) ⊗ SU(2) and compute eigenvalues.
"""
function _detect_SU2_kron_2fold(mat)
    size(mat) == (4, 4) || return nothing
    
    B11 = mat[1:2, 1:2]
    B12 = mat[1:2, 3:4]
    B21 = mat[3:4, 1:2]
    B22 = mat[3:4, 3:4]
    
    # Strategy A: Check if any diagonal block is directly SU(2)
    for Bkk in [B11, B22]
        Bkk_simplified = Symbolics.simplify.(Bkk)
        if _is_SU2_trig(Bkk_simplified)
            U2 = Bkk_simplified
            c2, s2 = _SU2_trace(U2)
            
            # Find a non-zero element in U2 to extract U1 entries
            ref_i, ref_j = 1, 1
            if _issymzero_trig(real(U2[1, 1])) && _issymzero_trig(imag(U2[1, 1]))
                ref_i, ref_j = 1, 2
            end
            
            U2_ref = U2[ref_i, ref_j]
            
            U1_11 = Symbolics.simplify(B11[ref_i, ref_j] / U2_ref)
            U1_12 = Symbolics.simplify(B12[ref_i, ref_j] / U2_ref)
            U1_21 = Symbolics.simplify(B21[ref_i, ref_j] / U2_ref)
            U1_22 = Symbolics.simplify(B22[ref_i, ref_j] / U2_ref)
            
            U1 = [U1_11 U1_12; U1_21 U1_22]
            U1_simplified = Symbolics.simplify.(U1)
            
            if _is_SU2_trig(U1_simplified)
                c1, s1 = _SU2_trace(U1_simplified)
                return _compute_SU2_kron_eigenvalues_from_traces(c1, s1, c2, s2)
            end
        end
    end
    
    # Strategy B: Use determinant
    det_B11 = det(B11)
    det_B11_simplified = trig_simplify(Symbolics.simplify(det_B11))
    
    sqrt_det = _try_symbolic_sqrt(det_B11_simplified)
    if !isnothing(sqrt_det)
        U1_11 = sqrt_det
        
        if !_issymzero_trig(real(U1_11)) || !_issymzero_trig(imag(U1_11))
            U2 = Symbolics.simplify.(B11 / U1_11)
            
            if _is_SU2_trig(U2)
                c2, s2 = _SU2_trace(U2)
                
                ref_i, ref_j = 1, 1
                if _issymzero_trig(real(U2[1, 1])) && _issymzero_trig(imag(U2[1, 1]))
                    ref_i, ref_j = 1, 2
                end
                U2_ref = U2[ref_i, ref_j]
                
                U1_12 = Symbolics.simplify(B12[ref_i, ref_j] / U2_ref)
                U1_21 = Symbolics.simplify(B21[ref_i, ref_j] / U2_ref)
                U1_22 = Symbolics.simplify(B22[ref_i, ref_j] / U2_ref)
                
                U1 = [U1_11 U1_12; U1_21 U1_22]
                
                if _is_SU2_trig(U1)
                    c1, s1 = _SU2_trace(U1)
                    return _compute_SU2_kron_eigenvalues_from_traces(c1, s1, c2, s2)
                end
            end
        end
    end
    
    # Strategy C: Orthogonality-based trace extraction
    result = _try_SU2_kron_from_traces(mat, [[B11, B12], [B21, B22]])
    if !isnothing(result)
        return result
    end
    
    return nothing
end

"""
    _detect_SU2_kron_3fold(mat)

Detect if 8×8 matrix is SU(2) ⊗ SU(2) ⊗ SU(2) and compute eigenvalues.
"""
function _detect_SU2_kron_3fold(mat)
    size(mat) == (8, 8) || return nothing
    
    # Try to factor as 4 ⊗ 2
    kron_info = _try_kronecker_factorization(mat, 4, 2)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        B_simplified = Symbolics.simplify.(B)
        if _is_SU2_trig(B_simplified)
            c_B, s_B = _SU2_trace(B_simplified)
            
            vals_A = _detect_SU2_kron_2fold(A)
            
            if !isnothing(vals_A)
                eigenvalues = []
                for λ_A in vals_A
                    λ_B_plus = c_B + im*s_B
                    λ_B_minus = c_B - im*s_B
                    push!(eigenvalues, simplify_eigenvalue(trig_simplify(Symbolics.simplify(λ_A * λ_B_plus))))
                    push!(eigenvalues, simplify_eigenvalue(trig_simplify(Symbolics.simplify(λ_A * λ_B_minus))))
                end
                return eigenvalues
            end
        end
    end
    
    # Try to factor as 2 ⊗ 4
    kron_info = _try_kronecker_factorization(mat, 2, 4)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        A_simplified = Symbolics.simplify.(A)
        if _is_SU2_trig(A_simplified)
            c_A, s_A = _SU2_trace(A_simplified)
            
            vals_B = _detect_SU2_kron_2fold(B)
            
            if !isnothing(vals_B)
                eigenvalues = []
                λ_A_plus = c_A + im*s_A
                λ_A_minus = c_A - im*s_A
                for λ_B in vals_B
                    push!(eigenvalues, simplify_eigenvalue(trig_simplify(Symbolics.simplify(λ_A_plus * λ_B))))
                    push!(eigenvalues, simplify_eigenvalue(trig_simplify(Symbolics.simplify(λ_A_minus * λ_B))))
                end
                return eigenvalues
            end
        end
    end
    
    return nothing
end

# ============================================================================
# Kronecker Products - Eigenvalue Computation
# ============================================================================

"""
    _compute_SU2_kron_eigenvalues_from_traces(c1, s1, c2, s2)

Compute eigenvalues of U1 ⊗ U2 where:
- U1 has eigenvalues e^{±iθ₁/2} with cos(θ₁/2) = c1, sin(θ₁/2) = s1
- U2 has eigenvalues e^{±iθ₂/2} with cos(θ₂/2) = c2, sin(θ₂/2) = s2
"""
function _compute_SU2_kron_eigenvalues_from_traces(c1, s1, c2, s2)
    # Real and imaginary parts before simplification
    re1, im1 = c1*c2 - s1*s2, c1*s2 + s1*c2      # e^{i(θ₁+θ₂)/2}
    re2, im2 = c1*c2 + s1*s2, c1*s2 - s1*c2      # e^{i(θ₁-θ₂)/2}
    re3, im3 = c1*c2 + s1*s2, -(c1*s2 - s1*c2)   # e^{i(-θ₁+θ₂)/2}
    re4, im4 = c1*c2 - s1*s2, -(c1*s2 + s1*c2)   # e^{i(-θ₁-θ₂)/2}
    
    λ1 = trig_simplify(re1) + im*trig_simplify(im1)
    λ2 = trig_simplify(re2) + im*trig_simplify(im2)
    λ3 = trig_simplify(re3) + im*trig_simplify(im3)
    λ4 = trig_simplify(re4) + im*trig_simplify(im4)
    
    return [simplify_eigenvalue(λ) for λ in [λ1, λ2, λ3, λ4]]
end

"""
    _try_SU2_kron_from_traces(mat, blocks)

Extract SU(2) Kronecker eigenvalues using orthogonality-based trace extraction.
"""
function _try_SU2_kron_from_traces(mat, blocks)
    B11, B12 = blocks[1]
    B21, B22 = blocks[2]
    
    tr_B11 = tr(B11)
    tr_B12 = tr(B12)
    
    norm_sq_B11 = real(tr_B11)^2 + imag(tr_B11)^2
    norm_sq_B12 = real(tr_B12)^2 + imag(tr_B12)^2
    
    sum_norm_sq = trig_simplify(Symbolics.simplify(norm_sq_B11 + norm_sq_B12))
    
    sqrt_sum = _try_symbolic_sqrt(sum_norm_sq)
    if isnothing(sqrt_sum)
        sum_simplified = aggressive_simplify(sum_norm_sq)
        sqrt_sum = _try_symbolic_sqrt(sum_simplified)
        isnothing(sqrt_sum) && return nothing
    end
    
    tr_U2_magnitude = sqrt_sum
    cos2_half = trig_simplify(tr_U2_magnitude / 2)
    
    sin2_sq = 1 - cos2_half^2
    if !(sin2_sq isa Num) && sin2_sq isa Number && real(sin2_sq) < 0
        return nothing
    end
    sin2_half = aggressive_simplify(sqrt(sin2_sq))
    
    tr_K = tr(mat)
    
    if !(tr_U2_magnitude isa Num) && tr_U2_magnitude isa Number && isapprox(tr_U2_magnitude, 0, atol=1e-10)
        return nothing
    end
    
    tr_U1 = trig_simplify(Symbolics.simplify(tr_K / tr_U2_magnitude))
    cos1_half = trig_simplify(real(tr_U1) / 2)
    
    sin1_sq = 1 - cos1_half^2
    if !(sin1_sq isa Num) && sin1_sq isa Number && real(sin1_sq) < 0
        return nothing
    end
    sin1_half = aggressive_simplify(sqrt(sin1_sq))
    
    return _compute_SU2_kron_eigenvalues_from_traces(cos1_half, sin1_half, cos2_half, sin2_half)
end
