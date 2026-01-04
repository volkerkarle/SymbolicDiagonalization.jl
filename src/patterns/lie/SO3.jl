# ============================================================================
# SO(3) - Special Orthogonal Group in 3D (3D Rotations)
# ============================================================================
#
# This file contains all SO(3) related functionality:
# - Constructors: SO3_Rx(θ), SO3_Ry(θ), SO3_Rz(θ)
# - Detection: _is_SO3(A)
# - Eigenvalues: _SO3_eigenvalues(A)
# - Kronecker products: detection and eigenvalue computation
#
# SO(3) is the group of 3×3 rotation matrices with eigenvalues:
#   {1, e^{iθ}, e^{-iθ}} where cos(θ) = (tr(A) - 1) / 2
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Constructors
# ============================================================================

"""
    SO3_Rx(θ) -> Matrix{Num}

Construct rotation matrix around x-axis:

    [1    0       0   ]
    [0  cos(θ) -sin(θ)]
    [0  sin(θ)  cos(θ)]

Eigenvalues are `1, cos(θ) ± i·sin(θ)`.
"""
function SO3_Rx(θ)
    c, s = cos(θ), sin(θ)
    return [1 0 0; 0 c -s; 0 s c]
end

"""
    SO3_Ry(θ) -> Matrix{Num}

Construct rotation matrix around y-axis:

    [ cos(θ)  0  sin(θ)]
    [   0     1    0   ]
    [-sin(θ)  0  cos(θ)]

Eigenvalues are `1, cos(θ) ± i·sin(θ)`.
"""
function SO3_Ry(θ)
    c, s = cos(θ), sin(θ)
    return [c 0 s; 0 1 0; -s 0 c]
end

"""
    SO3_Rz(θ) -> Matrix{Num}

Construct rotation matrix around z-axis:

    [cos(θ) -sin(θ)  0]
    [sin(θ)  cos(θ)  0]
    [  0       0     1]

Eigenvalues are `1, cos(θ) ± i·sin(θ)`.
"""
function SO3_Rz(θ)
    c, s = cos(θ), sin(θ)
    return [c -s 0; s c 0; 0 0 1]
end

# ============================================================================
# Detection
# ============================================================================

"""
    _is_SO3(A)

Check if A is a 3×3 rotation matrix (SO(3)).
Returns true if A ∈ SO(3), false otherwise.
"""
function _is_SO3(A)
    size(A) == (3, 3) || return false
    return !isnothing(_is_special_orthogonal(A))
end

"""
    _is_SO3_trig(A)

Check if A is an SO(3) matrix using trig-aware simplification.
This is useful for composed rotations like Euler angles.
"""
function _is_SO3_trig(A)
    size(A) == (3, 3) || return false
    
    # Check orthogonality with trig simplification
    ATA = A' * A
    n = 3
    for i in 1:n, j in 1:n
        expected = i == j ? 1 : 0
        diff = trig_simplify(ATA[i, j] - expected)
        if !_issymzero(diff)
            return false
        end
    end
    
    # Check det = 1
    d = trig_simplify(det(A))
    if !_issymzero(d - 1)
        return false
    end
    
    return true
end

# ============================================================================
# Eigenvalues
# ============================================================================

"""
    _SO3_eigenvalues(A)

Compute eigenvalues of an SO(3) matrix (3D rotation).

Using Rodrigues' formula:
- One eigenvalue is always 1 (the rotation axis)
- Other two are e^{±iθ} where cos(θ) = (tr(A) - 1) / 2
"""
function _SO3_eigenvalues(A)
    _is_SO3(A) || return nothing
    
    tr_A = tr(A)
    cos_theta = (tr_A - 1) / 2
    
    # sin(θ) = sqrt(1 - cos²(θ))
    sin_theta_sq = 1 - cos_theta^2
    sin_theta = aggressive_simplify(sqrt(sin_theta_sq))
    
    lambda1 = 1
    lambda2 = cos_theta + im * sin_theta
    lambda3 = cos_theta - im * sin_theta
    
    return [lambda1, simplify_eigenvalue(lambda2), simplify_eigenvalue(lambda3)]
end

# ============================================================================
# Eigenvectors
# ============================================================================

"""
    _SO3_axis_type(A)

Detect if A is a rotation around a principal axis (x, y, or z).

Returns:
- :Rx if rotation around x-axis: [1 0 0; 0 c -s; 0 s c]
- :Ry if rotation around y-axis: [c 0 s; 0 1 0; -s 0 c]
- :Rz if rotation around z-axis: [c -s 0; s c 0; 0 0 1]
- nothing otherwise (general rotation)

For axis-aligned rotations, we have closed-form eigenvectors.
For general rotations, eigenvector computation requires symbolic nullspace.
"""
function _SO3_axis_type(A)
    size(A) == (3, 3) || return nothing
    
    # Check for Rz: [c -s 0; s c 0; 0 0 1]
    if _issymzero(A[1, 3]) && _issymzero(A[2, 3]) && 
       _issymzero(A[3, 1]) && _issymzero(A[3, 2]) && _issymzero(A[3, 3] - 1)
        # Verify 2x2 block is SO(2)
        if !isnothing(_is_SO2(A[1:2, 1:2]))
            return :Rz
        end
    end
    
    # Check for Rx: [1 0 0; 0 c -s; 0 s c]
    if _issymzero(A[1, 1] - 1) && _issymzero(A[1, 2]) && _issymzero(A[1, 3]) &&
       _issymzero(A[2, 1]) && _issymzero(A[3, 1])
        # Verify 2x2 block is SO(2)
        if !isnothing(_is_SO2(A[2:3, 2:3]))
            return :Rx
        end
    end
    
    # Check for Ry: [c 0 s; 0 1 0; -s 0 c]
    if _issymzero(A[2, 2] - 1) && _issymzero(A[1, 2]) && _issymzero(A[3, 2]) &&
       _issymzero(A[2, 1]) && _issymzero(A[2, 3])
        # Verify it's actually Ry structure: A[1,3] = sin(θ), A[3,1] = -sin(θ)
        if _issymzero(A[1, 1] - A[3, 3]) && _issymzero(A[1, 3] + A[3, 1])
            return :Ry
        end
    end
    
    return nothing
end

"""
    _SO3_axis_eigenvectors(axis_type)

Return the eigenvectors for an axis-aligned SO(3) rotation.

For eigenvalues {1, e^{-iθ}, e^{iθ}}:

Rz (rotation around z): 
  - λ=1: [0, 0, 1] (z-axis)
  - λ=e^{-iθ}: [1, i, 0]/√2 
  - λ=e^{iθ}: [1, -i, 0]/√2

Rx (rotation around x):
  - λ=1: [1, 0, 0] (x-axis)
  - λ=e^{-iθ}: [0, 1, i]/√2
  - λ=e^{iθ}: [0, 1, -i]/√2

Ry (rotation around y):
  - λ=1: [0, 1, 0] (y-axis)
  - λ=e^{-iθ}: [1, 0, i]/√2
  - λ=e^{iθ}: [1, 0, -i]/√2

Note: Returns unnormalized eigenvectors (without 1/√2 factor).
"""
function _SO3_axis_eigenvectors(axis_type)
    if axis_type == :Rz
        # [1, i, 0] for e^{-iθ}, [1, -i, 0] for e^{iθ}
        return [[0, 0, 1], [1, im, 0], [1, -im, 0]]
    elseif axis_type == :Rx
        # [0, 1, i] for e^{-iθ}, [0, 1, -i] for e^{iθ}
        return [[1, 0, 0], [0, 1, im], [0, 1, -im]]
    elseif axis_type == :Ry
        # [1, 0, i] for e^{-iθ}, [1, 0, -i] for e^{iθ}
        return [[0, 1, 0], [1, 0, im], [1, 0, -im]]
    end
    return nothing
end

"""
    _SO3_eigenpairs(A)

Compute eigenvalue-eigenvector pairs for an SO(3) matrix.

For axis-aligned rotations (Rx, Ry, Rz), returns closed-form eigenpairs.
For general rotations, returns nothing (requires nullspace computation).

Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing.

# Example
```julia
@variables θ
R = SO3_Rz(θ)
pairs = _SO3_eigenpairs(R)
# pairs[1] = (1, [[0, 0, 1]])
# pairs[2] = (cos(θ) - im*sin(θ), [[1, im, 0]])
# pairs[3] = (cos(θ) + im*sin(θ), [[1, -im, 0]])
```
"""
function _SO3_eigenpairs(A)
    size(A) == (3, 3) || return nothing
    
    # First check axis type - this is fast and doesn't require full SO(3) verification
    axis_type = _SO3_axis_type(A)
    
    # For general rotations, we don't have closed-form eigenvectors
    # The rotation axis (eigenvector for λ=1) would need to be computed
    # from the antisymmetric part of A, which is more complex symbolically
    isnothing(axis_type) && return nothing
    
    # Verify it's actually SO(3) using trig-aware check
    _is_SO3(A) || _is_SO3_trig(A) || return nothing
    
    # Compute eigenvalues
    tr_A = tr(A)
    cos_theta = (tr_A - 1) / 2
    sin_theta_sq = 1 - cos_theta^2
    sin_theta = aggressive_simplify(sqrt(sin_theta_sq))
    
    # Eigenvalues: 1, e^{-iθ}, e^{iθ}
    lambda1 = 1
    lambda2 = simplify_eigenvalue(cos_theta - im * sin_theta)  # e^{-iθ}
    lambda3 = simplify_eigenvalue(cos_theta + im * sin_theta)  # e^{iθ}
    vals = [lambda1, lambda2, lambda3]
    
    # Get eigenvectors for this axis type
    vecs = _SO3_axis_eigenvectors(axis_type)
    isnothing(vecs) && return nothing
    
    # Return as vector of (eigenvalue, [eigenvector]) tuples
    return [(vals[1], [vecs[1]]), (vals[2], [vecs[2]]), (vals[3], [vecs[3]])]
end

# ============================================================================
# Kronecker Products - Detection
# ============================================================================

"""
    _detect_SO3_kronecker_product(mat)

Detect if mat is a Kronecker product of SO(3) rotation matrices.

For a 9×9 matrix, check if it's R₁ ⊗ R₂ where R₁, R₂ ∈ SO(3).
For a 27×27 matrix, check if it's R₁ ⊗ R₂ ⊗ R₃ where Rᵢ ∈ SO(3).

Returns a vector of eigenvalues if detected, nothing otherwise.
"""
function _detect_SO3_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Only for symbolic matrices
    _is_numeric_matrix(mat) && return nothing
    
    # Check if dimension is a power of 3 (9, 27, 81, ...)
    N < 9 && return nothing
    log3_N = log(3, N)
    isinteger(round(log3_N)) || return nothing
    k = Int(round(log3_N))
    k < 2 && return nothing  # Need at least 3⊗3
    
    if k == 2
        return _detect_SO3_kron_2fold(mat)
    elseif k == 3
        return _detect_SO3_kron_3fold(mat)
    end
    
    return nothing
end

"""
    _detect_SO3_kron_2fold(mat)

Detect SO(3) ⊗ SO(3) structure in a 9×9 matrix.
"""
function _detect_SO3_kron_2fold(mat)
    size(mat) == (9, 9) || return nothing
    
    # Extract 3×3 blocks: blocks[i][j] = R1[i,j] * R2 for R1 ⊗ R2
    blocks = [[mat[3(i-1)+1:3i, 3(j-1)+1:3j] for j in 1:3] for i in 1:3]
    
    # Strategy 1: Check diagonal blocks for SO(3) structure
    for diag_idx in [3, 2, 1]  # Try in order: Rz, Ry, Rx
        Bkk = blocks[diag_idx][diag_idx]
        if _is_SO3_trig(Bkk)
            R2_candidate = Bkk
            
            ref_idx = diag_idx
            R2_ref = R2_candidate[ref_idx, ref_idx]
            
            if _issymzero(R2_ref)
                R2_ref = R2_candidate[1, 1]
                ref_idx = 1
            end
            
            if !_issymzero(R2_ref)
                R1_entries = [Symbolics.simplify(blocks[i][j][ref_idx, ref_idx] / R2_ref) for i in 1:3, j in 1:3]
                R1_candidate = R1_entries
                
                if _is_SO3_trig(R1_candidate)
                    return _compute_SO3_kron_eigenvalues(R1_candidate, R2_candidate)
                end
            end
        end
    end
    
    # Strategy 2: Try general Kronecker factorization
    kron_info = _try_kronecker_factorization(mat, 3, 3)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        A_simplified = Symbolics.simplify.(A)
        B_simplified = Symbolics.simplify.(B)
        
        if _is_SO3_trig(A_simplified) && _is_SO3_trig(B_simplified)
            return _compute_SO3_kron_eigenvalues(A_simplified, B_simplified)
        end
        
        # Check for scaled orthogonal matrices
        AtA = Symbolics.simplify.(A' * A)
        BtB = Symbolics.simplify.(B' * B)
        
        A_scale_sq = _check_scalar_identity(AtA)
        B_scale_sq = _check_scalar_identity(BtB)
        
        if !isnothing(A_scale_sq) && !isnothing(B_scale_sq)
            # Can still use trace-based approach
        end
    end
    
    # Strategy 3: Use trace-based computation
    result = _try_SO3_kron_from_traces(mat, blocks)
    if !isnothing(result)
        return result
    end
    
    return nothing
end

"""
    _detect_SO3_kron_3fold(mat)

Detect SO(3) ⊗ SO(3) ⊗ SO(3) structure in a 27×27 matrix.
"""
function _detect_SO3_kron_3fold(mat)
    size(mat) == (27, 27) || return nothing
    
    # First try to factor as 3 ⊗ 9
    kron_info = _try_kronecker_factorization(mat, 3, 9)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        A_simplified = Symbolics.simplify.(A)
        if _is_SO3_trig(A_simplified)
            vals_A = _SO3_eigenvalues(A_simplified)
            vals_B = _detect_SO3_kron_2fold(B)
            
            if !isnothing(vals_A) && !isnothing(vals_B)
                eigenvalues = []
                for λ in vals_A, μ in vals_B
                    push!(eigenvalues, simplify_eigenvalue(Symbolics.simplify(λ * μ)))
                end
                return eigenvalues
            end
        end
    end
    
    # Try to factor as 9 ⊗ 3
    kron_info = _try_kronecker_factorization(mat, 9, 3)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        B_simplified = Symbolics.simplify.(B)
        if _is_SO3_trig(B_simplified)
            vals_B = _SO3_eigenvalues(B_simplified)
            vals_A = _detect_SO3_kron_2fold(A)
            
            if !isnothing(vals_A) && !isnothing(vals_B)
                eigenvalues = []
                for λ in vals_A, μ in vals_B
                    push!(eigenvalues, simplify_eigenvalue(Symbolics.simplify(λ * μ)))
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
    _compute_SO3_kron_eigenvalues(R1, R2)

Compute eigenvalues of R1 ⊗ R2 where R1, R2 ∈ SO(3).

For Kronecker product, eigenvalues are all 9 products of individual eigenvalues.
"""
function _compute_SO3_kron_eigenvalues(R1, R2)
    vals_R1 = _SO3_eigenvalues(R1)
    vals_R2 = _SO3_eigenvalues(R2)
    
    if isnothing(vals_R1) || isnothing(vals_R2)
        return nothing
    end
    
    eigenvalues = []
    for λ in vals_R1, μ in vals_R2
        product = λ * μ
        if product isa Complex
            re = Symbolics.simplify(real(product))
            im_part = Symbolics.simplify(imag(product))
            push!(eigenvalues, re + im * im_part)
        else
            push!(eigenvalues, Symbolics.simplify(product))
        end
    end
    return eigenvalues
end

"""
    _compute_SO3_kron_eigenvalues_from_traces(tr_R1, tr_R2)

Compute eigenvalues of SO(3) ⊗ SO(3) from the traces of the factors.
"""
function _compute_SO3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
    cos_theta1 = Symbolics.simplify((tr_R1 - 1) / 2)
    cos_theta2 = Symbolics.simplify((tr_R2 - 1) / 2)
    
    sin_theta1_sq = Symbolics.simplify(1 - cos_theta1^2)
    sin_theta2_sq = Symbolics.simplify(1 - cos_theta2^2)
    
    sin_theta1 = sqrt(sin_theta1_sq)
    sin_theta2 = sqrt(sin_theta2_sq)
    
    λ = [1, cos_theta1 + im * sin_theta1, cos_theta1 - im * sin_theta1]
    μ = [1, cos_theta2 + im * sin_theta2, cos_theta2 - im * sin_theta2]
    
    eigenvalues = []
    for l in λ, m in μ
        product = l * m
        if product isa Complex
            re = Symbolics.simplify(real(product))
            im_part = Symbolics.simplify(imag(product))
            push!(eigenvalues, re + im * im_part)
        else
            push!(eigenvalues, Symbolics.simplify(product))
        end
    end
    
    return eigenvalues
end

"""
    _try_SO3_kron_from_traces(mat, blocks)

Try to compute SO(3) ⊗ SO(3) eigenvalues using trace information.
"""
function _try_SO3_kron_from_traces(mat, blocks)
    # tr(K) from the matrix
    tr_K = Symbolics.simplify(sum(mat[i, i] for i in 1:9))
    
    # Strategy A: Find a diagonal block that is directly SO(3)
    for k in 1:3
        Bkk = blocks[k][k]
        if _is_SO3_trig(Bkk)
            tr_R2 = Symbolics.simplify(sum(Bkk[i, i] for i in 1:3))
            tr_R1 = Symbolics.simplify(tr_K / tr_R2)
            return _compute_SO3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
        end
    end
    
    # Strategy B: Use determinant to extract scale factor
    for k in 1:3
        Bkk = blocks[k][k]
        det_Bkk = trig_simplify(det(Bkk))
        
        if _issymzero(det_Bkk)
            continue
        end
        
        c_k = _try_symbolic_cbrt(det_Bkk)
        if !isnothing(c_k)
            tr_Bkk = Symbolics.simplify(sum(Bkk[i, i] for i in 1:3))
            tr_R2 = trig_simplify(tr_Bkk / c_k)
            tr_R1 = trig_simplify(tr_K / tr_R2)
            return _compute_SO3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
        end
    end
    
    # Strategy C: Use orthogonality to extract tr(R₂)² from block traces
    tr_row1 = [Symbolics.simplify(sum(blocks[1][j][i, i] for i in 1:3)) for j in 1:3]
    tr_R2_sq = trig_simplify(sum(t^2 for t in tr_row1))
    tr_R2 = sqrt(tr_R2_sq)
    tr_R1 = trig_simplify(tr_K / tr_R2)
    
    return _compute_SO3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
end

# ============================================================================
# Helper Functions (used by Kronecker detection)
# ============================================================================

"""
    _check_scalar_identity(M)

Check if M is a scalar multiple of the identity matrix.
Returns the scalar if so, nothing otherwise.
"""
function _check_scalar_identity(M)
    n = size(M, 1)
    size(M, 2) == n || return nothing
    
    diag_val = Symbolics.simplify(M[1, 1])
    
    for i in 2:n
        if !_issymzero_trig(M[i, i] - diag_val)
            return nothing
        end
    end
    
    for i in 1:n, j in 1:n
        if i != j && !_issymzero_trig(M[i, j])
            return nothing
        end
    end
    
    return diag_val
end

"""
    _try_symbolic_sqrt(expr)

Try to compute sqrt(expr) symbolically if expr is a perfect square.
"""
function _try_symbolic_sqrt(expr)
    unwrapped = Symbolics.unwrap(expr)
    
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (^)
        args = Symbolics.arguments(unwrapped)
        if length(args) == 2 && isequal(Symbolics.value(args[2]), 2)
            return Num(args[1])
        end
    end
    
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (*)
        args = Symbolics.arguments(unwrapped)
        sqrt_factors = []
        for arg in args
            if Symbolics.iscall(arg) && Symbolics.operation(arg) === (^)
                sub_args = Symbolics.arguments(arg)
                if length(sub_args) == 2 && isequal(Symbolics.value(sub_args[2]), 2)
                    push!(sqrt_factors, sub_args[1])
                else
                    return nothing
                end
            else
                val = Symbolics.value(arg)
                if val isa Number && isreal(val) && real(val) >= 0
                    s = sqrt(real(val))
                    if s^2 == real(val)
                        push!(sqrt_factors, s)
                    else
                        return nothing
                    end
                else
                    return nothing
                end
            end
        end
        return Num(prod(sqrt_factors))
    end
    
    if expr isa Real && !(expr isa Num)
        if expr >= 0
            s = sqrt(expr)
            if s^2 == expr
                return s
            end
        end
    end
    
    return nothing
end

"""
    _try_symbolic_cbrt(expr)

Try to compute cbrt(expr) (cube root) symbolically if expr is a perfect cube.
"""
function _try_symbolic_cbrt(expr)
    unwrapped = Symbolics.unwrap(expr)
    
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (^)
        args = Symbolics.arguments(unwrapped)
        if length(args) == 2 && isequal(Symbolics.value(args[2]), 3)
            return Num(args[1])
        end
    end
    
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (*)
        args = Symbolics.arguments(unwrapped)
        cbrt_factors = []
        for arg in args
            if Symbolics.iscall(arg) && Symbolics.operation(arg) === (^)
                sub_args = Symbolics.arguments(arg)
                if length(sub_args) == 2 && isequal(Symbolics.value(sub_args[2]), 3)
                    push!(cbrt_factors, sub_args[1])
                else
                    return nothing
                end
            else
                return nothing
            end
        end
        return Num(prod(cbrt_factors))
    end
    
    if expr isa Real && !(expr isa Num)
        c = cbrt(expr)
        if c^3 ≈ expr
            return c
        end
    end
    
    return nothing
end
