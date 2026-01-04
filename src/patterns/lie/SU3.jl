# ============================================================================
# SU(3) - Special Unitary Group in 3D
# ============================================================================
#
# This file contains all SU(3) related functionality:
# - Gell-Mann matrices: gellmann_1 through gellmann_8
# - Constructors: SU3_diagonal, SU3_diagonal_trig
# - Detection: _is_SU3(A)
# - Eigenvalues: _SU3_eigenvalues(A)
# - Kronecker products: SU3_kron, SU3_kron_eigenvalues, detection
#
# SU(3) is the group of 3×3 special unitary matrices.
# For diagonal SU(3): eigenvalues are e^{iθ₁}, e^{iθ₂}, e^{-i(θ₁+θ₂)}
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Gell-Mann Matrices
# ============================================================================

"""
    gellmann_1() -> Matrix{Int}

Gell-Mann matrix λ₁:
    [0 1 0]
    [1 0 0]
    [0 0 0]
"""
gellmann_1() = [0 1 0; 1 0 0; 0 0 0]

"""
    gellmann_2() -> Matrix{Complex{Int}}

Gell-Mann matrix λ₂:
    [0 -i 0]
    [i  0 0]
    [0  0 0]
"""
gellmann_2() = [0 -im 0; im 0 0; 0 0 0]

"""
    gellmann_3() -> Matrix{Int}

Gell-Mann matrix λ₃:
    [1  0 0]
    [0 -1 0]
    [0  0 0]
"""
gellmann_3() = [1 0 0; 0 -1 0; 0 0 0]

"""
    gellmann_4() -> Matrix{Int}

Gell-Mann matrix λ₄:
    [0 0 1]
    [0 0 0]
    [1 0 0]
"""
gellmann_4() = [0 0 1; 0 0 0; 1 0 0]

"""
    gellmann_5() -> Matrix{Complex{Int}}

Gell-Mann matrix λ₅:
    [0 0 -i]
    [0 0  0]
    [i 0  0]
"""
gellmann_5() = [0 0 -im; 0 0 0; im 0 0]

"""
    gellmann_6() -> Matrix{Int}

Gell-Mann matrix λ₆:
    [0 0 0]
    [0 0 1]
    [0 1 0]
"""
gellmann_6() = [0 0 0; 0 0 1; 0 1 0]

"""
    gellmann_7() -> Matrix{Complex{Int}}

Gell-Mann matrix λ₇:
    [0 0  0]
    [0 0 -i]
    [0 i  0]
"""
gellmann_7() = [0 0 0; 0 0 -im; 0 im 0]

"""
    gellmann_8() -> Matrix{Float64}

Gell-Mann matrix λ₈:
    [1/√3   0     0   ]
    [  0  1/√3    0   ]
    [  0    0  -2/√3  ]
"""
gellmann_8() = [1/sqrt(3) 0 0; 0 1/sqrt(3) 0; 0 0 -2/sqrt(3)]

"""
    gellmann_matrices() -> Vector{Matrix}

Return all 8 Gell-Mann matrices [λ₁, λ₂, ..., λ₈].
"""
gellmann_matrices() = [gellmann_1(), gellmann_2(), gellmann_3(), gellmann_4(), 
                        gellmann_5(), gellmann_6(), gellmann_7(), gellmann_8()]

# ============================================================================
# Constructors
# ============================================================================

"""
    SU3_diagonal(θ₁, θ₂) -> Matrix

Construct a diagonal SU(3) matrix (maximal torus element):

    [e^{iθ₁}      0            0       ]
    [   0      e^{iθ₂}         0       ]
    [   0         0      e^{-i(θ₁+θ₂)}]

The constraint e^{i(θ₁+θ₂+θ₃)} = 1 (det=1) gives θ₃ = -(θ₁+θ₂).

Eigenvalues are: `e^{iθ₁}, e^{iθ₂}, e^{-i(θ₁+θ₂)}`.
"""
function SU3_diagonal(θ₁, θ₂)
    return Diagonal([exp(im*θ₁), exp(im*θ₂), exp(-im*(θ₁ + θ₂))])
end

"""
    SU3_diagonal_trig(θ₁, θ₂) -> Matrix

Construct a diagonal SU(3) matrix using trigonometric form:

    [cos(θ₁)+i·sin(θ₁)         0                    0           ]
    [       0          cos(θ₂)+i·sin(θ₂)            0           ]
    [       0                  0          cos(θ₁+θ₂)-i·sin(θ₁+θ₂)]

This form is better for symbolic manipulation and trig simplification.
"""
function SU3_diagonal_trig(θ₁, θ₂)
    c1, s1 = cos(θ₁), sin(θ₁)
    c2, s2 = cos(θ₂), sin(θ₂)
    c12, s12 = cos(θ₁ + θ₂), sin(θ₁ + θ₂)
    return Diagonal([c1 + im*s1, c2 + im*s2, c12 - im*s12])
end

# ============================================================================
# Detection
# ============================================================================

"""
    _is_SU3(A)

Check if A is a 3×3 special unitary matrix (SU(3)).
Returns true if A ∈ SU(3), false otherwise.
"""
function _is_SU3(A)
    size(A) == (3, 3) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _is_SU3_trig(A)

Check if A is an SU(3) matrix using trig-aware simplification.
"""
function _is_SU3_trig(A)
    size(A) == (3, 3) || return false
    
    # Check unitarity: A * A^H = I with trig simplification
    AAH = A * adjoint(A)
    for i in 1:3, j in 1:3
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
    _SU3_eigenvalues(A)

Compute eigenvalues of an SU(3) matrix using the cubic formula.

For SU(3), eigenvalues are on the unit circle with product = 1.
The characteristic polynomial is: λ³ - tr(A)·λ² + e₂·λ - 1 = 0
"""
function _SU3_eigenvalues(A)
    _is_SU3(A) || return nothing
    
    t1 = tr(A)
    t2 = tr(A * A)
    e2 = (t1^2 - t2) / 2
    
    # Characteristic polynomial: λ³ - t1·λ² + e2·λ - 1 = 0
    coeffs = [-1, e2, -t1, 1]
    return symbolic_roots(coeffs)
end

# ============================================================================
# Eigenvectors
# ============================================================================

"""
    _SU3_eigenpairs(A)

Compute eigenvalue-eigenvector pairs for an SU(3) matrix.

For diagonal SU(3) matrices, eigenvalues are directly on the diagonal
and eigenvectors are the standard basis vectors [1,0,0], [0,1,0], [0,0,1].

For general SU(3) matrices, returns nothing (requires nullspace computation).

Returns Vector{Tuple{eigenvalue, Vector{eigenvector}}} or nothing.

# Example
```julia
@variables θ₁ θ₂
U = SU3_diagonal_trig(θ₁, θ₂)
pairs = _SU3_eigenpairs(U)
# pairs[1] = (e^{iθ₁}, [[1, 0, 0]])
# pairs[2] = (e^{iθ₂}, [[0, 1, 0]])
# pairs[3] = (e^{-i(θ₁+θ₂)}, [[0, 0, 1]])
```
"""
function _SU3_eigenpairs(A)
    size(A) == (3, 3) || return nothing
    
    # Check if matrix is diagonal first (fast check)
    if !_is_diagonal_matrix(A)
        # For non-diagonal SU(3), we need nullspace computation
        # TODO: Could add special cases for other structured SU(3) matrices
        return nothing
    end
    
    # Verify it's actually SU(3) using trig-aware check
    # (The standard _is_SU3 may fail on trig expressions)
    _is_SU3(A) || _is_SU3_trig(A) || return nothing
    
    # For diagonal SU(3), eigenvalues are on the diagonal
    # and eigenvectors are standard basis vectors
    λ1 = A[1, 1]
    λ2 = A[2, 2]
    λ3 = A[3, 3]
    
    v1 = [1, 0, 0]
    v2 = [0, 1, 0]
    v3 = [0, 0, 1]
    
    return [(λ1, [v1]), (λ2, [v2]), (λ3, [v3])]
end

# ============================================================================
# Kronecker Products - Constructors
# ============================================================================

"""
    SU3_kron_eigenvalues(angles1::Tuple, angles2::Tuple) -> Vector

Compute eigenvalues of SU(3)⊗SU(3) for diagonal SU(3) matrices.

For U₁ = diag(e^{iα₁}, e^{iα₂}, e^{-i(α₁+α₂)}) and 
    U₂ = diag(e^{iβ₁}, e^{iβ₂}, e^{-i(β₁+β₂)}):

The 9 eigenvalues of U₁⊗U₂ are all products of individual eigenvalues.
"""
function SU3_kron_eigenvalues(angles1::Tuple, angles2::Tuple)
    θ₁, θ₂ = angles1
    φ₁, φ₂ = angles2
    
    # SU(3) eigenvalue phases
    λ1_phases = [θ₁, θ₂, -(θ₁ + θ₂)]
    λ2_phases = [φ₁, φ₂, -(φ₁ + φ₂)]
    
    eigenvalues = []
    for α in λ1_phases
        for β in λ2_phases
            combined_phase = α + β
            push!(eigenvalues, cos(combined_phase) + im * sin(combined_phase))
        end
    end
    
    return eigenvalues
end

"""
    SU3_kron(angles1::Tuple, angles2::Tuple) -> Matrix

Construct the Kronecker product of two diagonal SU(3) matrices.
"""
function SU3_kron(angles1::Tuple, angles2::Tuple)
    U1 = SU3_diagonal_trig(angles1...)
    U2 = SU3_diagonal_trig(angles2...)
    return kron(U1, U2)
end

# ============================================================================
# Kronecker Products - Detection
# ============================================================================

"""
    _detect_SU3_kronecker_product(mat)

Detect if a matrix is a Kronecker product of SU(3) matrices and compute eigenvalues.

Returns eigenvalues if detected as SU(3)⊗SU(3), nothing otherwise.

Key disambiguation from SO(3)⊗SO(3):
- Both produce 9×9 matrices
- SU(3)⊗SU(3) is complex unitary
- SO(3)⊗SO(3) is real orthogonal
"""
function _detect_SU3_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # SU(3)⊗SU(3) produces 9×9
    N == 9 || return nothing
    
    # Key disambiguation: SU(3)⊗SU(3) has complex entries
    if !_has_complex_entries(mat)
        return nothing
    end
    
    return _detect_SU3_kron_2fold(mat)
end

"""
    _detect_SU3_kron_2fold(mat)

Detect if 9×9 matrix is SU(3) ⊗ SU(3) and compute eigenvalues.
"""
function _detect_SU3_kron_2fold(mat)
    size(mat) == (9, 9) || return nothing
    
    # Strategy 0: Check if matrix is diagonal
    if _is_diagonal_matrix(mat)
        eigenvalues = []
        for i in 1:9
            val = mat[i, i]
            re_simplified = trig_simplify(real(val))
            im_simplified = trig_simplify(imag(val))
            push!(eigenvalues, re_simplified + im * im_simplified)
        end
        return eigenvalues
    end
    
    # Extract 3×3 blocks
    function get_block(i, j)
        ri = (i-1)*3+1 : i*3
        cj = (j-1)*3+1 : j*3
        return mat[ri, cj]
    end
    
    blocks = [get_block(i, j) for i in 1:3, j in 1:3]
    
    # Strategy A: Check if any diagonal block is directly SU(3)
    for k in 1:3
        Bkk = blocks[k, k]
        Bkk_simplified = Symbolics.simplify.(Bkk)
        if _is_SU3_trig(Bkk_simplified)
            U2 = Bkk_simplified
            tr_U2 = tr(U2)
            
            # Find a non-zero element in U₂
            ref_i, ref_j = 1, 1
            for ri in 1:3, rj in 1:3
                if !_issymzero_trig(real(U2[ri, rj])) || !_issymzero_trig(imag(U2[ri, rj]))
                    ref_i, ref_j = ri, rj
                    break
                end
            end
            
            U2_ref = U2[ref_i, ref_j]
            
            # Extract U₁ matrix elements
            U1 = Matrix{Any}(undef, 3, 3)
            for i in 1:3, j in 1:3
                U1[i, j] = Symbolics.simplify(blocks[i, j][ref_i, ref_j] / U2_ref)
            end
            
            U1_simplified = Symbolics.simplify.(U1)
            
            if _is_SU3_trig(U1_simplified)
                tr_U1 = tr(U1_simplified)
                return _compute_SU3_kron_eigenvalues_from_traces(tr_U1, tr_U2)
            end
        end
    end
    
    # Strategy B: Use determinant of diagonal blocks
    B11 = blocks[1, 1]
    det_B11 = det(B11)
    det_B11_simplified = trig_simplify(Symbolics.simplify(det_B11))
    
    cbrt_det = _try_symbolic_cbrt(det_B11_simplified)
    if !isnothing(cbrt_det)
        U1_11 = cbrt_det
        
        if !_issymzero_trig(real(U1_11)) || !_issymzero_trig(imag(U1_11))
            U2 = Symbolics.simplify.(B11 / U1_11)
            
            if _is_SU3_trig(U2)
                tr_U2 = tr(U2)
                
                ref_i, ref_j = 1, 1
                for ri in 1:3, rj in 1:3
                    if !_issymzero_trig(real(U2[ri, rj])) || !_issymzero_trig(imag(U2[ri, rj]))
                        ref_i, ref_j = ri, rj
                        break
                    end
                end
                U2_ref = U2[ref_i, ref_j]
                
                U1 = Matrix{Any}(undef, 3, 3)
                for i in 1:3, j in 1:3
                    U1[i, j] = Symbolics.simplify(blocks[i, j][ref_i, ref_j] / U2_ref)
                end
                
                if _is_SU3_trig(U1)
                    tr_U1 = tr(U1)
                    return _compute_SU3_kron_eigenvalues_from_traces(tr_U1, tr_U2)
                end
            end
        end
    end
    
    # Strategy C: Orthogonality-based trace extraction
    result = _try_SU3_kron_from_traces(mat, blocks)
    if !isnothing(result)
        return result
    end
    
    return nothing
end

# ============================================================================
# Kronecker Products - Eigenvalue Computation
# ============================================================================

"""
    _compute_SU3_kron_eigenvalues_from_traces(tr_U1, tr_U2)

Compute eigenvalues of U₁ ⊗ U₂ for SU(3) matrices from their traces.
"""
function _compute_SU3_kron_eigenvalues_from_traces(tr_U1, tr_U2)
    # For SU(3), characteristic polynomial: λ³ - tr·λ² + tr*·λ - 1 = 0
    tr_U1_conj = conj(tr_U1)
    tr_U2_conj = conj(tr_U2)
    
    coeffs_U1 = [-1, tr_U1_conj, -tr_U1, 1]
    coeffs_U2 = [-1, tr_U2_conj, -tr_U2, 1]
    
    λs = symbolic_roots(coeffs_U1)
    μs = symbolic_roots(coeffs_U2)
    
    eigenvalues = []
    for λ in λs
        for μ in μs
            push!(eigenvalues, simplify_eigenvalue(trig_simplify(Symbolics.simplify(λ * μ))))
        end
    end
    
    return eigenvalues
end

"""
    _try_SU3_kron_from_traces(mat, blocks)

Extract SU(3) Kronecker eigenvalues using orthogonality-based trace extraction.
"""
function _try_SU3_kron_from_traces(mat, blocks)
    # Compute Σⱼ |tr(B₁ⱼ)|² for first row of blocks
    sum_norm_sq = 0
    for j in 1:3
        tr_B1j = tr(blocks[1, j])
        norm_sq = real(tr_B1j)^2 + imag(tr_B1j)^2
        sum_norm_sq = sum_norm_sq + norm_sq
    end
    
    sum_norm_sq = trig_simplify(Symbolics.simplify(sum_norm_sq))
    
    sqrt_sum = _try_symbolic_sqrt(sum_norm_sq)
    if isnothing(sqrt_sum)
        sum_simplified = aggressive_simplify(sum_norm_sq)
        sqrt_sum = _try_symbolic_sqrt(sum_simplified)
        isnothing(sqrt_sum) && return nothing
    end
    
    tr_U2_magnitude = sqrt_sum
    
    if !(tr_U2_magnitude isa Num) && tr_U2_magnitude isa Number && isapprox(tr_U2_magnitude, 0, atol=1e-10)
        return nothing
    end
    
    # Check if matrix is diagonal
    is_diag = true
    for i in 1:9, j in 1:9
        if i != j
            if !_issymzero_trig(real(mat[i, j])) || !_issymzero_trig(imag(mat[i, j]))
                is_diag = false
                break
            end
        end
    end
    
    if is_diag
        eigenvalues = [simplify_eigenvalue(mat[i, i]) for i in 1:9]
        return eigenvalues
    end
    
    return nothing
end
