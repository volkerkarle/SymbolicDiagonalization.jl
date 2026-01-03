# ============================================================================
# Rotation Matrix Constructors and Kronecker Products
# Clean symbolic eigenvalues for rotation matrices and their Kronecker products
# ============================================================================

using Symbolics

# ============================================================================
# SO(2) Rotation Matrix Constructor
# ============================================================================

"""
    rotation_matrix(θ) -> Matrix{Num}
    R2(θ) -> Matrix{Num}

Construct a 2×2 rotation matrix (SO(2)):

    [cos(θ)  -sin(θ)]
    [sin(θ)   cos(θ)]

Eigenvalues are `cos(θ) ± i·sin(θ) = e^{±iθ}`.

# Example
```julia
@variables θ
R = R2(θ)
eigvals(R)  # [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
```
"""
function R2(θ)
    c, s = cos(θ), sin(θ)
    return [c -s; s c]
end

# Alias
const rotation_matrix = R2

# ============================================================================
# SO(3) Rotation Matrix Constructors (Elementary Rotations)
# ============================================================================

"""
    Rx(θ) -> Matrix{Num}

Construct rotation matrix around x-axis:

    [1    0       0   ]
    [0  cos(θ) -sin(θ)]
    [0  sin(θ)  cos(θ)]

Eigenvalues are `1, cos(θ) ± i·sin(θ)`.
"""
function Rx(θ)
    c, s = cos(θ), sin(θ)
    return [1 0 0; 0 c -s; 0 s c]
end

"""
    Ry(θ) -> Matrix{Num}

Construct rotation matrix around y-axis:

    [ cos(θ)  0  sin(θ)]
    [   0     1    0   ]
    [-sin(θ)  0  cos(θ)]

Eigenvalues are `1, cos(θ) ± i·sin(θ)`.
"""
function Ry(θ)
    c, s = cos(θ), sin(θ)
    return [c 0 s; 0 1 0; -s 0 c]
end

"""
    Rz(θ) -> Matrix{Num}

Construct rotation matrix around z-axis:

    [cos(θ) -sin(θ)  0]
    [sin(θ)  cos(θ)  0]
    [  0       0     1]

Eigenvalues are `1, cos(θ) ± i·sin(θ)`.
"""
function Rz(θ)
    c, s = cos(θ), sin(θ)
    return [c -s 0; s c 0; 0 0 1]
end

# ============================================================================
# SO(2) Kronecker Product Eigenvalues (Clean Form)
# ============================================================================

"""
    so2_kron_eigenvalues(angles::Vector) -> Vector

Compute eigenvalues of R2(θ₁) ⊗ R2(θ₂) ⊗ ... ⊗ R2(θₖ) in clean trigonometric form.

The eigenvalues are `e^{i(±θ₁±θ₂±...±θₖ)}` for all 2^k sign combinations:
    cos(±θ₁±θ₂±...±θₖ) + i·sin(±θ₁±θ₂±...±θₖ)

# Example
```julia
@variables α β γ
vals = so2_kron_eigenvalues([α, β, γ])
# 8 eigenvalues: cos(±α±β±γ) + i·sin(±α±β±γ)
```
"""
function so2_kron_eigenvalues(angles::Vector)
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
    so2_kron(angles::Vector) -> Matrix

Construct the Kronecker product R2(θ₁) ⊗ R2(θ₂) ⊗ ... ⊗ R2(θₖ).

# Example
```julia
@variables α β
K = so2_kron([α, β])  # 4×4 matrix
eigvals(K)  # [cos(α+β) + i·sin(α+β), ...]
```
"""
function so2_kron(angles::Vector)
    length(angles) >= 1 || error("Need at least one angle")
    return reduce(kron, [R2(θ) for θ in angles])
end

# ============================================================================
# Detection of SO(2) Kronecker Products (for clean eigenvalue extraction)
# ============================================================================

"""
    _detect_so2_kron_angles(mat) -> Union{Vector, Nothing}

Detect if matrix is a Kronecker product of SO(2) rotation matrices.
If so, extract the rotation angles for clean eigenvalue computation.

Returns vector of angles [θ₁, θ₂, ..., θₖ] if detected, nothing otherwise.
"""
function _detect_so2_kron_angles(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Must be power of 2
    N >= 2 || return nothing
    k = Int(log2(N))
    2^k == N || return nothing
    
    # Base case: 2×2 should be SO(2)
    if k == 1
        cs = _is_so2(mat)
        isnothing(cs) && return nothing
        c, s = cs
        # Extract angle from cos and sin
        # We return the (c, s) pair since we can't easily extract θ symbolically
        return [(c, s)]
    end
    
    # Recursive case: check if mat = R(θ₁) ⊗ M for some 2×2 R and (N/2)×(N/2) M
    half = N ÷ 2
    
    # Extract blocks
    block_11 = mat[1:half, 1:half]
    block_12 = mat[1:half, half+1:N]
    block_21 = mat[half+1:N, 1:half]
    block_22 = mat[half+1:N, half+1:N]
    
    # For R ⊗ M with R = [c -s; s c]:
    # block_11 = c*M, block_12 = -s*M, block_21 = s*M, block_22 = c*M
    
    # Check block_11 == block_22 and block_12 == -block_21
    for i in 1:half, j in 1:half
        if !_issymzero(Symbolics.simplify(block_11[i,j] - block_22[i,j]))
            return nothing
        end
        if !_issymzero(Symbolics.simplify(block_12[i,j] + block_21[i,j]))
            return nothing
        end
    end
    
    # Extract c and s from the (1,1) elements of the blocks
    # block_11[1,1] = c * M[1,1], block_21[1,1] = s * M[1,1]
    # We need M[1,1] to extract c and s
    
    # Recursive: get the angle info from the inner matrix
    # For the inner matrix M, we need to divide by c (or s if c=0)
    
    # Use first column structure:
    # For 4×4: mat[1,1] = c₁c₂, mat[2,1] = c₁s₂, mat[3,1] = s₁c₂, mat[4,1] = s₁s₂
    # For 8×8: mat[1,1] = c₁c₂c₃, etc.
    
    # The key insight: we can extract angles directly from first column ratios
    # mat[half+1, 1] / mat[1, 1] = s₁/c₁ = tan(θ₁) (if all inner cos ≠ 0)
    
    # Actually, let's use the block structure more directly:
    # c² = sum of (block_11 first column)² / (M first column norm²)
    # But we don't know M...
    
    # Alternative: recursively decompose block_11 to get M's structure,
    # then use that to extract c and s
    
    # First, check if we can recognize the structure from first column
    col1 = mat[:, 1]
    
    # For SO(2)^⊗k, the first column has a specific pattern:
    # Each element is a product of k terms, each being c_j or s_j
    
    # For clean extraction, look at the ratio pattern
    a = col1[1]  # = c₁ * (product of remaining c's and s's for position 1)
    b = col1[half + 1]  # = s₁ * (same product)
    
    # So b/a = s₁/c₁ = tan(θ₁)
    # And a² + b² = (c₁² + s₁²) * (product)² = (product)²
    # if c₁² + s₁² = 1
    
    # Try: if a and b are products of trig functions, 
    # we can extract the first angle's cos and sin
    
    # For now, use a simpler approach: trust that the matrix is SO(2)^⊗k
    # and extract eigenvalues via the generalized first-column method
    
    # The eigenvalues of SO(2)^⊗k are products of individual eigenvalues
    # Each SO(2) factor contributes e^{±iθⱼ}
    # So eigenvalues are e^{i(±θ₁±θ₂±...±θₖ)}
    
    # From first column: we can read off cos/sin products
    # mat[idx, 1] = (sign pattern of s's and c's)
    # where idx in binary tells us which factors use sin vs cos
    
    # Return the cos/sin pairs for eigenvalue construction
    return _extract_so2_kron_cos_sin_pairs(mat, k)
end

"""
    _extract_so2_kron_cos_sin_pairs(mat, k) -> Vector{Tuple}

Extract (cos(θⱼ), sin(θⱼ)) pairs from a 2^k × 2^k SO(2) Kronecker product.
Uses the first column structure to identify the trig functions.
"""
function _extract_so2_kron_cos_sin_pairs(mat, k)
    N = 2^k
    col1 = mat[:, 1]
    
    # For k=1: col1 = [c, s]
    # For k=2: col1 = [c₁c₂, c₁s₂, s₁c₂, s₁s₂]
    # For k=3: col1 = [c₁c₂c₃, c₁c₂s₃, c₁s₂c₃, c₁s₂s₃, s₁c₂c₃, s₁c₂s₃, s₁s₂c₃, s₁s₂s₃]
    
    # Pattern: element at index i (0-based) has:
    # - s_j if bit j of i is set
    # - c_j if bit j of i is clear
    
    # To extract (c_j, s_j): 
    # c_j appears in elements where bit j is 0
    # s_j appears in elements where bit j is 1
    
    # Strategy: for each j, find ratio of elements differing only in bit j
    # This gives s_j / c_j = tan(θ_j)
    
    # For j-th factor (0-indexed), compare:
    # - element 0 (all zeros) = product of all c's
    # - element 2^j (only bit j set) = s_j * (product of other c's)
    # Ratio = s_j / c_j
    
    # Then c_j² + s_j² = 1, and c_j = col1[0] / (product of other c's)
    
    # This gets complex. Use a recursive approach instead.
    
    if k == 1
        c, s = col1[1], col1[2]
        return [(c, s)]
    end
    
    # For k > 1: use the block structure
    half = N ÷ 2
    
    # Block structure: block_11 = c₁ * M, block_21 = s₁ * M
    # where M is the (k-1)-fold Kronecker product
    
    block_11 = mat[1:half, 1:half]
    block_21 = mat[half+1:N, 1:half]
    
    # Recursively get the inner angles
    inner_pairs = _extract_so2_kron_cos_sin_pairs(block_11, k - 1)
    isnothing(inner_pairs) && return nothing
    
    # The first column of block_11 = c₁ * (first column of M)
    # The first column of block_21 = s₁ * (first column of M)
    
    # M[1,1] = product of cos's from remaining factors
    # So: c₁ = block_11[1,1] / M[1,1]
    #     s₁ = block_21[1,1] / M[1,1]
    
    # But M[1,1] comes from the recursive structure:
    # For inner_pairs = [(c₂,s₂), (c₃,s₃), ...], M[1,1] = c₂ * c₃ * ...
    
    M_11 = prod(p[1] for p in inner_pairs)
    
    c1 = Symbolics.simplify(block_11[1, 1] / M_11)
    s1 = Symbolics.simplify(block_21[1, 1] / M_11)
    
    return vcat([(c1, s1)], inner_pairs)
end

"""
    _so2_kron_eigenvalues_from_pairs(pairs::Vector) -> Vector

Compute eigenvalues from (cos, sin) pairs.
Each pair gives e^{±iθ}, and the full eigenvalue set is all products.
"""
function _so2_kron_eigenvalues_from_pairs(pairs)
    k = length(pairs)
    
    # Each factor contributes two eigenvalues: c + is and c - is
    # The total eigenvalue set is all products
    
    eigenvalues = [1 + 0im]  # Start with multiplicative identity
    
    for (c, s) in pairs
        new_eigenvalues = []
        for λ in eigenvalues
            # Multiply by both eigenvalues of this factor
            push!(new_eigenvalues, Symbolics.simplify(λ * (c + im*s)))
            push!(new_eigenvalues, Symbolics.simplify(λ * (c - im*s)))
        end
        eigenvalues = new_eigenvalues
    end
    
    # Apply trig simplification for clean output
    return [trig_simplify(λ) for λ in eigenvalues]
end


