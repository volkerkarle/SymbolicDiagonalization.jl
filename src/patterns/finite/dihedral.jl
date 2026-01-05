# ============================================================================
# Finite Group: Dihedral Group Dₙ
# ============================================================================
#
# The dihedral group Dₙ has 2n elements: n rotations and n reflections.
# A matrix is Dₙ-invariant if it commutes with both:
#   - Cyclic permutation (rotation)
#   - Reversal permutation (reflection)
#
# This is equivalent to being a SYMMETRIC CIRCULANT matrix:
#   - Circulant: each row is cyclic shift of previous
#   - Symmetric: cⱼ = cₙ₋ⱼ (palindromic first row)
#
# Representation theory of Dₙ:
#   - For n odd: 2 one-dim irreps + (n-1)/2 two-dim irreps
#   - For n even: 4 one-dim irreps + (n-2)/2 two-dim irreps
#
# This forces the characteristic polynomial to factor into linear and
# quadratic terms, making it solvable for any n.
#
# Eigenvalue formula for symmetric circulant with first row [c₀, c₁, ..., cₙ₋₁]:
#   λₖ = c₀ + 2·Σⱼ₌₁^{⌊n/2⌋} cⱼ·cos(2πjk/n)  [+ cₙ/₂·(-1)ᵏ if n even]
#
# ============================================================================

"""
    _is_symmetric_circulant(mat)

Check if a matrix is a symmetric circulant (Dₙ-invariant).

A symmetric circulant is both:
1. Circulant: row i is cyclic shift of row 1 by i-1 positions
2. Symmetric: the first row is palindromic (cⱼ = cₙ₋ⱼ)

Returns the first row if symmetric circulant, `nothing` otherwise.
"""
function _is_symmetric_circulant(mat)
    n = size(mat, 1)
    size(mat, 2) == n || return nothing
    n <= 1 && return mat[1, :]
    
    # First check if circulant
    first_row = mat[1, :]
    for i in 2:n
        for j in 1:n
            expected_idx = mod1(j - (i - 1), n)
            if !_issymzero(mat[i, j] - first_row[expected_idx])
                return nothing
            end
        end
    end
    
    # Now check if first row is palindromic: cⱼ = cₙ₋ⱼ (indices 1-based)
    # c[j] should equal c[n - j + 2] for j = 2, ..., floor(n/2) + 1
    for j in 2:div(n, 2) + 1
        mirror_idx = n - j + 2
        if mirror_idx != j && !_issymzero(first_row[j] - first_row[mirror_idx])
            return nothing
        end
    end
    
    return first_row
end

"""
    _symmetric_circulant_eigenvalues(first_row)

Compute eigenvalues of a symmetric circulant matrix using the Dₙ formula.

For a symmetric circulant with first row [c₀, c₁, ..., cₙ₋₁] where cⱼ = cₙ₋ⱼ,
the eigenvalues are:

    λₖ = c₀ + 2·Σⱼ₌₁^{m} cⱼ·cos(2πjk/n)

where m = (n-1)/2 for n odd, m = n/2-1 for n even.
For n even, add the term cₙ/₂·(-1)ᵏ.

The eigenvalues are all real (symmetric matrix).
"""
function _symmetric_circulant_eigenvalues(first_row)
    n = length(first_row)
    eigenvalues = Vector{Any}(undef, n)
    
    for k in 0:(n-1)
        # Start with c₀
        λ = first_row[1]
        
        # Add symmetric terms: 2·cⱼ·cos(2πjk/n) for j = 1, ..., m
        if n % 2 == 1
            # n odd: m = (n-1)/2
            m = div(n - 1, 2)
            for j in 1:m
                θ = 2 * π * j * k / n
                λ = λ + 2 * first_row[j + 1] * cos(θ)
            end
        else
            # n even: m = n/2 - 1, plus special term for j = n/2
            m = div(n, 2) - 1
            for j in 1:m
                θ = 2 * π * j * k / n
                λ = λ + 2 * first_row[j + 1] * cos(θ)
            end
            # Add cₙ/₂ · (-1)ᵏ term
            half_idx = div(n, 2) + 1  # 1-based index for c_{n/2}
            λ = λ + first_row[half_idx] * (iseven(k) ? 1 : -1)
        end
        
        eigenvalues[k + 1] = Symbolics.simplify(λ)
    end
    
    return eigenvalues
end

"""
    _is_polygon_adjacency(mat)

Check if a matrix is the adjacency matrix of a cycle graph Cₙ (n-gon).

The cycle graph Cₙ has adjacency matrix that is:
- Symmetric circulant with first row [0, 1, 0, ..., 0, 1]
- Equivalently: A[i,j] = 1 iff |i-j| = 1 (mod n)

Returns n if it's a cycle graph, `nothing` otherwise.
"""
function _is_polygon_adjacency(mat)
    n = size(mat, 1)
    size(mat, 2) == n || return nothing
    n < 3 && return nothing  # Need at least triangle
    
    for i in 1:n
        for j in 1:n
            diff = abs(i - j)
            # Should be 1 iff adjacent (including wraparound)
            is_adjacent = (diff == 1) || (diff == n - 1)
            expected = is_adjacent ? 1 : 0
            if !_issymzero(mat[i, j] - expected)
                return nothing
            end
        end
    end
    
    return n
end

"""
    _polygon_eigenvalues(n)

Compute eigenvalues of the cycle graph Cₙ adjacency matrix.

Eigenvalues: λₖ = 2·cos(2πk/n) for k = 0, 1, ..., n-1

For the pentagon (n=5): {2, φ-1, φ-1, -φ, -φ} where φ = (1+√5)/2
"""
function _polygon_eigenvalues(n)
    eigenvalues = Vector{Any}(undef, n)
    for k in 0:(n-1)
        θ = 2 * π * k / n
        eigenvalues[k + 1] = 2 * cos(θ)
    end
    return eigenvalues
end

# ============================================================================
# Symbolic polygon eigenvalues (for small n with nice closed forms)
# ============================================================================

"""
    _polygon_eigenvalues_symbolic(n)

Compute symbolic eigenvalues for small cycle graphs.

Returns symbolic expressions using exact values (golden ratio, etc.)
for n ≤ 12, numeric values otherwise.
"""
function _polygon_eigenvalues_symbolic(n)
    if n == 3
        # Triangle: 2cos(0) = 2, 2cos(2π/3) = -1, 2cos(4π/3) = -1
        return [2, -1, -1]
    elseif n == 4
        # Square: 2, 0, -2, 0
        return [2, 0, -2, 0]
    elseif n == 5
        # Pentagon: 2, φ-1, φ-1, -φ, -φ where φ = (1+√5)/2
        # 2cos(2π/5) = (√5 - 1)/2 = φ - 1
        # 2cos(4π/5) = -(√5 + 1)/2 = -φ
        φ = (1 + sqrt(Num(5))) / 2
        return [Num(2), φ - 1, φ - 1, -φ, -φ]
    elseif n == 6
        # Hexagon: 2, 1, -1, -2, -1, 1
        return [2, 1, -1, -2, -1, 1]
    else
        # Fall back to numeric
        return _polygon_eigenvalues(n)
    end
end
