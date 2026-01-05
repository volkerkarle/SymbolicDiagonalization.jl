# ============================================================================
# Finite Group: Quaternion Group Q₈
# ============================================================================
#
# The quaternion group Q₈ = {±1, ±i, ±j, ±k} has 8 elements with:
#   i² = j² = k² = ijk = -1
#
# Q₈ is non-abelian but has simple representation theory:
#   - 4 one-dimensional irreps (on quotient Q₈/{±1} ≅ Z₂ × Z₂)
#   - 1 two-dimensional irrep (the standard quaternion representation)
#
# Character table of Q₈:
#   Class:   {1}  {-1}  {±i}  {±j}  {±k}
#   χ₁:       1     1     1     1     1    (trivial)
#   χ₂:       1     1     1    -1    -1
#   χ₃:       1     1    -1     1    -1
#   χ₄:       1     1    -1    -1     1
#   χ₅:       2    -2     0     0     0    (2D irrep)
#
# For Q₈-invariant matrices (commuting with the regular representation),
# eigenvalues can be computed from character sums with multiplicity from
# the character table.
#
# Common applications:
#   - Spin systems with quaternionic structure
#   - Signal processing with quaternion filters
#   - 4D rotation decomposition
# ============================================================================

"""
    _quaternion_unit_matrices()

Return the 2×2 matrix representations of quaternion units {1, i, j, k}.

Using the standard representation:
- 1 = [1 0; 0 1]
- i = [i 0; 0 -i]  (using complex i)
- j = [0 1; -1 0]
- k = [0 i; i 0]

These satisfy i² = j² = k² = ijk = -I.
"""
function _quaternion_unit_matrices()
    I2 = Complex{Int}[1 0; 0 1]
    Qi = Complex{Int}[im 0; 0 -im]
    Qj = Complex{Int}[0 1; -1 0]
    Qk = Complex{Int}[0 im; im 0]
    return (I2, Qi, Qj, Qk)
end

"""
    _is_quaternion_matrix(mat)

Check if a 2×2 matrix can be written as a quaternion: a·1 + b·i + c·j + d·k.

A quaternion matrix has the form:
```
[a + bi    c + di ]
[-c + di   a - bi ]
```
where a, b, c, d are real (or symbolic) scalars.

Returns (a, b, c, d) if the matrix is a quaternion, `nothing` otherwise.
"""
function _is_quaternion_matrix(mat)
    size(mat) == (2, 2) || return nothing
    
    # Extract elements
    m11 = mat[1, 1]  # a + bi
    m12 = mat[1, 2]  # c + di
    m21 = mat[2, 1]  # -c + di
    m22 = mat[2, 2]  # a - bi
    
    # Check quaternion structure: m11 = conj(m22), m12 = -conj(m21)
    # For symbolic: real(m11) = real(m22), imag(m11) = -imag(m22)
    #               real(m12) = -real(m21), imag(m12) = imag(m21)
    
    if !_issymzero(m11 + m22 - 2*real(m11)) || !_issymzero(m11 - m22 - 2im*imag(m11))
        # Check m22 = conj(m11) = real(m11) - i*imag(m11)
        if !_issymzero(real(m11) - real(m22)) || !_issymzero(imag(m11) + imag(m22))
            return nothing
        end
    end
    
    if !_issymzero(m12 + m21 - 2im*imag(m12)) || !_issymzero(m12 - m21 - 2*real(m12))
        # Check m21 = -conj(m12) = -real(m12) + i*imag(m12)
        if !_issymzero(real(m12) + real(m21)) || !_issymzero(imag(m12) - imag(m21))
            return nothing
        end
    end
    
    # Extract quaternion components
    a = real(m11)
    b = imag(m11)
    c = real(m12)
    d = imag(m12)
    
    return (a, b, c, d)
end

"""
    _quaternion_eigenvalues(a, b, c, d)

Compute eigenvalues of a quaternion matrix q = a + bi + cj + dk.

For a quaternion with norm |q| = √(a² + b² + c² + d²), the eigenvalues are:
    λ₁ = a + √(b² + c² + d²)·i = a + |Im(q)|·i
    λ₂ = a - √(b² + c² + d²)·i = a - |Im(q)|·i

These are complex conjugates with real part a and imaginary part ±|Im(q)|.

For a pure quaternion (a = 0): eigenvalues are ±|q|·i (purely imaginary).
For a real scalar (b = c = d = 0): eigenvalue is a (with multiplicity 2).
"""
function _quaternion_eigenvalues(a, b, c, d)
    # |Im(q)| = √(b² + c² + d²)
    im_norm_sq = b^2 + c^2 + d^2
    
    # Check if it's a real scalar
    if _issymzero(im_norm_sq)
        return [a, a]
    end
    
    im_norm = sqrt(im_norm_sq)
    λ1 = a + im*im_norm
    λ2 = a - im*im_norm
    
    return [Symbolics.simplify(λ1), Symbolics.simplify(λ2)]
end

"""
    _is_block_quaternion(mat)

Check if a matrix is block-diagonal with each block being a quaternion matrix.

For n×n matrix where n is even, checks if:
- Matrix is block diagonal with 2×2 blocks
- Each 2×2 block is a valid quaternion matrix

Returns vector of (a, b, c, d) tuples for each block, or `nothing`.
"""
function _is_block_quaternion(mat)
    n = size(mat, 1)
    n == size(mat, 2) || return nothing
    n % 2 == 0 || return nothing
    n_blocks = div(n, 2)
    
    # Check block diagonal structure
    for i in 1:n_blocks
        row_start = 2*(i-1) + 1
        row_end = 2*i
        
        for j in 1:n_blocks
            col_start = 2*(j-1) + 1
            col_end = 2*j
            
            block = mat[row_start:row_end, col_start:col_end]
            
            if i == j
                # Diagonal block: should be quaternion
                quat = _is_quaternion_matrix(block)
                if isnothing(quat)
                    return nothing
                end
            else
                # Off-diagonal block: should be zero
                if !all(_issymzero, block)
                    return nothing
                end
            end
        end
    end
    
    # Extract quaternion parameters for each block
    quaternions = Vector{NTuple{4, Any}}()
    for i in 1:n_blocks
        row_start = 2*(i-1) + 1
        row_end = 2*i
        block = mat[row_start:row_end, row_start:row_end]
        quat = _is_quaternion_matrix(block)
        push!(quaternions, quat)
    end
    
    return quaternions
end

"""
    _block_quaternion_eigenvalues(quaternions)

Compute eigenvalues of a block-diagonal quaternion matrix.

Given vector of (a, b, c, d) tuples for each 2×2 quaternion block,
returns all eigenvalues.
"""
function _block_quaternion_eigenvalues(quaternions)
    all_eigenvalues = Vector{Any}()
    
    for (a, b, c, d) in quaternions
        block_eigs = _quaternion_eigenvalues(a, b, c, d)
        append!(all_eigenvalues, block_eigs)
    end
    
    return all_eigenvalues
end

"""
    _Q8_group_elements()

Return the 8 elements of Q₈ in order: 1, -1, i, -i, j, -j, k, -k
Each element is represented as (a, b, c, d) where q = a·1 + b·i + c·j + d·k
"""
function _Q8_group_elements()
    return [
        (1, 0, 0, 0),   # 1
        (-1, 0, 0, 0),  # -1
        (0, 1, 0, 0),   # i
        (0, -1, 0, 0),  # -i
        (0, 0, 1, 0),   # j
        (0, 0, -1, 0),  # -j
        (0, 0, 0, 1),   # k
        (0, 0, 0, -1),  # -k
    ]
end

"""
    _Q8_multiplication_table()

Return the Q₈ group multiplication table as an 8×8 matrix.
Entry (i,j) gives the index of the product g_i · g_j in the element list.

Element ordering: 1, -1, i, -i, j, -j, k, -k (indices 1-8)
"""
function _Q8_multiplication_table()
    # Multiplication table for Q₈ = {1, -1, i, -i, j, -j, k, -k}
    # Using indices 1-8 for the 8 elements
    return [
        # 1   -1   i   -i   j   -j   k   -k
        1    2    3    4    5    6    7    8;   # 1·x
        2    1    4    3    6    5    8    7;   # (-1)·x
        3    4    2    1    7    8    6    5;   # i·x
        4    3    1    2    8    7    5    6;   # (-i)·x
        5    6    8    7    2    1    3    4;   # j·x
        6    5    7    8    1    2    4    3;   # (-j)·x
        7    8    5    6    4    3    2    1;   # k·x
        8    7    6    5    3    4    1    2;   # (-k)·x
    ]
end

"""
    _is_Q8_regular_representation(mat)

Check if an 8×8 matrix is a Q₈-invariant matrix in the regular representation.

A matrix M in the group algebra ℂ[Q₈] has the form:
    M = Σ_{g ∈ Q₈} c_g · R_g

where R_g is the regular representation matrix for element g (permutation matrix
from right multiplication by g⁻¹).

Such matrices are characterized by having constant values along "Q₈-orbits":
for each g ∈ Q₈, M[h, k] = M[hg, kg] for all h, k.

Returns a tuple of 8 coefficients (c₁, c₋₁, cᵢ, c₋ᵢ, cⱼ, c₋ⱼ, cₖ, c₋ₖ) if the
matrix is Q₈-invariant, or `nothing` otherwise.
"""
function _is_Q8_regular_representation(mat)
    size(mat) == (8, 8) || return nothing
    
    # In the regular representation, a Q₈-invariant matrix M satisfies:
    # M[h, k] depends only on the "quotient" h⁻¹·k
    #
    # So M[i, j] = c_{g_i⁻¹ · g_j} for all i, j
    #
    # This means we need to check that M[i,j] = M[i',j'] whenever
    # g_i⁻¹·g_j = g_{i'}⁻¹·g_{j'}
    
    mult_table = _Q8_multiplication_table()
    elements = _Q8_group_elements()
    
    # Compute inverse indices: for Q₈, inverse of (a,b,c,d) is (a,-b,-c,-d)
    # except for ±1 which are self-inverse
    inverse_idx = [1, 2, 4, 3, 6, 5, 8, 7]  # inv(1)=1, inv(-1)=-1, inv(i)=-i, etc.
    
    # Build a map from (i,j) -> index of g_i⁻¹·g_j
    quotient_idx = zeros(Int, 8, 8)
    for i in 1:8, j in 1:8
        inv_i = inverse_idx[i]
        quotient_idx[i, j] = mult_table[inv_i, j]
    end
    
    # Extract the coefficient for each group element by finding a representative entry
    # For each group element g (index 1-8), find the first (i,j) such that g_i⁻¹·g_j = g
    coeffs = Vector{Any}(undef, 8)
    for g in 1:8
        # Find first (i,j) with quotient_idx[i,j] == g
        found = false
        for i in 1:8, j in 1:8
            if quotient_idx[i, j] == g
                coeffs[g] = mat[i, j]
                found = true
                break
            end
        end
        if !found
            return nothing  # Should never happen
        end
    end
    
    # Now verify all entries match the expected pattern
    for i in 1:8, j in 1:8
        g = quotient_idx[i, j]
        expected = coeffs[g]
        if !_issymzero(mat[i, j] - expected)
            return nothing
        end
    end
    
    return tuple(coeffs...)
end

"""
    _Q8_regular_eigenvalues(coeffs)

Compute eigenvalues of a Q₈-invariant matrix given its group algebra coefficients.

For M = Σ c_g · R_g in ℂ[Q₈], the eigenvalues come from the character table of Q₈:

Character table of Q₈:
   Class:   {1}  {-1}  {±i}  {±j}  {±k}
   χ₁:       1     1     1     1     1    (trivial)
   χ₂:       1     1     1    -1    -1
   χ₃:       1     1    -1     1    -1
   χ₄:       1     1    -1    -1     1
   χ₅:       2    -2     0     0     0    (2D irrep)

For each irrep χ with dimension d, we get eigenvalue:
   λ_χ = (1/d) · Σ_{g ∈ G} c_g · χ(g)

with multiplicity d.

So:
   λ₁ = c₁ + c₋₁ + cᵢ + c₋ᵢ + cⱼ + c₋ⱼ + cₖ + c₋ₖ        (mult 1)
   λ₂ = c₁ + c₋₁ + cᵢ + c₋ᵢ - cⱼ - c₋ⱼ - cₖ - c₋ₖ        (mult 1)
   λ₃ = c₁ + c₋₁ - cᵢ - c₋ᵢ + cⱼ + c₋ⱼ - cₖ - c₋ₖ        (mult 1)
   λ₄ = c₁ + c₋₁ - cᵢ - c₋ᵢ - cⱼ - c₋ⱼ + cₖ + c₋ₖ        (mult 1)
   λ₅ = c₁ - c₋₁                                          (mult 4)
"""
function _Q8_regular_eigenvalues(coeffs)
    c1, cm1, ci, cmi, cj, cmj, ck, cmk = coeffs
    
    # 1D irreps
    λ1 = c1 + cm1 + ci + cmi + cj + cmj + ck + cmk  # trivial
    λ2 = c1 + cm1 + ci + cmi - cj - cmj - ck - cmk  # sign on j,k
    λ3 = c1 + cm1 - ci - cmi + cj + cmj - ck - cmk  # sign on i,k
    λ4 = c1 + cm1 - ci - cmi - cj - cmj + ck + cmk  # sign on i,j
    
    # 2D irrep (mult 4 because 2² = 4 in regular representation)
    # χ₅(1) = 2, χ₅(-1) = -2, χ₅(±i) = χ₅(±j) = χ₅(±k) = 0
    # λ₅ = (1/2) · (2·c₁ + (-2)·c₋₁) = c₁ - c₋₁
    λ5 = c1 - cm1
    
    # Simplify all eigenvalues
    eigenvalues = [
        Symbolics.simplify(λ1),
        Symbolics.simplify(λ2),
        Symbolics.simplify(λ3),
        Symbolics.simplify(λ4),
        Symbolics.simplify(λ5),
        Symbolics.simplify(λ5),
        Symbolics.simplify(λ5),
        Symbolics.simplify(λ5),
    ]
    
    return eigenvalues
end

"""
    Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)

Construct an 8×8 Q₈-invariant matrix from group algebra coefficients.

Returns M = c₁·R₁ + c₋₁·R₋₁ + cᵢ·Rᵢ + c₋ᵢ·R₋ᵢ + cⱼ·Rⱼ + c₋ⱼ·R₋ⱼ + cₖ·Rₖ + c₋ₖ·R₋ₖ

where R_g is the permutation matrix for right multiplication by g⁻¹.
"""
function Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)
    coeffs = (c1, cm1, ci, cmi, cj, cmj, ck, cmk)
    mult_table = _Q8_multiplication_table()
    inverse_idx = [1, 2, 4, 3, 6, 5, 8, 7]
    
    # Compute quotient indices
    quotient_idx = zeros(Int, 8, 8)
    for i in 1:8, j in 1:8
        inv_i = inverse_idx[i]
        quotient_idx[i, j] = mult_table[inv_i, j]
    end
    
    # Build matrix
    T = promote_type(typeof(c1), typeof(cm1), typeof(ci), typeof(cmi),
                     typeof(cj), typeof(cmj), typeof(ck), typeof(cmk))
    mat = Matrix{T}(undef, 8, 8)
    for i in 1:8, j in 1:8
        mat[i, j] = coeffs[quotient_idx[i, j]]
    end
    
    return mat
end

# ============================================================================
# Quaternion algebra for symbolic diagonalization
# ============================================================================

"""
    _symbolic_quaternion_norm_squared(a, b, c, d)

Compute |q|² = a² + b² + c² + d² for quaternion q = a + bi + cj + dk.
"""
function _symbolic_quaternion_norm_squared(a, b, c, d)
    return Symbolics.simplify(a^2 + b^2 + c^2 + d^2)
end

"""
    _symbolic_quaternion_inverse(a, b, c, d)

Compute the inverse of quaternion q = a + bi + cj + dk.

q⁻¹ = q̄ / |q|² = (a - bi - cj - dk) / (a² + b² + c² + d²)

Returns (a', b', c', d') such that q⁻¹ = a' + b'i + c'j + d'k.
"""
function _symbolic_quaternion_inverse(a, b, c, d)
    norm_sq = _symbolic_quaternion_norm_squared(a, b, c, d)
    
    if _issymzero(norm_sq)
        error("Cannot invert zero quaternion")
    end
    
    a_inv = Symbolics.simplify(a / norm_sq)
    b_inv = Symbolics.simplify(-b / norm_sq)
    c_inv = Symbolics.simplify(-c / norm_sq)
    d_inv = Symbolics.simplify(-d / norm_sq)
    
    return (a_inv, b_inv, c_inv, d_inv)
end

"""
    _symbolic_quaternion_multiply(q1, q2)

Multiply two quaternions q1 = (a₁, b₁, c₁, d₁) and q2 = (a₂, b₂, c₂, d₂).

Using: (a + bi + cj + dk)(a' + b'i + c'j + d'k) =
  (aa' - bb' - cc' - dd') + (ab' + ba' + cd' - dc')i +
  (ac' - bd' + ca' + db')j + (ad' + bc' - cb' + da')k

Returns (a, b, c, d) for the product.
"""
function _symbolic_quaternion_multiply(q1, q2)
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    
    a = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c = a1*c2 - b1*d2 + c1*a2 + d1*b2
    d = a1*d2 + b1*c2 - c1*b2 + d1*a2
    
    return (Symbolics.simplify(a), Symbolics.simplify(b), 
            Symbolics.simplify(c), Symbolics.simplify(d))
end
