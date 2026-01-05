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
    _is_Q8_regular_representation(mat)

Check if an 8×8 matrix is invariant under the regular representation of Q₈.

The regular representation of Q₈ acts on an 8-dimensional space (one basis
element per group element). A Q₈-invariant matrix commutes with this action.

Due to Q₈'s representation theory, such matrices decompose as:
- 4 copies of 1D trivial-like irreps
- 1 copy of 2D quaternion irrep (with multiplicity 2)

This is a specialized check for the 8×8 case.

Returns `true` if the matrix appears to be Q₈-regular-invariant.
"""
function _is_Q8_regular_representation(mat)
    size(mat) == (8, 8) || return false
    
    # For a matrix to be Q₈-regular invariant, it must have very specific
    # structure based on the group multiplication table.
    # This is complex to check in general, so we use a simpler heuristic:
    # Check if eigenvalues match the expected pattern from character theory.
    
    # Q₈ has 5 conjugacy classes, so Q₈-invariant matrices should have
    # at most 5 distinct eigenvalues (in the regular representation context).
    
    # This is a placeholder - full implementation would need group
    # algebra structure checks.
    return false
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
