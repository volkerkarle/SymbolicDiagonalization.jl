# ============================================================================
# Structure Detection and Matrix Property Analysis
# Functions to detect matrix properties and exploitable structure
# ============================================================================

"""
    _detect_structure(mat)

Detect the structure of a matrix for eigenvalue computation optimization.
Returns one of: :diagonal, :triangular, :hermitian, :symmetric, :unitary, :none
"""
function _detect_structure(mat)
    _is_diagonal(mat) && return :diagonal
    _is_triangular(mat) && return :triangular
    _is_hermitian(mat) && return :hermitian
    _is_symmetric(mat) && return :symmetric
    _is_unitary(mat) && return :unitary
    return :none
end

# ============================================================================
# Matrix Property Checks
# Boolean functions to check matrix properties
# ============================================================================

"""
    _is_diagonal(mat)

Check if a matrix is diagonal (all off-diagonal elements are zero).
"""
function _is_diagonal(mat)
    m, n = size(mat)
    for i in 1:m, j in 1:n
        i == j && continue
        !_issymzero(mat[i, j]) && return false
    end
    return true
end

"""
    _is_triangular(mat)

Check if a matrix is triangular (upper or lower).
"""
function _is_triangular(mat)
    m, n = size(mat)
    upper = true
    lower = true
    for i in 1:m, j in 1:n
        if i < j && !_issymzero(mat[i, j])
            upper = false
        elseif i > j && !_issymzero(mat[i, j])
            lower = false
        end
    end
    return upper || lower
end

"""
    _is_symmetric(mat)

Check if a matrix is symmetric (A = A^T).
"""
function _is_symmetric(mat)
    m, n = size(mat)
    m == n || return false
    for i in 1:m, j in i+1:n
        !_issymzero(mat[i, j] - mat[j, i]) && return false
    end
    return true
end

"""
    _is_hermitian(mat)

Check if a matrix is Hermitian (A = A^H).
"""
function _is_hermitian(mat)
    m, n = size(mat)
    m == n || return false
    for i in 1:m, j in i+1:n
        !_issymzero(mat[i, j] - conj(mat[j, i])) && return false
    end
    return true
end

"""
    _is_unitary(mat)

Check if a matrix is unitary (A * A^H = I).
"""
function _is_unitary(mat)
    m, n = size(mat)
    m == n || return false
    U = mat * adjoint(mat)
    for i in 1:m, j in 1:n
        target = i == j ? one(eltype(U)) : zero(eltype(U))
        !_issymzero(U[i, j] - target) && return false
    end
    return true
end

# ============================================================================
# Block Structure Detection
# ============================================================================

"""
    _block_split(mat)

Detect a block-diagonal split; returns the split index or nothing.
"""
function _block_split(mat)
    n = size(mat, 1)
    n <= 1 && return nothing
    for k in 1:n-1
        if all(_issymzero, mat[1:k, k+1:n]) && all(_issymzero, mat[k+1:n, 1:k])
            return k
        end
    end
    return nothing
end

"""
    _detect_multiple_blocks(mat)

Detect all block-diagonal structure in a matrix. Returns a vector of block ranges
[(start1, end1), (start2, end2), ...] or nothing if not block-diagonal.

This generalizes _block_split to handle multiple blocks. The algorithm finds maximal
contiguous blocks where each block has no interaction (zero off-diagonal blocks) with
other blocks.

Example:
```
[A 0 0]
[0 B 0]  →  [(1, 2), (3, 3), (4, 5)]
[0 0 C]
```
"""
function _detect_multiple_blocks(mat)
    n = size(mat, 1)
    n <= 1 && return nothing
    
    # Use a greedy algorithm: find the smallest independent block starting from position 1,
    # then recursively find blocks in the remaining matrix
    blocks = Tuple{Int,Int}[]
    current_start = 1
    
    while current_start <= n
        # Find the smallest k such that mat[current_start:k, *] and mat[*, current_start:k]
        # form an independent block (no interaction with indices k+1:n)
        block_end = n  # Default: extend to the end if no boundary found
        
        for k in current_start:(n-1)
            # Check if there's a block boundary after position k
            # i.e., mat[current_start:k, k+1:n] and mat[k+1:n, current_start:k] are zero
            upper_right_zero = all(_issymzero, mat[i, j] for i in current_start:k, j in (k+1):n)
            lower_left_zero = all(_issymzero, mat[i, j] for i in (k+1):n, j in current_start:k)
            
            if upper_right_zero && lower_left_zero
                # Found a boundary! This is the end of the current block
                block_end = k
                break
            end
        end
        
        push!(blocks, (current_start, block_end))
        current_start = block_end + 1
    end
    
    # Only return if we found actual blocks (more than just the whole matrix)
    return length(blocks) > 1 ? blocks : nothing
end

# ============================================================================
# Persymmetric Structure Detection
# ============================================================================

"""
    _is_persymmetric(mat)

Check if a matrix is persymmetric: Q[i,j] == Q[n+1-j, n+1-i]
"""
function _is_persymmetric(mat)
    n = size(mat, 1)
    for i in 1:n, j in 1:n
        if !isequal(Symbolics.simplify(mat[i, j] - mat[n+1-j, n+1-i]), 0)
            return false
        end
    end
    return true
end

"""
    _persymmetric_split(mat)

For symmetric persymmetric matrices, compute a transformation that reduces
the eigenvalue problem. Returns `(P, S, D)` where:
- `P` is the transformation matrix
- `S` and `D` are the symmetric and skew-symmetric parts
- Eigenvalues of `mat` = eigenvalues of `S` ∪ eigenvalues of `D`

If the matrix is not symmetric persymmetric, returns `nothing`.
"""
function _persymmetric_split(mat)
    n = size(mat, 1)
    n % 2 != 0 && return nothing  # Only works for even dimensions
    # Check if symmetric or Hermitian (both work with this decomposition)
    if !(_is_symmetric(mat) || _is_hermitian(mat))
        return nothing
    end
    !_is_persymmetric(mat) && return nothing
    
    # For symmetric persymmetric matrices, we can use the fact that they
    # commute with the exchange matrix J. The eigenspaces of J give us
    # a block decomposition.
    
    # Exchange matrix J (anti-identity)
    J = zeros(eltype(mat), n, n)
    for i in 1:n
        J[i, n+1-i] = one(eltype(mat))
    end
    
    # Build transformation matrix P = [(I+J)/√2, (I-J)/√2]
    # But for symbolic computation, avoid √2 and just use [I+J, I-J]
    # The eigenvalues are the same, eigenvectors are just scaled
    I_n = Matrix{eltype(mat)}(I, n, n)
    
    # For a symmetric persymmetric matrix Q:
    # Q*J = J*Q (they commute)
    # So Q can be decomposed using eigenvectors of J
    
    # J has eigenvalues +1 and -1 (for even n)
    # Eigenvectors of eigenvalue +1: (e_i + e_{n+1-i})/√2
    # Eigenvectors of eigenvalue -1: (e_i - e_{n+1-i})/√2
    
    # Compute Q in the transformed basis
    # This gives us two smaller blocks
    
    # For n=4: P = [e1+e4, e2+e3, e1-e4, e2-e3]
    # In this basis, Q becomes block diagonal
    
    half = div(n, 2)
    P = zeros(eltype(mat), n, n)
    
    # First half of columns: e_i + e_{n+1-i} for i=1:half
    for i in 1:half
        P[i, i] = one(eltype(mat))
        P[n+1-i, i] = one(eltype(mat))
    end
    
    # Second half of columns: e_i - e_{n+1-i} for i=1:half  
    for i in 1:half
        P[i, half+i] = one(eltype(mat))
        P[n+1-i, half+i] = -one(eltype(mat))
    end
    
    # Transform: Q_new = P^T * Q * P
    # Note: P^T * P = 2I (since we didn't normalize by 1/√2), so the eigenvalues
    # of Q_transformed are 2× the eigenvalues of Q. We need to divide by 2.
    Q_transformed = transpose(P) * mat * P
    
    # Extract the two blocks
    block1 = Q_transformed[1:half, 1:half]
    block2 = Q_transformed[half+1:end, half+1:end]
    
    # Divide by 2 to account for P not being normalized
    # (P^T * P = 2I instead of I, so eigenvalues are scaled by 2)
    block1 = Symbolics.simplify.(block1 ./ 2)
    block2 = Symbolics.simplify.(block2 ./ 2)
    
    # For Hermitian matrices, the blocks are real but may have Complex{Num} type
    # Convert to real Num if all imaginary parts are zero
    if eltype(mat) <: Complex
        # Check if blocks are actually real and convert if so
        block1_real = all(x -> isequal(Symbolics.simplify(imag(x)), 0), block1)
        block2_real = all(x -> isequal(Symbolics.simplify(imag(x)), 0), block2)
        
        if block1_real
            block1 = real.(block1)
        end
        if block2_real
            block2 = real.(block2)
        end
    end
    
    return (block1, block2, P)
end

# ============================================================================
# Complexity Analysis
# ============================================================================

"""
    _count_symbolic_vars(A)

Count distinct symbolic variables appearing in matrix entries.
"""
function _count_symbolic_vars(A)
    vars = Set{Any}()
    for elem in A
        if elem isa Num
            try
                syms = Symbolics.get_variables(elem)
                for s in syms
                    push!(vars, s)
                end
            catch
            end
        end
    end
    return length(vars)
end

"""
    _check_complexity(A; threshold = 5, quiet = false)

Warn about complexity based on number of distinct symbolic variables.
"""
function _check_complexity(A; threshold = 5, quiet = false)
    n_vars = _count_symbolic_vars(A)
    if !quiet && n_vars > threshold
        @warn "Matrix contains $n_vars symbolic variables (threshold: $threshold). Operations may be slow or timeout with many parameters. Consider structured matrices, partial substitution, or expand=false."
    end
    return n_vars
end

# ============================================================================
# Adjugate Matrix Computation
# Helper functions for computing adjugate matrices (used in eigenvector computation)
# ============================================================================

"""
    _adjugate_vectors(M)

Compute eigenvectors using the adjugate matrix method.
Only used for small matrices (n ≤ MAX_ADJUGATE_SIZE) due to O(n!) complexity.
"""
function _adjugate_vectors(M)
    n = size(M, 1)
    # Adjugate computation has O(n!) complexity, only use for small matrices
    n <= MAX_ADJUGATE_SIZE || return Vector{Vector{eltype(M)}}()
    adj = _adjugate(M)
    vecs = Vector{Vector{eltype(M)}}()
    for j in 1:n
        col = Symbolics.simplify.(adj[:, j])
        all(_issymzero, col) && continue
        push!(vecs, col)
        # Return immediately after finding the first non-zero column to avoid
        # collecting multiple linearly dependent eigenvectors from the adjugate.
        return vecs
    end
    return vecs
end

"""
    _adjugate(M)

Compute the adjugate (classical adjoint) of matrix M.
"""
function _adjugate(M)
    n = size(M, 1)
    adj = Matrix{eltype(M)}(undef, n, n)
    if n == 1
        adj[1, 1] = one(eltype(M))
        return adj
    end
    for i in 1:n, j in 1:n
        minor_det = _minor_det(M, i, j)
        adj[j, i] = (-one(eltype(M)))^(i + j) * minor_det
    end
    return Symbolics.simplify.(adj)
end

"""
    _minor_det(M, row, col)

Compute the (row, col) minor determinant of matrix M.
"""
function _minor_det(M, row, col)
    n = size(M, 1)
    if n == 1
        return one(eltype(M))
    elseif n == 2
        r = setdiff(1:2, row)
        c = setdiff(1:2, col)
        return M[r[1], c[1]]
    elseif n == 3
        rows = setdiff(1:3, row)
        cols = setdiff(1:3, col)
        a, b = rows
        c, d = cols
        return Symbolics.simplify(M[a, c] * M[b, d] - M[a, d] * M[b, c])
    else
        error("minor_det only implemented for n ≤ 3")
    end
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    _diag_error_msg(struct_hint, found, n)

Generate an error message for non-diagonalizable matrices.
"""
function _diag_error_msg(struct_hint, found, n)
    base = "matrix is not diagonalizable: found $found independent eigenvectors for size $n"
    if struct_hint in (:hermitian, :symmetric, :unitary)
        return "$base (note: $struct_hint matrices should be diagonalizable; check symbolic simplification or assumptions)"
    end
    return base
end

"""
    _std_basis(T, n, idx)

Create a standard basis vector of type T with 1 at index idx.
"""
function _std_basis(T, n, idx)
    v = zeros(T, n)
    v[idx] = one(T)
    return v
end

"""
    _build_eigenvalue_result(vals, λ, expand)

Build the standard return tuple for eigenvalue computations.
Returns `(vals, poly, λ)` where `poly` is the characteristic polynomial.
"""
function _build_eigenvalue_result(vals, λ, expand)
    poly_raw = prod(λ .- vals)
    poly = expand ? Symbolics.expand(poly_raw) : poly_raw
    return vals, poly, λ
end
