# ============================================================================
# Public API - Main Entry Points
# ============================================================================

"""
    symbolic_eigenvalues(A; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)

Return `(vals, poly, λ)` where `vals` are symbolic eigenvalues (closed-form up
to degree 4), `poly` is the characteristic polynomial, and `λ` is the symbol
used. This skips eigenvector computation.

`structure` can hint matrix type: `:auto` (detect), `:hermitian`,
`:symmetric`, `:unitary`, or `:none`. Structured inputs may bypass the quartic
solver when diagonal/triangular and may influence messaging, but they do not
override the degree ≤ 4 closed-form limitation.

For matrices larger than 4×4, the function attempts to detect exploitable structure:
- **Block-diagonal**: Matrices with independent blocks are decomposed recursively
- **Persymmetric**: Symmetric/Hermitian matrices with Q[i,j] = Q[n+1-j,n+1-i] can be
  split into two half-sized blocks (works for even dimensions)
- **Diagonal/Triangular**: Eigenvalues read directly from diagonal

If no structure is found allowing decomposition into degree ≤ 4 subproblems, an error
is raised (Abel-Ruffini theorem: no general closed-form for degree ≥ 5).

`expand` controls whether to expand the characteristic polynomial for display.
`complexity_threshold` warns when the number of distinct symbolic variables
exceeds this threshold; set to `nothing` to disable.
`timeout` sets maximum computation time in seconds (default: 300, set to `nothing` to disable).
`max_terms` limits expression complexity during simplification (default: 10000).
"""
function symbolic_eigenvalues(A; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)
    λ = isnothing(var) ? _fresh_lambda() : var
    mat = Matrix(A)
    struct_hint = structure === :auto ? _detect_structure(mat) : structure
    if !isnothing(complexity_threshold)
        _check_complexity(mat; threshold = complexity_threshold)
    end
    
    # Try persymmetric splitting before block-diagonal (for even dimensions)
    if size(mat, 1) % 2 == 0 && struct_hint in (:symmetric, :hermitian)
        persym_split = _persymmetric_split(mat)
        if !isnothing(persym_split)
            block1, block2, _ = persym_split
            # Check if blocks are solvable (degree ≤ 4)
            n1, n2 = size(block1, 1), size(block2, 1)
            if n1 <= 4 && n2 <= 4
                vals1, poly1, λ = symbolic_eigenvalues(block1; var = λ, structure = :symmetric, expand = expand, complexity_threshold = nothing, timeout = timeout, max_terms = max_terms)
                vals2, poly2, _ = symbolic_eigenvalues(block2; var = λ, structure = :symmetric, expand = expand, complexity_threshold = nothing, timeout = timeout, max_terms = max_terms)
                vals = vcat(vals1, vals2)
                poly = expand ? Symbolics.expand(poly1 * poly2) : poly1 * poly2
                return vals, poly, λ
            else
                @info "Persymmetric structure detected, but blocks are too large ($(n1)×$(n1) and $(n2)×$(n2)). Closed-form solutions only exist for degrees ≤ 4."
            end
        end
    end
    
    # Check for special 5×5 tridiagonal pattern
    tridiag_pattern = _detect_special_5x5_tridiagonal(mat)
    if !isnothing(tridiag_pattern)
        a, b, d = tridiag_pattern
        # Eigenvalues: a - √(2b² + d²), a - b, a, a + b, a + √(2b² + d²)
        sqrt_term = sqrt(2*b^2 + d^2)
        vals = [a - sqrt_term, a - b, a, a + b, a + sqrt_term]
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    
    # Block-diagonal shortcut: handle each block separately to avoid high-degree polynomials
    # First try multiple block detection for better decomposition
    multiple_blocks = _detect_multiple_blocks(mat)
    if !isnothing(multiple_blocks)
        # Check if all blocks are solvable
        block_sizes = [stop - start + 1 for (start, stop) in multiple_blocks]
        all_solvable = all(s -> s <= 4, block_sizes)
        
        if all_solvable
            # Solve each block recursively
            all_vals = []
            all_polys = []
            for (start, stop) in multiple_blocks
                block = mat[start:stop, start:stop]
                vals_block, poly_block, λ = symbolic_eigenvalues(block; var = λ, structure = struct_hint, expand = expand, complexity_threshold = nothing, timeout = timeout, max_terms = max_terms)
                append!(all_vals, vals_block)
                push!(all_polys, poly_block)
            end
            poly = expand ? Symbolics.expand(prod(all_polys)) : prod(all_polys)
            return all_vals, poly, λ
        else
            # Found blocks but some are too large
            max_size = maximum(block_sizes)
            @warn "Block-diagonal structure detected with $(length(multiple_blocks)) blocks of sizes $(block_sizes), but largest block ($(max_size)×$(max_size)) exceeds degree 4. Closed-form solutions only exist for degrees ≤ 4."
        end
    end
    
    # Fallback to simple binary split for compatibility
    split = _block_split(mat)
    if !isnothing(split)
        left = mat[1:split, 1:split]
        right = mat[split+1:end, split+1:end]
        vals_left, poly_left, λ = symbolic_eigenvalues(left; var = λ, structure = struct_hint, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
        vals_right, poly_right, _ = symbolic_eigenvalues(right; var = λ, structure = struct_hint, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
        vals = vcat(vals_left, vals_right)
        poly = expand ? Symbolics.expand(poly_left * poly_right) : poly_left * poly_right
        return vals, poly, λ
    end
    if _is_diagonal(mat) || _is_triangular(mat)
        # Shortcut avoids the quartic solver for common structured inputs.
        vals = diag(mat)
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    
    # Before attempting quartic solver, check if matrix is too large
    n = size(mat, 1)
    if n > 4
        error("Cannot compute closed-form eigenvalues for $(n)×$(n) matrix: " *
              "closed-form solutions only exist for degrees ≤ 4 (Abel-Ruffini theorem). " *
              "No exploitable structure (block-diagonal, persymmetric, diagonal, or triangular) was found. " *
              "Consider: (1) rearranging to expose block structure, " *
              "(2) using numeric methods, or (3) computing the characteristic polynomial only.")
    end
    
    poly, coeffs, λ = characteristic_polynomial(A; var = λ)
    vals = symbolic_roots(coeffs; timeout = timeout, max_terms = max_terms)
    # Do not re-expand poly here; characteristic_polynomial already expands as needed
    return vals, poly, λ
end

"""
    symbolic_eigenpairs(A; var = nothing, compute_vectors = true, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)

Return eigenvalues and (optionally) eigenvectors of `A` symbolically.
Yields a vector of `(λ, vectors)` tuples where each `vectors` entry is a basis
for the nullspace of `A - λI`. Set `compute_vectors = false` to avoid the
potentially expensive nullspace computation; in that case, vectors will be
empty arrays.

Diagonal and triangular matrices are detected automatically and bypass the
quartic solver by using their diagonal entries as eigenvalues.

`expand` controls whether to expand the characteristic polynomial for display.
`complexity_threshold` warns when the number of distinct symbolic variables
exceeds this threshold; set to `nothing` to disable.
`timeout` sets maximum computation time in seconds (default: 300, set to `nothing` to disable).
`max_terms` limits expression complexity during simplification (default: 10000).
"""
function symbolic_eigenpairs(A; var = nothing, compute_vectors = true, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)
    λ = isnothing(var) ? _fresh_lambda() : var
    mat = Matrix(A)
    n = size(mat, 1)

    struct_hint = structure === :auto ? _detect_structure(mat) : structure
    if !isnothing(complexity_threshold)
        _check_complexity(mat; threshold = complexity_threshold)
    end
    split = _block_split(mat)
    if !isnothing(split)
        left = mat[1:split, 1:split]
        right = mat[split+1:end, split+1:end]
        pairs_left, poly_left, λ = symbolic_eigenpairs(left; var = λ, compute_vectors = compute_vectors, structure = struct_hint, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
        pairs_right, poly_right, _ = symbolic_eigenpairs(right; var = λ, compute_vectors = compute_vectors, structure = struct_hint, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
        pairs = Vector{Tuple{Any, Vector}}()
        for (val, vecs) in pairs_left
            padded = compute_vectors ? [vcat(vec, zeros(eltype(vec), n - split)) for vec in vecs] : Vector{Any}()
            push!(pairs, (val, padded))
        end
        for (val, vecs) in pairs_right
            padded = compute_vectors ? [vcat(zeros(eltype(vec), split), vec) for vec in vecs] : Vector{Any}()
            push!(pairs, (val, padded))
        end
        poly = expand ? Symbolics.expand(poly_left * poly_right) : poly_left * poly_right
        return pairs, poly, λ
    end

    vals, poly, λ = symbolic_eigenvalues(mat; var = λ, structure = structure, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
    if !compute_vectors
        pairs = [(v, Vector{Any}()) for v in vals]
        return pairs, poly, λ
    end

    I_n = Matrix(I, n, n)
    pairs = Vector{Tuple{Any, Vector}}()
    
    # Note: Eigenvector computation for each eigenvalue is independent and could
    # theoretically be parallelized. However, Symbolics.jl uses task-local storage
    # for hashconsing which is not thread-safe. Attempting to use Threads.@threads
    # here causes crashes due to concurrent access to the hashcons cache.
    # See: https://github.com/JuliaSymbolics/SymbolicUtils.jl/issues/
    for v in vals
        shifted = Symbolics.simplify.(mat .- v .* I_n)
        # Try adjugate-based eigenvectors for small matrices to avoid brittle pivots.
        vecs = _adjugate_vectors(shifted)
        if isempty(vecs)
            vecs = _nullspace(shifted)
        end
        push!(pairs, (v, vecs))
    end
    return pairs, poly, λ
end

"""
    symbolic_diagonalize(A; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)

Attempt to diagonalize `A` symbolically. Returns `(P, D, pairs)` where `P`
collects eigenvectors as columns, `D` is diagonal with eigenvalues, and `pairs`
matches the output of `symbolic_eigenpairs`. An error is thrown if an
insufficient number of linearly independent eigenvectors is found.

`expand` controls whether to expand the characteristic polynomial during intermediate calls.
`complexity_threshold` warns when the number of distinct symbolic variables exceeds this threshold (set to `nothing` to disable).
`timeout` sets maximum computation time in seconds (default: 300, set to `nothing` to disable).
`max_terms` limits expression complexity during simplification (default: 10000).
"""
function symbolic_diagonalize(A; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)
    pairs, _, _ = symbolic_eigenpairs(A; var = var, compute_vectors = true, structure = structure, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
    n = size(A, 1)
    struct_hint = structure === :auto ? _detect_structure(A) : structure
    eigenvectors = Vector{Any}()
    eigenvalues = Any[]
    for (val, vecs) in pairs
        for v in vecs
            push!(eigenvectors, v)
            push!(eigenvalues, val)
        end
    end
    length(eigenvectors) < n && error(_diag_error_msg(struct_hint, length(eigenvectors), n))

    # If we have exactly n eigenvectors and all eigenvalues are distinct,
    # we can skip the expensive RREF computation since distinct eigenvalues
    # guarantee linear independence of their eigenvectors.
    if length(eigenvectors) == n && length(unique(eigenvalues)) == n
        P = reduce(hcat, eigenvectors)
        D = Matrix{eltype(P)}(Diagonal(eigenvalues))
        return P, D, pairs
    end

    # For cases with repeated eigenvalues or more than n eigenvectors,
    # use RREF to check linear independence and select n independent ones.
    B = reduce(hcat, eigenvectors)
    _, pivots = _rref(B)
    length(pivots) < n && error(_diag_error_msg(struct_hint, length(pivots), n))

    # Pick the first n independent eigenvectors (in RREF order) to build P, D.
    selected = pivots[1:n]
    P = reduce(hcat, eigenvectors[selected])
    D = Matrix{eltype(P)}(Diagonal(eigenvalues[selected]))
    return P, D, pairs
end

# ============================================================================
# Structure Detection and Analysis
# Functions to detect matrix properties and exploitable structure
# ============================================================================

function _detect_structure(mat)
    _is_diagonal(mat) && return :diagonal
    _is_triangular(mat) && return :triangular
    _is_hermitian(mat) && return :hermitian
    _is_symmetric(mat) && return :symmetric
    _is_unitary(mat) && return :unitary
    return :none
end

# Count distinct symbolic variables appearing in matrix entries
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

# Warn about complexity based on number of distinct symbolic variables
function _check_complexity(A; threshold = 5, quiet = false)
    n_vars = _count_symbolic_vars(A)
    if !quiet && n_vars > threshold
        @warn "Matrix contains $n_vars symbolic variables (threshold: $threshold). Operations may be slow or timeout with many parameters. Consider structured matrices, partial substitution, or expand=false."
    end
    return n_vars
end

function _diag_error_msg(struct_hint, found, n)
    base = "matrix is not diagonalizable: found $found independent eigenvectors for size $n"
    if struct_hint in (:hermitian, :symmetric, :unitary)
        return "$base (note: $struct_hint matrices should be diagonalizable; check symbolic simplification or assumptions)"
    end
    return base
end

# Detect a block-diagonal split; returns the split index or nothing.
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

# Check if a matrix is persymmetric: Q[i,j] == Q[n+1-j, n+1-i]
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
# Adjugate Matrix Computation
# Helper functions for computing adjugate matrices (used in eigenvector computation)
# ============================================================================

function _adjugate_vectors(M)
    n = size(M, 1)
    n <= 3 || return Vector{Vector{eltype(M)}}()
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
# Matrix Property Checks
# Boolean functions to check matrix properties (diagonal, triangular, symmetric, etc.)
# ============================================================================

function _is_diagonal(mat)
    m, n = size(mat)
    for i in 1:m, j in 1:n
        i == j && continue
        !_issymzero(mat[i, j]) && return false
    end
    return true
end

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

function _is_symmetric(mat)
    m, n = size(mat)
    m == n || return false
    for i in 1:m, j in i+1:n
        !_issymzero(mat[i, j] - mat[j, i]) && return false
    end
    return true
end

function _is_hermitian(mat)
    m, n = size(mat)
    m == n || return false
    for i in 1:m, j in i+1:n
        !_issymzero(mat[i, j] - conj(mat[j, i])) && return false
    end
    return true
end

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
# Special Pattern Detection
# Detection of special matrix patterns with closed-form eigenvalues
# ============================================================================

"""
Detect specific 5×5 tridiagonal patterns with closed-form eigenvalues.

Pattern 1:                    Pattern 2:
    [a  b  0  0  0]              [a  b  0  0  0]
    [b  a  d  0  0]              [b  a  b  0  0]
    [0  d  a  b  0]     OR       [0  b  a  d  0]
    [0  0  b  a  b]              [0  0  d  a  b]
    [0  0  0  b  a]              [0  0  0  b  a]
    
Both patterns have identical closed-form eigenvalues:
    λ₁ = a - √(2b² + d²)
    λ₂ = a - b
    λ₃ = a
    λ₄ = a + b
    λ₅ = a + √(2b² + d²)
    
Returns (a, b, d) if pattern matches, nothing otherwise.
"""
function _detect_special_5x5_tridiagonal(mat)
    size(mat) == (5, 5) || return nothing
    
    # Check that matrix is tridiagonal (zeros beyond ±1 diagonals)
    for i in 1:5, j in 1:5
        if abs(i - j) > 1 && !_issymzero(mat[i, j])
            return nothing
        end
    end
    
    # Extract diagonal and off-diagonal elements
    a = mat[1, 1]
    
    # Check if all diagonal elements are equal
    for i in 2:5
        if !_issymzero(mat[i, i] - a)
            return nothing
        end
    end
    
    # Check if symmetric
    if !_is_symmetric(mat)
        return nothing
    end
    
    # Extract off-diagonal pattern: should be [b, d, b, b] or [b, b, d, b]
    off1 = mat[1, 2]
    off2 = mat[2, 3]
    off3 = mat[3, 4]
    off4 = mat[4, 5]
    
    # Pattern 1: [b, d, b, b] - positions 1, 3, 4 are the same
    if _issymzero(off1 - off3) && _issymzero(off1 - off4)
        b = off1
        d = off2
        return (a, b, d)
    end
    
    # Pattern 2: [b, b, d, b] - positions 1, 2, 4 are the same
    if _issymzero(off1 - off2) && _issymzero(off1 - off4)
        b = off1
        d = off3
        return (a, b, d)
    end
    
    return nothing
end

# ============================================================================
# Utility Functions
# Helper functions for basis vectors, zero checks, and variable management
# ============================================================================

function _std_basis(T, n, idx)
    v = zeros(T, n)
    v[idx] = one(T)
    return v
end

"""
    LinearAlgebra.eigen(A; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)

Compute the eigenvalue decomposition of a symbolic matrix `A`, returning an `Eigen` object
with fields `values` (eigenvalues as a vector) and `vectors` (eigenvectors as columns of a matrix).

This is a familiar interface matching `LinearAlgebra.eigen` for numeric matrices.

# Examples
```julia
using SymbolicDiagonalization, Symbolics, LinearAlgebra

@variables a b
A = [a 0; 0 b]
E = eigen(A)
E.values  # Vector of eigenvalues: [a, b]
E.vectors # Matrix with eigenvectors as columns
```

# Keyword Arguments
- `var`: Variable to use for the characteristic polynomial (default: auto-generated λ)
- `structure`: Hint matrix type (`:auto`, `:hermitian`, `:symmetric`, `:unitary`, `:none`)
- `expand`: Whether to expand the characteristic polynomial (default: `true`)
- `complexity_threshold`: Warn if symbolic variables exceed this count (default: 5, set to `nothing` to disable)
- `timeout`: Maximum computation time in seconds (default: 300, set to `nothing` to disable)
- `max_terms`: Maximum expression size during simplification (default: 10000)

# Throws
- `ErrorException`: If the matrix is not diagonalizable (insufficient independent eigenvectors)
- `ComputationTimeoutError`: If computation exceeds timeout
- `ExpressionComplexityError`: If expressions grow too large
"""
function LinearAlgebra.eigen(A::Matrix{<:Union{Num, Complex{Num}}}; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)
    V, D, _ = symbolic_diagonalize(A; var = var, structure = structure, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
    values = [D[i, i] for i in 1:size(D, 1)]
    return LinearAlgebra.Eigen(values, V)
end

"""
    LinearAlgebra.eigvals(A; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)

Compute only the eigenvalues of a symbolic matrix `A`, returning them as a vector.

This is a familiar interface matching `LinearAlgebra.eigvals` for numeric matrices.

# Examples
```julia
using SymbolicDiagonalization, Symbolics, LinearAlgebra

@variables a b
A = [a 0; 0 b]
λ = eigvals(A)  # Returns [a, b]
```

# Keyword Arguments
- `var`: Variable to use for the characteristic polynomial (default: auto-generated λ)
- `structure`: Hint matrix type (`:auto`, `:hermitian`, `:symmetric`, `:unitary`, `:none`)
- `expand`: Whether to expand the characteristic polynomial (default: `true`)
- `complexity_threshold`: Warn if symbolic variables exceed this count (default: 5, set to `nothing` to disable)
- `timeout`: Maximum computation time in seconds (default: 300, set to `nothing` to disable)
- `max_terms`: Maximum expression size during simplification (default: 10000)

# Throws
- `ComputationTimeoutError`: If computation exceeds timeout
- `ExpressionComplexityError`: If expressions grow too large
"""
function LinearAlgebra.eigvals(A::Matrix{<:Union{Num, Complex{Num}}}; var = nothing, structure = :auto, expand = true, complexity_threshold = 5, timeout = 300, max_terms = 10000)
    vals, _, _ = symbolic_eigenvalues(A; var = var, structure = structure, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
    return vals
end
