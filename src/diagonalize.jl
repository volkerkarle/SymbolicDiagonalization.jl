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

    # Check for hypercube graph Q_n (any size 2^n)
    hypercube_dim = _is_hypercube_graph(mat)
    if !isnothing(hypercube_dim)
        vals = _hypercube_eigenvalues(hypercube_dim)
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    
    # Check for strongly regular graph srg(n,k,λ,μ) (any size n)
    srg_params = _is_strongly_regular_graph(mat)
    if !isnothing(srg_params)
        n, k, λ_param, μ = srg_params
        vals = _strongly_regular_eigenvalues(n, k, λ_param, μ)
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    

        _check_complexity(mat; threshold = complexity_threshold)
    end
    
    # Check for circulant matrix (any size n)
    if _is_circulant(mat)
        vals = _circulant_eigenvalues(mat)
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    
    # Check for block circulant matrix (any size n×n with k×k blocks)
    block_circ_info = _is_block_circulant(mat)
    if !isnothing(block_circ_info)
        n_blocks, block_size, blocks = block_circ_info
        try
            vals = _block_circulant_eigenvalues(mat, n_blocks, block_size, blocks; 
                                               var=λ, timeout=timeout, max_terms=max_terms)
            poly_raw = prod(λ .- vals)
            poly = expand ? Symbolics.expand(poly_raw) : poly_raw
            return vals, poly, λ
        catch e
            # If block eigenvalue computation fails, fall through to other methods
            @debug "Block circulant structure detected but eigenvalue computation failed" exception=e
        end
    end
    
    # Check for Toeplitz tridiagonal (any size n)
    toeplitz_params = _is_toeplitz_tridiagonal(mat)
    if !isnothing(toeplitz_params)
        a, b, c = toeplitz_params
        n = size(mat, 1)
        vals = _toeplitz_tridiagonal_eigenvalues(n, a, b, c)
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    
    # Check for anti-diagonal matrix (any size n)
    if _is_antidiagonal(mat)
        vals = _antidiagonal_eigenvalues(mat)
        if !isnothing(vals)
            poly_raw = prod(λ .- vals)
            poly = expand ? Symbolics.expand(poly_raw) : poly_raw
            return vals, poly, λ
        end
    end
    
    # Check for permutation matrix (any size n)
    if _is_permutation_matrix(mat)
        vals = _compute_permutation_eigenvalues(mat)
        poly_raw = prod(λ .- vals)
        poly = expand ? Symbolics.expand(poly_raw) : poly_raw
        return vals, poly, λ
    end
    
    # Check for Kronecker product A ⊗ B (any size m*n)
    kron_info = _is_kronecker_product(mat)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        try
            vals = _kronecker_eigenvalues(A, B, m, n; var=λ, timeout=timeout, max_terms=max_terms)
            poly_raw = prod(λ .- vals)
            poly = expand ? Symbolics.expand(poly_raw) : poly_raw
            return vals, poly, λ
        catch e
            # If Kronecker eigenvalue computation fails, fall through to other methods
            @debug "Kronecker product structure detected but eigenvalue computation failed" exception=e
        end
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

function _is_antidiagonal(mat)
    m, n = size(mat)
    m == n || return false
    # Check that only anti-diagonal entries (i + j = n + 1) are non-zero
    for i in 1:m, j in 1:n
        if i + j == n + 1
            continue  # Anti-diagonal entry, can be anything
        else
            !_issymzero(mat[i, j]) && return false
        end
    end
    return true
end

# ============================================================================
# Special Pattern Detection
# Detection of special matrix patterns with closed-form eigenvalues
# ============================================================================

"""
    _is_circulant(mat)

Check if a matrix is circulant: each row is a cyclic shift of the first row.
A circulant matrix has the form:
```
[c₀  c₁  c₂  ...  cₙ₋₁]
[cₙ₋₁ c₀  c₁  ...  cₙ₋₂]
[cₙ₋₂ cₙ₋₁ c₀  ...  cₙ₋₃]
[...              ...]
[c₁  c₂  c₃  ...  c₀  ]
```
"""
function _is_circulant(mat)
    n = size(mat, 1)
    n == size(mat, 2) || return false
    n <= 1 && return true  # Trivial case
    
    # Extract the first row
    first_row = mat[1, :]
    
    # Check that each row is a cyclic shift of the first row
    for i in 2:n
        for j in 1:n
            # Row i should have first_row shifted by (i-1) positions to the right
            # mat[i, j] should equal first_row[mod1(j - (i-1), n)]
            expected_idx = mod1(j - (i - 1), n)
            if !_issymzero(mat[i, j] - first_row[expected_idx])
                return false
            end
        end
    end
    
    return true
end

"""
    _circulant_eigenvalues(mat)

Compute eigenvalues of a circulant matrix using the DFT formula.

For a circulant matrix C = circ(c₀, c₁, ..., cₙ₋₁), the eigenvalues are:
    λⱼ = c₀ + c₁ωʲ + c₂ω²ʲ + ... + cₙ₋₁ω⁽ⁿ⁻¹⁾ʲ

where ω = exp(2πi/n) is the n-th primitive root of unity, for j = 0, 1, ..., n-1.

This is essentially evaluating the polynomial p(z) = c₀ + c₁z + ... + cₙ₋₁zⁿ⁻¹
at the n-th roots of unity.
"""
function _circulant_eigenvalues(mat)
    n = size(mat, 1)
    first_row = mat[1, :]
    
    # Compute eigenvalues using DFT formula
    # ω = exp(2πi/n)
    eigenvalues = Vector{Any}(undef, n)
    
    for j in 0:(n-1)
        # λⱼ = Σₖ cₖ ω^(jk) = Σₖ cₖ exp(2πijk/n)
        λ = first_row[1]  # c₀ term
        
        for k in 1:(n-1)
            # ω^(jk) = exp(2πijk/n)
            # We'll use: exp(iθ) = cos(θ) + i*sin(θ)
            θ = 2 * π * j * k / n
            ω_power = cos(θ) + im * sin(θ)
            λ += first_row[k+1] * ω_power
        end
        
        eigenvalues[j+1] = Symbolics.simplify(λ)
    end
    
    return eigenvalues
end

"""
    _is_block_circulant(mat)

Check if a matrix is block circulant: a matrix partitioned into blocks where 
each block row is a cyclic shift of the first block row.

A block circulant matrix has the form:
```
[A₀  A₁  A₂  ...  Aₙ₋₁]
[Aₙ₋₁ A₀  A₁  ...  Aₙ₋₂]
[Aₙ₋₂ Aₙ₋₁ A₀  ...  Aₙ₋₃]
[...              ...]
[A₁  A₂  A₃  ...  A₀  ]
```
where each Aᵢ is a k×k block.

Returns (n_blocks, block_size, blocks) if block circulant, nothing otherwise.
- n_blocks: number of blocks per row/column
- block_size: size of each square block (k×k)
- blocks: vector [A₀, A₁, ..., Aₙ₋₁] of the blocks from the first block row
"""
function _is_block_circulant(mat)
    m, n = size(mat)
    m == n || return nothing
    m <= 1 && return nothing  # Too small for meaningful block structure
    
    # Try different block sizes k where m is divisible by k
    # Start with smallest non-trivial block size
    for k in 2:div(m, 2)
        m % k == 0 || continue
        
        n_blocks = div(m, k)
        n_blocks >= 2 || continue  # Need at least 2 blocks
        
        # Extract blocks from the first block row
        blocks = Vector{Matrix}(undef, n_blocks)
        for block_idx in 1:n_blocks
            col_start = (block_idx - 1) * k + 1
            col_end = block_idx * k
            blocks[block_idx] = mat[1:k, col_start:col_end]
        end
        
        # Check if each block row is a cyclic shift
        is_block_circulant = true
        for block_row in 2:n_blocks
            row_start = (block_row - 1) * k + 1
            row_end = block_row * k
            
            for block_col in 1:n_blocks
                col_start = (block_col - 1) * k + 1
                col_end = block_col * k
                
                # This block should equal blocks[shift_idx] where shift_idx is determined by cyclic shift
                # Block row block_row, column block_col should be blocks[mod1(block_col - (block_row - 1), n_blocks)]
                shift_idx = mod1(block_col - (block_row - 1), n_blocks)
                expected_block = blocks[shift_idx]
                actual_block = mat[row_start:row_end, col_start:col_end]
                
                # Check if blocks match
                if !all(_issymzero(actual_block[i, j] - expected_block[i, j]) 
                       for i in 1:k, j in 1:k)
                    is_block_circulant = false
                    break
                end
            end
            
            is_block_circulant || break
        end
        
        if is_block_circulant
            return (n_blocks, k, blocks)
        end
    end
    
    return nothing
end

"""
    _is_numeric_matrix(mat)

Check if all entries in the matrix are numeric (not symbolic).
"""
function _is_numeric_matrix(mat)
    return all(x -> x isa Number, mat)
end

"""
    _block_circulant_eigenvalues(mat, n_blocks, block_size, blocks)

Compute eigenvalues of a block circulant matrix using block diagonalization.

For a block circulant matrix with blocks [A₀, A₁, ..., Aₙ₋₁], the eigenvalues
are the union of eigenvalues of the n matrices:

    Dⱼ = Σₖ₌₀ⁿ⁻¹ ωʲᵏ Aₖ

where ω = exp(2πi/n) is the n-th primitive root of unity, for j = 0, 1, ..., n-1.

This reduces the n·k × n·k eigenvalue problem to n separate k×k eigenvalue problems.
"""
function _block_circulant_eigenvalues(mat, n_blocks, block_size, blocks; var=nothing, timeout=300, max_terms=10000)
    n = n_blocks
    k = block_size
    
    # Check if the matrix is fully numeric
    is_numeric = _is_numeric_matrix(mat)
    
    all_eigenvalues = Vector{Any}()
    
    for j in 0:(n-1)
        # Compute Dⱼ = Σₖ ωʲᵏ Aₖ
        if is_numeric
            D_j = zeros(Complex{Float64}, k, k)
        else
            D_j = zeros(eltype(mat), k, k)
        end
        
        for block_idx in 0:(n-1)
            # ω^(j*block_idx) = exp(2πi*j*block_idx/n)
            θ = 2 * π * j * block_idx / n
            ω_power = cos(θ) + im * sin(θ)
            
            D_j .+= ω_power .* blocks[block_idx + 1]
        end
        
        # For numeric matrices, use direct eigenvalue computation
        if is_numeric
            vals_j = eigvals(D_j)
            append!(all_eigenvalues, vals_j)
        else
            # For symbolic matrices, simplify and recursively solve
            D_j = Symbolics.simplify.(D_j)
            
            # Recursively solve the k×k eigenvalue problem
            try
                vals_j, _, _ = symbolic_eigenvalues(D_j; var=var, structure=:auto, 
                                                   expand=false, complexity_threshold=nothing,
                                                   timeout=timeout, max_terms=max_terms)
                append!(all_eigenvalues, vals_j)
            catch e
                # If we can't solve the block, this approach won't work
                rethrow(e)
            end
        end
    end
    
    return all_eigenvalues
end

"""
    _is_kronecker_product(mat)

Detect if a matrix is a Kronecker product A ⊗ B.

A Kronecker product A ⊗ B (where A is m×m and B is n×n) produces an (mn)×(mn) matrix
with block structure where each block A[i,j]*B appears at position (i,j).

Returns (A, B, m, n) if pattern detected, nothing otherwise.

Note: This is a heuristic detection that tries different factorizations.
"""
function _is_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Try different factorizations: N = m * n
    # We need at least 2×2 for both factors
    for m in 2:div(N, 2)
        N % m == 0 || continue
        n = div(N, m)
        n >= 2 || continue
        
        # Extract the reference block B from the (1,1) position
        # This should be A[1,1] * B
        B_candidate = mat[1:n, 1:n]
        
        # Check if B_candidate is zero - if so, we need a different reference block
        if all(_issymzero.(B_candidate))
            # Find first non-zero block
            found_nonzero = false
            for i_block in 1:m, j_block in 1:m
                row_start = (i_block - 1) * n + 1
                col_start = (j_block - 1) * n + 1
                block = mat[row_start:row_start+n-1, col_start:col_start+n-1]
                if !all(_issymzero.(block))
                    B_candidate = block
                    found_nonzero = true
                    break
                end
            end
            found_nonzero || continue
        end
        
        # Try to extract A by looking at the pattern of blocks
        A = zeros(eltype(mat), m, m)
        for i_block in 1:m, j_block in 1:m
            row_start = (i_block - 1) * n + 1
            col_start = (j_block - 1) * n + 1
            block = mat[row_start:row_start+n-1, col_start:col_start+n-1]
            
            # Check if block is a scalar multiple of B_candidate
            # Find a non-zero element in B_candidate to determine the scalar
            scalar = nothing
            for i in 1:n, j in 1:n
                if !_issymzero(B_candidate[i, j])
                    # block should equal A[i_block, j_block] * B_candidate
                    scalar = block[i, j] / B_candidate[i, j]
                    break
                end
            end
            
            if isnothing(scalar)
                # B_candidate was zero, so block should be zero
                if !all(_issymzero.(block))
                    @goto next_factorization
                end
                A[i_block, j_block] = 0
            else
                # Verify that block = scalar * B_candidate everywhere
                for i in 1:n, j in 1:n
                    expected = scalar * B_candidate[i, j]
                    if !_issymzero(block[i, j] - expected)
                        @goto next_factorization
                    end
                end
                A[i_block, j_block] = scalar
            end
        end
        
        # Found a valid Kronecker factorization!
        return (A, B_candidate, m, n)
        
        @label next_factorization
    end
    
    return nothing
end

"""
    _kronecker_eigenvalues(A, B, m, n; var=nothing, timeout=nothing, max_terms=nothing)

Compute eigenvalues of Kronecker product A ⊗ B.

The eigenvalues of A ⊗ B are all products λᵢ * μⱼ where λᵢ are eigenvalues of A
and μⱼ are eigenvalues of B.

For m×m matrix A with eigenvalues λ₁,...,λₘ and n×n matrix B with eigenvalues μ₁,...,μₙ,
the (mn)×(mn) Kronecker product A ⊗ B has eigenvalues {λᵢμⱼ : i=1..m, j=1..n}.
"""
function _kronecker_eigenvalues(A, B, m, n; var=nothing, timeout=nothing, max_terms=nothing)
    # Check if A and B are numeric - if so, use direct eigvals
    is_A_numeric = _is_numeric_matrix(A)
    is_B_numeric = _is_numeric_matrix(B)
    
    # Compute eigenvalues of A
    λ_A = if is_A_numeric
        # For numeric matrices, use direct computation
        eigvals(A)
    elseif m <= 4
        vals_A, _, _ = symbolic_eigenvalues(A; var=var, structure=:auto, 
                                          expand=false, complexity_threshold=nothing,
                                          timeout=timeout, max_terms=max_terms)
        vals_A
    else
        # Try recursive detection for A
        vals_A, _, _ = symbolic_eigenvalues(A; var=var, structure=:auto, 
                                          expand=false, complexity_threshold=nothing,
                                          timeout=timeout, max_terms=max_terms)
        vals_A
    end
    
    # Compute eigenvalues of B
    λ_B = if is_B_numeric
        # For numeric matrices, use direct computation
        eigvals(B)
    elseif n <= 4
        vals_B, _, _ = symbolic_eigenvalues(B; var=var, structure=:auto, 
                                          expand=false, complexity_threshold=nothing,
                                          timeout=timeout, max_terms=max_terms)
        vals_B
    else
        # Try recursive detection for B
        vals_B, _, _ = symbolic_eigenvalues(B; var=var, structure=:auto, 
                                          expand=false, complexity_threshold=nothing,
                                          timeout=timeout, max_terms=max_terms)
        vals_B
    end
    
    # Compute all products λᵢ * μⱼ
    eigenvalues = []
    for λ in λ_A, μ in λ_B
        push!(eigenvalues, λ * μ)
    end
    
    return eigenvalues
end

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

"""
    _is_toeplitz_tridiagonal(mat)

Check if a matrix is a symmetric Toeplitz tridiagonal matrix with constant diagonals.
Returns (a, b, c) where:
- a is the main diagonal constant
- b is the subdiagonal constant  
- c is the superdiagonal constant

Returns nothing if the matrix is not of this form or if it's not symmetric.

A symmetric Toeplitz tridiagonal matrix has the form:
```
[a  b  0  0  ...]
[b  a  b  0  ...]
[0  b  a  b  ...]
[...          ...]
[0  0  0  b  a]
```

Note: We only handle the symmetric case (b = c) because the eigenvalue formula
is guaranteed to work for symmetric matrices. For asymmetric tridiagonal matrices,
the formula may not apply if the matrix is defective (non-diagonalizable).
"""
function _is_toeplitz_tridiagonal(mat)
    n = size(mat, 1)
    n == size(mat, 2) || return nothing
    n <= 1 && return nothing  # Need at least 2×2
    
    # Check that matrix is symmetric first
    if !_is_symmetric(mat)
        return nothing
    end
    
    # Check that matrix is tridiagonal (zeros beyond ±1 diagonals)
    for i in 1:n, j in 1:n
        if abs(i - j) > 1 && !_issymzero(mat[i, j])
            return nothing
        end
    end
    
    # Extract diagonal constants
    a = mat[1, 1]
    
    # Check main diagonal is constant
    for i in 2:n
        if !_issymzero(mat[i, i] - a)
            return nothing
        end
    end
    
    # Extract and check subdiagonal (if it exists)
    b = n > 1 ? mat[2, 1] : zero(eltype(mat))
    for i in 3:n
        if !_issymzero(mat[i, i-1] - b)
            return nothing
        end
    end
    
    # For symmetric case, superdiagonal equals subdiagonal
    c = b
    
    return (a, b, c)
end

"""
    _toeplitz_tridiagonal_eigenvalues(n, a, b, c)

Compute eigenvalues of an n×n symmetric Toeplitz tridiagonal matrix with constants (a, b, c).

For a symmetric tridiagonal matrix (b = c), the eigenvalues are:
    λₖ = a + 2b·cos(kπ/(n+1))  for k = 1, 2, ..., n

This formula comes from the theory of orthogonal polynomials and is exact.
"""
function _toeplitz_tridiagonal_eigenvalues(n, a, b, c)
    eigenvalues = Vector{Any}(undef, n)
    
    # For the symmetric case: λₖ = a + 2b·cos(kπ/(n+1))
    # (We verified b = c in the detection function)
    
    for k in 1:n
        θ = k * π / (n + 1)
        λₖ = a + 2 * b * cos(θ)
        eigenvalues[k] = Symbolics.simplify(λₖ)
    end
    
    return eigenvalues
end

"""
    _antidiagonal_eigenvalues(mat)

Compute eigenvalues of an anti-diagonal matrix (non-zero only on anti-diagonal).

An anti-diagonal matrix has the form:
```
[0  0  ... 0  a₁]
[0  0  ... a₂ 0 ]
[     ...       ]
[0  aₙ₋₁ ... 0 0]
[aₙ 0  ... 0  0 ]
```

For a symmetric anti-diagonal matrix, the eigenvalues are ±aₖ (with sign pattern
depending on whether n is even or odd). For asymmetric cases, we use the fact that
the matrix is similar to a diagonal matrix via a permutation.
"""
function _antidiagonal_eigenvalues(mat)
    n = size(mat, 1)
    # Extract anti-diagonal elements
    antidiag = [mat[i, n + 1 - i] for i in 1:n]
    
    # For symmetric anti-diagonal: eigenvalues come in ±pairs
    # For general case: we need to be more careful
    # The eigenvalues depend on the parity structure
    
    # For now, handle the symmetric case
    if _is_symmetric(mat)
        eigenvalues = Vector{Any}(undef, n)
        if n % 2 == 1
            # Odd dimension: one zero eigenvalue, rest come in pairs
            mid = (n + 1) ÷ 2
            eigenvalues[1] = antidiag[mid]
            idx = 2
            for i in 1:mid-1
                eigenvalues[idx] = antidiag[i]
                eigenvalues[idx+1] = -antidiag[i]
                idx += 2
            end
        else
            # Even dimension: all come in ±pairs
            idx = 1
            for i in 1:n÷2
                eigenvalues[idx] = antidiag[i]
                eigenvalues[idx+1] = -antidiag[i]
                idx += 2
            end
        end
        return eigenvalues
    else
        # For non-symmetric, the analysis is more complex
        # Return nothing to fall back to general method
        return nothing
    end
end

# ============================================================================
# Permutation Matrix Pattern Detection and Computation
# ============================================================================

"""
    _is_permutation_matrix(A)

Check if matrix A is a permutation matrix (exactly one 1 in each row and column, rest zeros).

Returns `true` if A is a permutation matrix, `false` otherwise.

# Theory
A permutation matrix P has exactly one 1 in each row and each column, with all other entries being 0.
The eigenvalues of P are roots of unity, determined by its cycle structure.
"""
function _is_permutation_matrix(A)
    n = size(A, 1)
    
    # Helper function to check if a value is 1
    function _is_one(x)
        return _issymzero(x - 1)
    end
    
    # Check each row has exactly one 1 and rest zeros
    for i in 1:n
        count_ones = 0
        for j in 1:n
            val = A[i, j]
            if _is_one(val)
                count_ones += 1
            elseif !_issymzero(val)
                return false  # Found non-zero, non-one entry
            end
        end
        if count_ones != 1
            return false  # Row doesn't have exactly one 1
        end
    end
    
    # Check each column has exactly one 1
    for j in 1:n
        count_ones = 0
        for i in 1:n
            val = A[i, j]
            if _is_one(val)
                count_ones += 1
            end
        end
        if count_ones != 1
            return false  # Column doesn't have exactly one 1
        end
    end
    
    return true
end

"""
    _permutation_to_cycles(A)

Decompose permutation matrix A into its cycle structure.

Returns a vector of cycle lengths (e.g., [3, 2, 2] for cycles of length 3, 2, and 2).

# Example
For a permutation that maps 1→2→3→1, 4→5→4, 6→7→6, returns [3, 2, 2].
"""
function _permutation_to_cycles(A)
    n = size(A, 1)
    visited = falses(n)
    cycles = Int[]
    
    # Helper function to check if a value is 1
    function _is_one(x)
        return _issymzero(x - 1)
    end
    
    for start in 1:n
        if visited[start]
            continue
        end
        
        # Follow the cycle from start
        cycle_length = 0
        current = start
        
        while !visited[current]
            visited[current] = true
            cycle_length += 1
            
            # Find where current maps to
            next = 0
            for j in 1:n
                if _is_one(A[current, j])
                    next = j
                    break
                end
            end
            
            current = next
        end
        
        push!(cycles, cycle_length)
    end
    
    return cycles
end

"""
    _compute_permutation_eigenvalues(A)

Compute eigenvalues of permutation matrix A using cycle decomposition.

# Theory
For a permutation matrix with cycles of lengths [k₁, k₂, ..., kₘ], the eigenvalues are:
- For each cycle of length k: the k-th roots of unity exp(2πij/k) for j = 0, 1, ..., k-1

This provides a closed-form solution for any size permutation matrix.

# Returns
Vector of eigenvalues (may include complex values for cycles of length > 2).
"""
function _compute_permutation_eigenvalues(A)
    cycles = _permutation_to_cycles(A)
    n = size(A, 1)
    eigenvalues = Vector{Any}(undef, n)
    
    idx = 1
    for cycle_length in cycles
        if cycle_length == 1
            # Fixed point: eigenvalue is 1
            eigenvalues[idx] = 1
            idx += 1
        elseif cycle_length == 2
            # 2-cycle: eigenvalues are 1 and -1
            eigenvalues[idx] = 1
            eigenvalues[idx + 1] = -1
            idx += 2
        else
            # k-cycle: eigenvalues are k-th roots of unity
            # exp(2πij/k) for j = 0, 1, ..., k-1
            for j in 0:(cycle_length - 1)
                if j == 0
                    # j=0: eigenvalue is 1
                    eigenvalues[idx] = 1
                else
                    # Use symbolic representation: exp(2*pi*im*j/k)
                    angle = 2 * π * j / cycle_length
                    eigenvalues[idx] = exp(im * angle)
                end
                idx += 1
            end
        end
    end
    
    return eigenvalues
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

# ============================================================================
# Group Theory Patterns - Hypercube Graphs
# ============================================================================

"""
    _is_hypercube_graph(mat)

Check if matrix is the adjacency matrix of a hypercube graph Q_n.

A hypercube graph Q_n has 2^n vertices, where each vertex is an n-bit binary string.
Two vertices are connected if their binary representations differ in exactly one bit.

The adjacency matrix has a recursive structure:
Q_n = [Q_{n-1}  I; I  Q_{n-1}]

Returns `n` (the dimension) if matrix is a hypercube graph, `nothing` otherwise.

# Properties
- Q_n has 2^n vertices
- Each vertex has degree n
- Eigenvalues: λ_k = n - 2k for k = 0, 1, ..., n
- Multiplicity of λ_k is binomial(n, k)
"""
function _is_hypercube_graph(mat)
    N = size(mat, 1)
    
    # Check if size is a power of 2
    n = round(Int, log2(N))
    if 2^n != N
        return nothing
    end
    
    # For Q_0, we have a 1×1 zero matrix
    if n == 0
        return _issymzero(mat[1,1]) ? 0 : nothing
    end
    
    # For Q_1, we have [0 1; 1 0]
    if n == 1
        if _issymzero(mat[1,1]) && _issymzero(mat[2,2]) &&
           _issymzero(mat[1,2] - 1) && _issymzero(mat[2,1] - 1)
            return 1
        else
            return nothing
        end
    end
    
    # For Q_n with n ≥ 2, check recursive structure: [Q_{n-1}  I; I  Q_{n-1}]
    half = div(N, 2)
    
    # Extract blocks
    top_left = mat[1:half, 1:half]
    top_right = mat[1:half, (half+1):N]
    bottom_left = mat[(half+1):N, 1:half]
    bottom_right = mat[(half+1):N, (half+1):N]
    
    # Check if top_right and bottom_left are identity matrices
    for i in 1:half
        for j in 1:half
            if i == j
                if !_issymzero(top_right[i,j] - 1) || !_issymzero(bottom_left[i,j] - 1)
                    return nothing
                end
            else
                if !_issymzero(top_right[i,j]) || !_issymzero(bottom_left[i,j])
                    return nothing
                end
            end
        end
    end
    
    # Check if top_left and bottom_right are identical and recursively hypercubes
    for i in 1:half
        for j in 1:half
            if !_issymzero(top_left[i,j] - bottom_right[i,j])
                return nothing
            end
        end
    end
    
    # Recursively verify top_left is Q_{n-1}
    sub_n = _is_hypercube_graph(top_left)
    if sub_n == n - 1
        return n
    else
        return nothing
    end
end

"""
    _hypercube_eigenvalues(n)

Compute eigenvalues of hypercube graph Q_n.

The eigenvalues of Q_n are:
    λ_k = n - 2k  for k = 0, 1, 2, ..., n

with multiplicity binomial(n, k) for each λ_k.

# Example
For Q_3 (8 vertices, n=3):
- λ_0 = 3 with multiplicity binomial(3,0) = 1
- λ_1 = 1 with multiplicity binomial(3,1) = 3
- λ_2 = -1 with multiplicity binomial(3,2) = 3
- λ_3 = -3 with multiplicity binomial(3,3) = 1
"""
function _hypercube_eigenvalues(n)
    eigenvalues = Vector{Int}(undef, 2^n)
    idx = 1
    
    for k in 0:n
        λ_k = n - 2*k
        mult = binomial(n, k)
        
        for _ in 1:mult
            eigenvalues[idx] = λ_k
            idx += 1
        end
    end
    
    return eigenvalues
end

# ============================================================================
# Group Theory Patterns - Strongly Regular Graphs
# ============================================================================

"""
    _is_strongly_regular_graph(mat)

Check if matrix is the adjacency matrix of a strongly regular graph srg(n,k,λ,μ).

A strongly regular graph with parameters (n, k, λ, μ) is a k-regular graph on n vertices where:
- Every vertex has degree k
- Every pair of adjacent vertices has exactly λ common neighbors
- Every pair of non-adjacent vertices has exactly μ common neighbors

Returns `(n, k, λ, μ)` if the matrix represents an srg, `nothing` otherwise.

# Eigenvalues
An srg(n,k,λ,μ) has exactly 3 distinct eigenvalues:
- k with multiplicity 1
- r = (1/2)[λ - μ + √Δ] with multiplicity f
- s = (1/2)[λ - μ - √Δ] with multiplicity g

where Δ = (λ-μ)² + 4(k-μ) and f + g = n - 1.
"""
function _is_strongly_regular_graph(mat)
    n = size(mat, 1)
    
    if n < 2
        return nothing
    end
    
    # Check it's a 0-1 matrix with zero diagonal
    for i in 1:n
        if !_issymzero(mat[i,i])
            return nothing
        end
        for j in 1:n
            val = mat[i,j]
            if !(_issymzero(val) || _issymzero(val - 1))
                return nothing
            end
        end
    end
    
    # Check symmetry
    if !_is_symmetric(mat)
        return nothing
    end
    
    # Check k-regularity (all row sums equal)
    k = sum(mat[1, :])
    for i in 2:n
        if sum(mat[i, :]) != k
            return nothing
        end
    end
    
    # Compute λ: number of common neighbors for an adjacent pair
    # Find first edge
    edge_i, edge_j = 0, 0
    found_edge = false
    for i in 1:n
        for j in (i+1):n
            if _issymzero(mat[i,j] - 1)
                edge_i, edge_j = i, j
                found_edge = true
                break
            end
        end
        if found_edge
            break
        end
    end
    
    if !found_edge
        # No edges means empty graph, not strongly regular (unless n=1)
        return nothing
    end
    
    # Count common neighbors for this edge
    λ_param = 0
    for t in 1:n
        if t != edge_i && t != edge_j
            if _issymzero(mat[edge_i, t] - 1) && _issymzero(mat[edge_j, t] - 1)
                λ_param += 1
            end
        end
    end
    
    # Compute μ: number of common neighbors for a non-adjacent pair
    # Find first non-edge
    nonedge_i, nonedge_j = 0, 0
    found_nonedge = false
    for i in 1:n
        for j in (i+1):n
            if _issymzero(mat[i,j])
                nonedge_i, nonedge_j = i, j
                found_nonedge = true
                break
            end
        end
        if found_nonedge
            break
        end
    end
    
    if !found_nonedge
        # Complete graph
        μ = k  # In complete graph, all pairs are adjacent
    else
        # Count common neighbors for this non-edge
        μ = 0
        for t in 1:n
            if t != nonedge_i && t != nonedge_j
                if _issymzero(mat[nonedge_i, t] - 1) && _issymzero(mat[nonedge_j, t] - 1)
                    μ += 1
                end
            end
        end
    end
    
    # Verify all edges have λ common neighbors and all non-edges have μ common neighbors
    for i in 1:n
        for j in (i+1):n
            common = 0
            for t in 1:n
                if t != i && t != j
                    if _issymzero(mat[i, t] - 1) && _issymzero(mat[j, t] - 1)
                        common += 1
                    end
                end
            end
            
            if _issymzero(mat[i,j] - 1)  # Adjacent
                if common != λ_param
                    return nothing
                end
            elseif _issymzero(mat[i,j])  # Non-adjacent
                if common != μ
                    return nothing
                end
            end
        end
    end
    
    return (n, k, λ_param, μ)
end

"""
    _strongly_regular_eigenvalues(n, k, λ, μ)

Compute eigenvalues of strongly regular graph srg(n,k,λ,μ).

A strongly regular graph has exactly 3 distinct eigenvalues:
- k with multiplicity 1 (corresponding to the all-ones eigenvector)
- r = (1/2)[λ - μ + √Δ] with multiplicity f
- s = (1/2)[λ - μ - √Δ] with multiplicity g

where:
- Δ = (λ-μ)² + 4(k-μ)
- f = (k(k-λ-1)) / (k-r)
- g = (k(k-λ-1)) / (k-s)
- f + g = n - 1
"""
function _strongly_regular_eigenvalues(n, k, λ, μ)
    # Compute discriminant
    Δ = (λ - μ)^2 + 4*(k - μ)
    
    if Δ < 0
        error("Invalid strongly regular graph parameters: discriminant is negative")
    end
    
    sqrt_Δ = sqrt(Δ)
    
    # Compute the three eigenvalues
    r = (λ - μ + sqrt_Δ) / 2
    s = (λ - μ - sqrt_Δ) / 2
    
    # Compute multiplicities
    # Using formulas: f = (-k - (n-1)*s)/(r-s), g = n - 1 - f
    # But we can also use f + g = n - 1
    if k ≈ r
        f = n - 1
        g = 0
    elseif k ≈ s
        f = 0
        g = n - 1
    else
        f = round(Int, (-k - (n - 1)*s) / (r - s))
        g = n - 1 - f
    end
    
    # Build eigenvalue vector
    eigenvalues = [k]  # Multiplicity 1
    for _ in 1:f
        push!(eigenvalues, r)
    end
    for _ in 1:g
        push!(eigenvalues, s)
    end
    
    return eigenvalues
end
