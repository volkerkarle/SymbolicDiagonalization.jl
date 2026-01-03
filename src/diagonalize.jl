# ============================================================================
# Public API - Main Entry Points for Symbolic Diagonalization
# ============================================================================

# VALID_STRUCTURES is defined in SymbolicDiagonalization.jl

"""
    _validate_matrix_input(A, structure) -> Matrix

Validate input matrix for eigenvalue computation.
Checks that A is convertible to a Matrix, is square, and has a valid structure hint.
Returns the converted matrix.
"""
function _validate_matrix_input(A, structure)
    mat = try
        Matrix(A)
    catch e
        msg = e isa Exception && hasproperty(e, :msg) ? e.msg : string(e)
        throw(ArgumentError("Input must be convertible to a Matrix: $msg"))
    end
    
    m, n = size(mat)
    m == n || throw(ArgumentError("Matrix must be square, got $(m)×$(n)"))
    structure in VALID_STRUCTURES || throw(ArgumentError(
        "structure must be one of $(VALID_STRUCTURES), got $(structure)"))
    
    return mat
end

"""
    symbolic_eigenvalues(A; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)

Return `(vals, poly, λ)` where `vals` are symbolic eigenvalues (closed-form up
to degree 4), `poly` is the characteristic polynomial, and `λ` is the symbol
used. This skips eigenvector computation.

`structure` can hint matrix type: `:auto` (detect), `:hermitian`,
`:symmetric`, `:unitary`, `:general`, or `:none`. Structured inputs may bypass the quartic
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
function symbolic_eigenvalues(A; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)
    mat = _validate_matrix_input(A, structure)
    
    λ = isnothing(var) ? _fresh_lambda() : var
    struct_hint = structure === :auto ? _detect_structure(mat) : structure
    
    # Warn about complexity if threshold is set
    if !isnothing(complexity_threshold)
        _check_complexity(mat; threshold = complexity_threshold)
    end

    # Check for Lie group structure FIRST (SO(n), SU(n), Sp(2n), etc.)
    # These have exact closed-form eigenvalue formulas with no floating-point artifacts.
    # Only for symbolic matrices - numeric matrices are handled efficiently by LinearAlgebra
    if eltype(mat) <: Num || eltype(mat) <: Complex{Num}
        lie_vals = _lie_group_eigenvalues(mat)
        if !isnothing(lie_vals)
            return _build_eigenvalue_result(lie_vals, λ, expand)
        end
    end
    
    # Check for Lie algebra representations (spin-j of so(3)/su(2), etc.)
    # These are elements of Lie algebra representations with constrained eigenvalue spectra.
    # Only for symbolic matrices.
    if eltype(mat) <: Num || eltype(mat) <: Complex{Num}
        algebra_vals = _lie_algebra_eigenvalues(mat)
        if !isnothing(algebra_vals)
            return _build_eigenvalue_result(algebra_vals, λ, expand)
        end
    end

    # Check for hypercube graph Q_n (any size 2^n)
    hypercube_dim = _is_hypercube_graph(mat)
    if !isnothing(hypercube_dim)
        vals = _hypercube_eigenvalues(hypercube_dim)
        return _build_eigenvalue_result(vals, λ, expand)
    end
    
    # Check for strongly regular graph srg(n,k,λ,μ) (any size n)
    srg_params = _is_strongly_regular_graph(mat)
    if !isnothing(srg_params)
        n, k, λ_param, μ = srg_params
        vals = _strongly_regular_eigenvalues(n, k, λ_param, μ)
        return _build_eigenvalue_result(vals, λ, expand)
    end
    
    # Check for circulant matrix (any size n)
    if _is_circulant(mat)
        vals = _circulant_eigenvalues(mat)
        return _build_eigenvalue_result(vals, λ, expand)
    end
    
    # Check for block circulant matrix (any size n×n with k×k blocks)
    block_circ_info = _is_block_circulant(mat)
    if !isnothing(block_circ_info)
        n_blocks, block_size, blocks = block_circ_info
        try
            vals = _block_circulant_eigenvalues(mat, n_blocks, block_size, blocks; 
                                               var=λ, timeout=timeout, max_terms=max_terms)
            return _build_eigenvalue_result(vals, λ, expand)
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
        return _build_eigenvalue_result(vals, λ, expand)
    end
    
    # Check for anti-diagonal matrix (any size n)
    if _is_antidiagonal(mat)
        vals = _antidiagonal_eigenvalues(mat)
        if !isnothing(vals)
            return _build_eigenvalue_result(vals, λ, expand)
        end
    end
    
    # Check for permutation matrix (any size n)
    if _is_permutation_matrix(mat)
        vals = _compute_permutation_eigenvalues(mat)
        return _build_eigenvalue_result(vals, λ, expand)
    end
    
    # Check for Kronecker product of rotation matrices (SO(2)^⊗k)
    # This gives clean eigenvalues e^{i(±θ₁±θ₂±...)} without messy factorization
    rotation_kron_result = _detect_rotation_kronecker_product(mat)
    if !isnothing(rotation_kron_result)
        vals = rotation_kron_result
        return _build_eigenvalue_result(vals, λ, expand)
    end
    
    # Check for Kronecker product A ⊗ B (any size m*n)
    kron_info = _is_kronecker_product(mat)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        try
            vals = _kronecker_eigenvalues(A, B, m, n; var=λ, timeout=timeout, max_terms=max_terms)
            # For Kronecker products, avoid expanding the characteristic polynomial
            # as it can be extremely large (product of m*n terms)
            # Just return the factored form: prod(λ - val_i)
            poly = prod(λ .- vals)
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
        return _build_eigenvalue_result(vals, λ, expand)
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
            all_vals = Any[]
            all_polys = Any[]
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
        return _build_eigenvalue_result(vals, λ, expand)
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
    
    # Diagonal shift trick with variable substitution:
    # For a 6-variable matrix [a b c; b d e; c e f], we can:
    # 1. Shift by f*I to get [a-f b c; b d-f e; c e 0]
    # 2. Substitute fresh variables: ã = a-f, d̃ = d-f → [ã b c; b d̃ e; c e 0] (5 vars)
    # 3. Solve for eigenvalues in terms of fresh variables
    # 4. Substitute back: ã → a-f, d̃ → d-f
    # 5. Add shift back: λ_orig = λ_shifted + f
    shift_info = _prepare_diagonal_shift(mat)
    if !isnothing(shift_info)
        shifted_mat, shift_val, back_subs = shift_info
        
        # Compute eigenvalues of the shifted matrix (with fresh variables)
        poly_shifted, coeffs_shifted, λ = characteristic_polynomial(shifted_mat; var = λ)
        vals_shifted = symbolic_roots(coeffs_shifted; timeout = timeout, max_terms = max_terms)
        
        # Apply back-substitution to replace fresh variables with original expressions
        vals_substituted = _apply_back_substitution(vals_shifted, back_subs)
        
        # Add shift back to get original eigenvalues
        vals = [v + shift_val for v in vals_substituted]
        
        # Reconstruct polynomial for original matrix
        poly = prod(λ .- vals)
        if expand
            poly = Symbolics.expand(poly)
        end
        return vals, poly, λ
    end
    
    poly, coeffs, λ = characteristic_polynomial(A; var = λ)
    vals = symbolic_roots(coeffs; timeout = timeout, max_terms = max_terms)
    # Do not re-expand poly here; characteristic_polynomial already expands as needed
    return vals, poly, λ
end

"""
    symbolic_eigenpairs(A; var = nothing, compute_vectors = true, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)

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
function symbolic_eigenpairs(A; var = nothing, compute_vectors = true, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)
    mat = _validate_matrix_input(A, structure)
    n = size(mat, 1)
    
    λ = isnothing(var) ? _fresh_lambda() : var

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
    symbolic_diagonalize(A; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)

Attempt to diagonalize `A` symbolically. Returns `(P, D, pairs)` where `P`
collects eigenvectors as columns, `D` is diagonal with eigenvalues, and `pairs`
matches the output of `symbolic_eigenpairs`. An error is thrown if an
insufficient number of linearly independent eigenvectors is found.

`expand` controls whether to expand the characteristic polynomial during intermediate calls.
`complexity_threshold` warns when the number of distinct symbolic variables exceeds this threshold (set to `nothing` to disable).
`timeout` sets maximum computation time in seconds (default: 300, set to `nothing` to disable).
`max_terms` limits expression complexity during simplification (default: 10000).
"""
function symbolic_diagonalize(A; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)
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
# LinearAlgebra Interface Extensions
# Familiar interfaces matching LinearAlgebra.eigen and LinearAlgebra.eigvals
# ============================================================================

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
- `structure`: Hint matrix type (`:auto`, `:hermitian`, `:symmetric`, `:unitary`, `:general`, `:none`)
- `expand`: Whether to expand the characteristic polynomial (default: `true`)
- `complexity_threshold`: Warn if symbolic variables exceed this count (default: $(DEFAULT_COMPLEXITY_THRESHOLD), set to `nothing` to disable)
- `timeout`: Maximum computation time in seconds (default: $(DEFAULT_TIMEOUT_SECONDS), set to `nothing` to disable)
- `max_terms`: Maximum expression size during simplification (default: $(DEFAULT_MAX_TERMS))

# Throws
- `ErrorException`: If the matrix is not diagonalizable (insufficient independent eigenvectors)
- `ComputationTimeoutError`: If computation exceeds timeout
- `ExpressionComplexityError`: If expressions grow too large
"""
function LinearAlgebra.eigen(A::Matrix{<:Union{Num, Complex{Num}}}; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)
    V, D, _ = symbolic_diagonalize(A; var = var, structure = structure, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
    values = [D[i, i] for i in 1:size(D, 1)]
    return LinearAlgebra.Eigen(values, V)
end

"""
    LinearAlgebra.eigvals(A; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)

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
- `structure`: Hint matrix type (`:auto`, `:hermitian`, `:symmetric`, `:unitary`, `:general`, `:none`)
- `expand`: Whether to expand the characteristic polynomial (default: `true`)
- `complexity_threshold`: Warn if symbolic variables exceed this count (default: $(DEFAULT_COMPLEXITY_THRESHOLD), set to `nothing` to disable)
- `timeout`: Maximum computation time in seconds (default: $(DEFAULT_TIMEOUT_SECONDS), set to `nothing` to disable)
- `max_terms`: Maximum expression size during simplification (default: $(DEFAULT_MAX_TERMS))

# Throws
- `ComputationTimeoutError`: If computation exceeds timeout
- `ExpressionComplexityError`: If expressions grow too large
"""
function LinearAlgebra.eigvals(A::Matrix{<:Union{Num, Complex{Num}}}; var = nothing, structure = :auto, expand = true, complexity_threshold = DEFAULT_COMPLEXITY_THRESHOLD, timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)
    vals, _, _ = symbolic_eigenvalues(A; var = var, structure = structure, expand = expand, complexity_threshold = complexity_threshold, timeout = timeout, max_terms = max_terms)
    return vals
end
