# ============================================================================
# Kronecker Product Pattern Detection and Eigenvalue Computation
# ============================================================================

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
        
        result = _try_kronecker_factorization(mat, m, n)
        if !isnothing(result)
            return result
        end
    end
    
    return nothing
end

"""
    _try_kronecker_factorization(mat, m, n)

Try to factorize mat as A ⊗ B where A is m×m and B is n×n.
Returns (A, B, m, n) if successful, nothing otherwise.
"""
function _try_kronecker_factorization(mat, m, n)
    # Extract all m×m blocks of size n×n
    blocks = Matrix{Any}(undef, m, m)
    for i_block in 1:m, j_block in 1:m
        row_start = (i_block - 1) * n + 1
        col_start = (j_block - 1) * n + 1
        blocks[i_block, j_block] = mat[row_start:row_start+n-1, col_start:col_start+n-1]
    end
    
    # Find a reference block (first non-zero block)
    ref_i, ref_j = 1, 1
    B_ref = blocks[1, 1]
    if all(_issymzero.(B_ref))
        found = false
        for i in 1:m, j in 1:m
            if !all(_issymzero.(blocks[i, j]))
                ref_i, ref_j = i, j
                B_ref = blocks[i, j]
                found = true
                break
            end
        end
        found || return nothing
    end
    
    # For a Kronecker product A ⊗ B:
    # blocks[i,j] = A[i,j] * B
    # 
    # So: blocks[i,j] / blocks[ref_i, ref_j] = A[i,j] / A[ref_i, ref_j]
    # 
    # We can extract A (up to a scalar) by computing these ratios
    # Then B = blocks[ref_i, ref_j] / A[ref_i, ref_j]
    
    # Extract A matrix by looking at ratios of corresponding elements
    A = Matrix{eltype(mat)}(undef, m, m)
    
    # Find a non-zero position in B_ref to use for scalar extraction
    ref_elem_i, ref_elem_j = 0, 0
    for i in 1:n, j in 1:n
        if !_issymzero(B_ref[i, j])
            ref_elem_i, ref_elem_j = i, j
            break
        end
    end
    ref_elem_i == 0 && return nothing
    
    # Extract A[i,j] = blocks[i,j][ref_elem] / B_ref[ref_elem]
    # But B_ref = A[ref_i, ref_j] * B, so:
    # A[i,j] = blocks[i,j][ref_elem] / B_ref[ref_elem] * A[ref_i, ref_j]
    # 
    # To avoid carrying A[ref_i, ref_j] through, we set A[ref_i, ref_j] = 1
    # and normalize later if needed
    
    for i_block in 1:m, j_block in 1:m
        block = blocks[i_block, j_block]
        if all(_issymzero.(block))
            A[i_block, j_block] = zero(eltype(mat))
        else
            # A[i,j] / A[ref_i, ref_j] = block[ref_elem] / B_ref[ref_elem]
            ratio = Symbolics.simplify(block[ref_elem_i, ref_elem_j] / B_ref[ref_elem_i, ref_elem_j])
            A[i_block, j_block] = ratio
        end
    end
    
    # Verify the factorization: check that all blocks satisfy block[i,j] = A[i,j] * B_ref
    for i_block in 1:m, j_block in 1:m
        block = blocks[i_block, j_block]
        expected = A[i_block, j_block] .* B_ref
        for i in 1:n, j in 1:n
            diff = Symbolics.simplify(block[i, j] - expected[i, j])
            if !_issymzero(diff)
                return nothing
            end
        end
    end
    
    # Success! Now try to simplify A and B_ref
    # A is currently normalized so A[ref_i, ref_j] = 1
    # B_ref = A_original[ref_i, ref_j] * B_original
    #
    # For eigenvalue computation, we want to recover cleaner factors
    # The eigenvalues of A ⊗ B = eigenvalues of (c*A) ⊗ (B/c) for any scalar c
    # So we can use A and B_ref directly - the product of eigenvalues will be correct
    
    # Try to extract cleaner factors by finding common structure
    A_clean, B_clean = _simplify_kronecker_factors(A, B_ref, m, n)
    
    return (A_clean, B_clean, m, n)
end

"""
    _simplify_kronecker_factors(A, B, m, n)

Try to simplify Kronecker factors A and B to have fewer symbolic dependencies.
This is particularly useful when the original product was kron(X, Y) where X and Y
are symbolic matrices with independent variables.
"""
function _simplify_kronecker_factors(A, B, m, n)
    # First, try to detect if this came from kron of two symbolic matrices
    # by analyzing the multiplicative structure of elements
    result = _detect_symbolic_kronecker_factors(A, B, m, n)
    if !isnothing(result)
        return result
    end
    
    # Try to find Lie group structure by rescaling
    # The factorization normalizes A[1,1] = 1, but the original might have been
    # a Lie group with A[1,1] = cos(θ). Try to recover this.
    result_lie = _try_recover_lie_group_factors(A, B, m, n)
    if !isnothing(result_lie)
        return result_lie
    end
    
    # Fallback: try to extract common factor from B
    # Skip this normalization for numeric (non-symbolic) matrices as it can cause type issues
    # and is not needed for numeric eigenvalue computation
    if _is_numeric_matrix(A) && _is_numeric_matrix(B)
        return A, B
    end
    
    common_factor = nothing
    for i in 1:n, j in 1:n
        if !_issymzero(B[i, j])
            common_factor = B[i, j]
            break
        end
    end
    
    if isnothing(common_factor)
        return A, B
    end
    
    # Try to divide all elements of B by this factor
    # Use a matrix type that can hold any result (Float64 or symbolic)
    B_normalized = Matrix{Any}(undef, n, n)
    for i in 1:n, j in 1:n
        if _issymzero(B[i, j])
            B_normalized[i, j] = zero(eltype(B))
        else
            ratio = Symbolics.simplify(B[i, j] / common_factor)
            B_normalized[i, j] = ratio
        end
    end
    
    # Correspondingly scale A
    A_scaled = Symbolics.simplify.(A .* common_factor)
    
    # Check if the normalized version is simpler (has fewer total variables)
    vars_original = _count_symbolic_vars(A) + _count_symbolic_vars(B)
    vars_normalized = _count_symbolic_vars(A_scaled) + _count_symbolic_vars(B_normalized)
    
    if vars_normalized < vars_original
        return A_scaled, B_normalized
    else
        return A, B
    end
end

"""
    _try_recover_lie_group_factors(A, B, m, n)

Try to recover Lie group structure in Kronecker factors.

When we factor kron(R1, R2) where R1, R2 are rotation matrices (SO(2), etc.),
the factorization algorithm normalizes A[1,1] = 1, giving:
- A with elements like sin(θ)/cos(θ) instead of sin(θ), cos(θ)
- B with cos(θ) absorbed into it

This function tries to find a scalar c such that c*A is a Lie group element.
For SO(2): if A[1,1] = 1 and A[2,1] = tan(θ), then c = cos(θ) gives a rotation.

Returns (A_clean, B_clean) if successful, nothing otherwise.
"""
function _try_recover_lie_group_factors(A, B, m, n)
    # For 2×2 factors, try to detect scaled rotation matrices
    if m == 2
        # Check if A looks like a rotation matrix divided by cos(θ)
        # A = [1, -tan(θ); tan(θ), 1] which equals R(θ)/cos(θ)
        # Key property: A[1,1] = A[2,2] = 1 and A[1,2] = -A[2,1]
        
        a11 = A[1, 1]
        a22 = A[2, 2]
        a12 = A[1, 2]
        a21 = A[2, 1]
        
        # Check if diagonal elements are equal (both should be 1 after normalization)
        if _issymzero(Symbolics.simplify(a11 - a22))
            # Check if off-diagonal elements are negatives
            if _issymzero(Symbolics.simplify(a12 + a21))
                # This looks like a scaled rotation!
                # A = [1, -t; t, 1] where t = tan(θ)
                # Original rotation had cos(θ) = 1/√(1 + t²) = 1/√(1 + a21²)
                # But we can get cos(θ) from B more directly
                
                # B[1,1] should be cos(θ) * cos(φ) where cos(φ) is the diagonal of the R2 factor
                # For now, try a simpler approach: extract what multiplies to give B
                
                # The scaling factor is in B. B = original_A[1,1] * original_B
                # For SO(2) ⊗ SO(2): B[1,1] = cos(θ) * cos(φ), B[2,2] = cos(θ) * cos(φ)
                
                # Try to find the common factor in B's diagonal
                b11 = B[1, 1]
                b22 = B[2, 2]
                
                # If B is also a scaled rotation, B[1,1] = B[2,2]
                if _issymzero(Symbolics.simplify(b11 - b22)) && !_issymzero(b11)
                    # Check if B looks like c * rotation for some c
                    # B = [c*cos(φ), -c*sin(φ); c*sin(φ), c*cos(φ)]
                    # We have B[1,1] = c*cos(φ), B[2,1] = c*sin(φ)
                    # So c² = B[1,1]² + B[2,1]² and B[1,1]/c = cos(φ)
                    
                    b21 = B[2, 1]
                    
                    # c = sqrt(B[1,1]² + B[2,1]²)
                    c_squared = Symbolics.simplify(b11^2 + b21^2)
                    
                    # For symbolic expressions, we need to be careful
                    # Try to see if c_squared simplifies to something nice
                    
                    # The scaling factor for A should make A * scale = rotation
                    # For A = [1, -t; t, 1], the correct scale is 1/√(1+t²) = cos(θ)
                    t = a21  # = tan(θ) in the rotation case
                    denom = Symbolics.simplify(1 + t^2)
                    
                    # Compute A_scaled = A / √(1+t²) which should be R(θ)
                    # But we want A_scaled * scale = original rotation
                    # So scale = √(1+t²) / something... this is getting complicated
                    
                    # Simpler: try scale = B[1,1] / (something that makes A[1,1] * scale = cos(θ))
                    # If the original was kron(R1, R2) with R1[1,1] = cos(θ), R2[1,1] = cos(φ)
                    # Then B[1,1] = cos(θ) * cos(φ)
                    # And A[1,1] = 1 = cos(θ)/cos(θ) (normalized)
                    # So to recover: A_clean = A * cos(θ), B_clean = B / cos(θ)
                    # But cos(θ) = B[1,1] / cos(φ) and we don't know cos(φ) separately...
                    
                    # Actually, for eigenvalues, we can use a trick:
                    # The product of eigenvalues of A_clean and B_clean equals product of A and B
                    # So even if we don't factor perfectly, we can just use:
                    # - A scaled by 1/sqrt(det(A)) to normalize det = 1
                    # - B scaled correspondingly
                    
                    det_A = Symbolics.simplify(a11 * a22 - a12 * a21)
                    if !_issymzero(det_A)
                        # Scale A to have det = 1: A_unit = A / sqrt(det_A)
                        # Then B_scaled = B * sqrt(det_A)
                        
                        # For A = [1, -t; t, 1], det = 1 + t² = 1/cos²(θ)
                        # sqrt(det_A) = 1/cos(θ), so A_unit = A * cos(θ) = rotation!
                        
                        # Check if det_A can be written as a perfect square
                        # det_A = 1 + t² where t = tan(θ) = sin(θ)/cos(θ)
                        # det_A = 1 + sin²/cos² = (cos² + sin²)/cos² = 1/cos²
                        
                        # We can't easily compute sqrt(det_A) symbolically in general,
                        # but we can check if A * sqrt(det_A) is a Lie group
                        
                        # For now, just return the original factors
                        # The eigenvalue computation will still work, just less pretty
                    end
                end
            end
        end
    end
    
    return nothing
end

"""
    _detect_symbolic_kronecker_factors(A, B, m, n)

Detect if A and B came from kron(X, Y) where X and Y are symbolic matrices.
If so, extract the cleaner X and Y factors.

When we have kron(X, Y) with symbolic X[i,j] and Y[k,l], the elements are products
like X[i,j]*Y[k,l]. The current extraction gives us:
- A with elements like X[i,j]/X[1,1]  
- B with elements like X[1,1]*Y[k,l]

We want to recover X and Y (or equivalent clean factors).

For nested Kronecker products like X⊗Y⊗Z, elements are products of 3+ terms.
We identify the "outer" factor and group the remaining terms.
"""
function _detect_symbolic_kronecker_factors(A, B, m, n)
    # Check if B[1,1] is a product of array-indexed variables
    B_11 = B[1, 1]
    _issymzero(B_11) && return nothing
    
    # Try to parse B[1,1] as a product of indexed variables
    factors = _extract_product_factors(B_11)
    isnothing(factors) && return nothing
    length(factors) < 2 && return nothing
    
    # For each factor, check if it's an array-indexed variable
    parsed_factors = []
    for f in factors
        result = _parse_array_access(f)
        if !isnothing(result)
            base, idx = result
            push!(parsed_factors, (f, base, idx))
        end
    end
    
    # We need at least 2 array-indexed factors
    length(parsed_factors) < 2 && return nothing
    
    # Try each factor as the "outer" A factor
    # The outer factor should have consistent indexing across the A matrix ratios
    for (outer_f, outer_base, outer_idx) in parsed_factors
        # Try to reconstruct A using this base
        A_check = _try_reconstruct_A_factor(A, outer_base, m)
        if !isnothing(A_check)
            # Now reconstruct B by extracting the outer factor from each element
            B_clean = _try_extract_inner_factor(B, outer_base, n)
            if !isnothing(B_clean)
                return (A_check, B_clean)
            end
        end
    end
    
    return nothing
end

"""
    _try_extract_inner_factor(B, outer_base, n)

Extract the inner factor from B elements by removing the outer_base component.
For B elements like outer_base[i,j] * inner[k,l] or outer_base[i,j] * inner1[k,l] * inner2[k,l],
return a matrix with just the inner parts.
"""
function _try_extract_inner_factor(B, outer_base, n)
    B_clean = Matrix{eltype(B)}(undef, n, n)
    
    for i in 1:n, j in 1:n
        if _issymzero(B[i, j])
            B_clean[i, j] = zero(eltype(B))
            continue
        end
        
        # Get the factors of B[i,j]
        factors = _extract_product_factors(B[i, j])
        if isnothing(factors)
            return nothing
        end
        
        # Find and remove the factor from outer_base
        inner_factors = []
        found_outer = false
        for f in factors
            base, _ = _parse_array_access(f)
            if !isnothing(base) && isequal(base, outer_base) && !found_outer
                found_outer = true
                # Skip this factor (it's the outer one)
            else
                push!(inner_factors, f)
            end
        end
        
        !found_outer && return nothing
        
        # Reconstruct the inner part
        if isempty(inner_factors)
            B_clean[i, j] = one(eltype(B))
        elseif length(inner_factors) == 1
            B_clean[i, j] = inner_factors[1] isa Num ? inner_factors[1] : Num(inner_factors[1])
        else
            # Multiple inner factors - multiply them together
            result = inner_factors[1]
            for k in 2:length(inner_factors)
                result = result * inner_factors[k]
            end
            B_clean[i, j] = result isa Num ? result : Num(result)
        end
    end
    
    return B_clean
end

"""
    _extract_product_factors(expr)

Extract factors from a product expression. Returns a vector of factors or nothing.
"""
function _extract_product_factors(expr)
    expr isa Num || return nothing
    unwrapped = Symbolics.unwrap(expr)
    
    # Check if it's a multiplication
    try
        if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (*)
            return Symbolics.arguments(unwrapped)
        end
    catch
        return nothing
    end
    return nothing
end

"""
    _parse_array_access(expr)

Parse an expression like X[i,j] and return (base_name, (i, j)).
"""
function _parse_array_access(expr)
    try
        # For array-indexed variables like A[1,1], the structure is getindex(A, 1, 1)
        unwrapped = expr isa Num ? Symbolics.unwrap(expr) : expr
        if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === getindex
            args = Symbolics.arguments(unwrapped)
            if length(args) >= 3
                base = args[1]
                indices = (args[2], args[3])
                return (base, indices)
            end
        end
    catch
        return nothing
    end
    return nothing
end

"""
    _try_reconstruct_factors(B, A_base, B_base, m, n)

Try to reconstruct clean B factor assuming B elements are of form A_base[1,1]*B_base[i,j].
"""
function _try_reconstruct_factors(B, A_base, B_base, m, n)
    B_clean = Matrix{eltype(B)}(undef, n, n)
    
    for i in 1:n, j in 1:n
        if _issymzero(B[i, j])
            B_clean[i, j] = zero(eltype(B))
            continue
        end
        
        # Try to extract B_base[i,j] from B[i,j]
        # B[i,j] should be A_base[...]*B_base[i,j]
        factors = _extract_product_factors(B[i, j])
        if isnothing(factors) || length(factors) != 2
            return nothing
        end
        
        # Find which factor is from B_base
        found = false
        for f in factors
            base, idx = _parse_array_access(f)
            if !isnothing(base) && isequal(base, B_base)
                # Reconstruct the clean indexed access
                B_clean[i, j] = f isa Num ? f : Num(f)
                found = true
                break
            end
        end
        !found && return nothing
    end
    
    return (nothing, B_clean)  # A_clean computed separately
end

"""
    _try_reconstruct_A_factor(A, A_base, m)

Try to reconstruct clean A factor from the ratio matrix.
A elements are like A_base[i,j]/A_base[1,1], we want A_base[i,j].

This function validates that the ratio matrix actually corresponds to ratios
of the given base by checking that A[i,j] = A_base[i,j]/A_base[ref] simplifies correctly.
"""
function _try_reconstruct_A_factor(A, A_base, m)
    A_clean = Matrix{eltype(A)}(undef, m, m)
    
    # Find the reference position (where A[i,j] = 1)
    ref_i, ref_j = 0, 0
    for i in 1:m, j in 1:m
        if !_issymzero(A[i, j]) && _issymzero(A[i, j] - 1)
            ref_i, ref_j = i, j
            break
        end
    end
    ref_i == 0 && return nothing
    
    for i in 1:m, j in 1:m
        if _issymzero(A[i, j])
            A_clean[i, j] = zero(eltype(A))
        else
            # Construct A_base[i,j] by directly indexing
            try
                A_clean[i, j] = Num(A_base[i, j])
            catch
                return nothing
            end
            
            # Validate: A[i,j] should equal A_base[i,j] / A_base[ref_i, ref_j]
            # Wrap in Num for proper comparison
            expected_ratio = Num(A_base[i, j]) / Num(A_base[ref_i, ref_j])
            actual_ratio = A[i, j]
            diff = Symbolics.simplify(expected_ratio - actual_ratio)
            if !_issymzero(diff)
                return nothing
            end
        end
    end
    
    return A_clean
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
    eigenvalues = Any[]
    for λ in λ_A, μ in λ_B
        push!(eigenvalues, λ * μ)
    end
    
    return eigenvalues
end

# ============================================================================
# Rotation Matrix Kronecker Product Detection
# ============================================================================

"""
    _detect_rotation_kronecker_product(mat)

Detect if mat is a Kronecker product of SO(2) rotation matrices.

For SO(2) ⊗ SO(2) (4×4 matrix), eigenvalues are e^{i(±θ±φ)}.
For SO(2)^⊗k (2^k × 2^k matrix), eigenvalues are e^{i(±θ₁±θ₂±...±θₖ)}.

Uses the cleaner extraction approach from rotations.jl that extracts (cos, sin)
pairs directly from matrix entries, avoiding sqrt(cos²) expressions.

Returns a vector of eigenvalues if detected, nothing otherwise.
"""
function _detect_rotation_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Only for symbolic matrices - numeric matrices go through generic Kronecker detection
    _is_numeric_matrix(mat) && return nothing
    
    # Check if dimension is a power of 2
    N < 2 && return nothing
    log2_N = log2(N)
    isinteger(log2_N) || return nothing
    k = Int(log2_N)
    
    # For 2×2 matrices, use Lie group detection (handled by _lie_group_eigenvalues)
    k < 2 && return nothing
    
    # Try to extract cos/sin pairs directly from the first column structure
    # This approach handles the scaling correctly and avoids sqrt(cos²) issues
    pairs = try
        _extract_so2_kron_pairs_from_first_column(mat, k)
    catch
        nothing
    end
    
    isnothing(pairs) && return nothing
    
    # Compute clean eigenvalues from the pairs
    return _so2_kron_eigenvalues_from_pairs_v2(pairs)
end

"""
    _extract_so2_kron_cos_sin_pairs_v2(mat, k) -> Vector{Tuple}

Extract (cos(θⱼ), sin(θⱼ)) pairs from a 2^k × 2^k SO(2) Kronecker product.

Uses recursive block decomposition to extract cos/sin directly without sqrt.
For mat = R(θ₁) ⊗ R(θ₂) ⊗ ... ⊗ R(θₖ):
- block_11 = cos(θ₁) * M, block_21 = sin(θ₁) * M
- where M = R(θ₂) ⊗ ... ⊗ R(θₖ)

Extracts c₁, s₁ by dividing block entries by M[1,1] = product of inner cosines.
"""
function _extract_so2_kron_cos_sin_pairs_v2(mat, k)
    N = 2^k
    
    if k == 1
        # Base case: 2×2 should be SO(2)
        c, s = mat[1, 1], mat[2, 1]
        # Verify it's actually SO(2): check structure [c -s; s c]
        if !_issymzero(Symbolics.simplify(mat[1, 2] + s)) ||
           !_issymzero(Symbolics.simplify(mat[2, 2] - c))
            return nothing
        end
        # Verify c² + s² = 1 with trig simplification
        if !_issymzero_trig(Symbolics.simplify(c^2 + s^2 - 1))
            return nothing
        end
        return [(c, s)]
    end
    
    # For k > 1: decompose as R₁ ⊗ M
    half = N ÷ 2
    
    block_11 = mat[1:half, 1:half]
    block_12 = mat[1:half, half+1:N]
    block_21 = mat[half+1:N, 1:half]
    block_22 = mat[half+1:N, half+1:N]
    
    # Verify block structure: block_11 = block_22, block_12 = -block_21
    for i in 1:half, j in 1:half
        if !_issymzero(Symbolics.simplify(block_11[i,j] - block_22[i,j]))
            return nothing
        end
        if !_issymzero(Symbolics.simplify(block_12[i,j] + block_21[i,j]))
            return nothing
        end
    end
    
    # Recursively extract pairs from the inner matrix structure
    # block_11 = c₁ * M, so we need to normalize by c₁ to get M
    # But we don't know c₁ yet...
    #
    # Key insight: recursively decompose block_11 as if it were c₁*M.
    # The inner pairs will be scaled by c₁, but we can recover c₁ from
    # the product of inner cosines.
    #
    # Actually, for block_11 = c₁ * M where M = R₂ ⊗ ... ⊗ Rₖ:
    # - block_11[1,1] = c₁ * M[1,1] = c₁ * c₂ * ... * cₖ
    # - block_21[1,1] = s₁ * M[1,1] = s₁ * c₂ * ... * cₖ
    #
    # If we recursively decompose block_11/c₁ = M, we get the inner pairs.
    # But we need c₁ first, which depends on knowing M[1,1]...
    #
    # Bootstrap: For k=2, we can directly read the pairs from the first column.
    # For k>2, use the recursive structure.
    
    # Get first column entries (standard Kronecker structure)
    # mat[1,1] = c₁ * c₂ * ... * cₖ (all cosines)
    # mat[half+1, 1] = s₁ * c₂ * ... * cₖ (first sin, rest cos)
    #
    # Ratio: mat[half+1,1] / mat[1,1] = s₁/c₁ = tan(θ₁)
    # This gives us tan(θ₁), but we want c₁ and s₁ directly.
    
    # Recursively get inner pairs from block_11 (which is c₁ * M)
    # We treat block_11 as-is and extract its "pairs", which will be scaled.
    inner_pairs = _extract_so2_kron_cos_sin_pairs_v2(block_11, k - 1)
    isnothing(inner_pairs) && return nothing
    
    # The inner_pairs from block_11 = c₁ * M give us (c₁*c₂, c₁*s₂), (c₁*c₃, c₁*s₃)...
    # Wait, that's not right either. Let me think again.
    #
    # For block_11 = c₁ * M where M is a rotation Kronecker product:
    # - block_11 is NOT a rotation Kronecker product (it's scaled by c₁)
    # - So the recursive call will fail unless c₁ = 1
    #
    # Better approach: Use the first column to extract all pairs directly.
    # 
    # For mat = R₁ ⊗ R₂ ⊗ ... ⊗ Rₖ, the first column is:
    # [c₁c₂...cₖ, c₁c₂...sₖ, c₁c₂...cₖ₋₁sₖ₋₁cₖ, ..., s₁s₂...sₖ]
    # 
    # Actually the pattern is: element at index i (0-based) has factor sⱼ if bit j is set.
    # mat[i+1, 1] = ∏ⱼ (if bit j of i is set then sⱼ else cⱼ)
    
    # Direct extraction using ratio method:
    # For each j, compare elements differing only in bit j:
    # mat[2^j + 1, 1] / mat[1, 1] = sⱼ / cⱼ (if all other bits are 0)
    # Then use c² + s² = 1 with the known ratio
    
    pairs = Tuple{Any, Any}[]
    
    for j in 0:(k-1)
        # Index with only bit j set (1-indexed)
        idx_s = (1 << j) + 1  # Element with sⱼ instead of cⱼ at position j
        
        # mat[1,1] = product of all cⱼ
        # mat[idx_s, 1] = sⱼ * (product of all other cⱼ)
        
        all_cos = mat[1, 1]
        one_sin = mat[idx_s, 1]
        
        # Ratio = sⱼ / cⱼ
        # We have: one_sin / all_cos = sⱼ / cⱼ (if other terms cancel)
        # But actually: all_cos = c₀c₁...cⱼ...cₖ₋₁
        #               one_sin = c₀c₁...sⱼ...cₖ₋₁
        # So one_sin / all_cos = sⱼ / cⱼ ✓
        
        # Now: tanⱼ = sⱼ/cⱼ, and cⱼ² + sⱼ² = 1
        # → cⱼ² (1 + tan²ⱼ) = 1 → cⱼ² = 1/(1+tan²ⱼ) = cos²(θⱼ)
        # → cⱼ = ±cos(θⱼ), sⱼ = ±sin(θⱼ)
        #
        # But this still requires sqrt! The key is that the matrix entries
        # already contain the trig functions directly (e.g., cos(α)*cos(β)).
        #
        # If mat = R(α) ⊗ R(β), then:
        # mat[1,1] = cos(α)cos(β)
        # mat[2,1] = cos(α)sin(β)
        # mat[3,1] = sin(α)cos(β)
        # mat[4,1] = sin(α)sin(β)
        #
        # We can factor these: GCD of column 1 entries at indices with bit j=0
        # gives us all the cⱼ products except cⱼ, times something.
        #
        # Alternative: Use the known product structure.
        # For j=0 (last factor): c₀ appears in mat[1,1], mat[3,1] (odd rows ÷2=even)
        # For j=1 (first factor): c₁ appears in mat[1,1], mat[2,1] (rows < half)
        #
        # Division approach: 
        # cⱼ = mat[1,1] / (product of other cᵢ where i≠j)
        # But we don't know the other cᵢ...
        #
        # Recursive factoring:
        # For k=2: mat[1,1] = c₀c₁, mat[2,1] = c₀s₁, mat[3,1] = s₀c₁, mat[4,1] = s₀s₁
        # c₀ = mat[1,1]/c₁, s₀ = mat[3,1]/c₁
        # c₁ = mat[1,1]/c₀, s₁ = mat[2,1]/c₀
        # This is circular unless we can identify the individual trig functions.
        #
        # KEY INSIGHT: If the matrix was constructed from symbolic trig functions,
        # the entries ARE products of cos(θᵢ) and sin(θᵢ). We just need to
        # recognize this structure and extract the angles.
        #
        # Since we've verified the block structure (block_11 = block_22, 
        # block_12 = -block_21), we know this is an SO(2) Kronecker product.
        # The eigenvalues can be computed directly from any valid (c,s) pairs
        # as long as they satisfy c² + s² = 1.
        #
        # SIMPLE APPROACH: Read (cⱼ, sⱼ) directly from the first column structure.
        # For factor j (0-indexed from right), the 2×2 block structure at level j
        # gives us cⱼ*something and sⱼ*something.
        #
        # Level 0: block size 1, elements 1 vs 2, 3 vs 4, etc.
        # Level 1: block size 2, elements 1-2 vs 3-4, 5-6 vs 7-8, etc.
        # etc.
        
        # For level j (0-indexed), the cos part is in indices with bit j = 0,
        # the sin part is in indices with bit j = 1.
        # 
        # We need the ratio sⱼ/cⱼ = mat[idx_with_bit_j_set, 1] / mat[idx_with_bit_j_clear, 1]
        # where both indices have the same bits except bit j.
        #
        # Simplest: use indices 0 and 2^j (base indices)
    end
    
    # Fallback to the original approach that works for k=2
    return _extract_so2_kron_pairs_from_first_column(mat, k)
end

"""
    _extract_so2_kron_pairs_from_first_column(mat, k) -> Vector{Tuple}

Extract (cos, sin) pairs from the first column of an SO(2) Kronecker product.
Works by recursively factoring the first column structure.
"""
function _extract_so2_kron_pairs_from_first_column(mat, k)
    N = 2^k
    col1 = mat[:, 1]
    
    if k == 1
        c, s = col1[1], col1[2]
        return [(c, s)]
    end
    
    half = N ÷ 2
    
    # For mat = R₁ ⊗ M:
    # col1[1:half] = c₁ * M[:,1]
    # col1[half+1:N] = s₁ * M[:,1]
    
    # M is the (k-1)-fold product, M[:,1] is its first column
    # M[1,1] = product of cosines of remaining angles
    
    # Recursively get inner pairs
    # We need M = mat[1:half, 1:half] / c₁, but we don't know c₁
    # However, we can work with the scaled version and recover c₁ at the end
    
    # Extract inner pairs from the structure of col1[1:half]
    # which equals c₁ * M[:,1]
    
    # For k=2: col1 = [c₁c₂, c₁s₂, s₁c₂, s₁s₂]
    # col1[1:2] = [c₁c₂, c₁s₂] = c₁ * [c₂, s₂]
    # col1[3:4] = [s₁c₂, s₁s₂] = s₁ * [c₂, s₂]
    #
    # From [c₁c₂, c₁s₂], we get c₂ = col1[1]/c₁, s₂ = col1[2]/c₁
    # From [s₁c₂, s₁s₂], we get c₂ = col1[3]/s₁, s₂ = col1[4]/s₁
    #
    # Cross-multiply: col1[1]*col1[4] = c₁c₂s₁s₂ = col1[2]*col1[3] ✓
    # This gives us the constraint but not the individual values.
    #
    # Use the block matrix M = block_11 / c₁:
    # M is orthogonal, so M^T M = I
    # block_11^T block_11 = c₁² I (since block_11 = c₁ M)
    # So c₁² = (block_11^T block_11)[1,1] = sum of col1[1:half].^2
    
    block_11 = mat[1:half, 1:half]
    block_21 = mat[half+1:N, 1:half]
    
    # c₁² = ||block_11[:,1]||² (first column of block_11)
    c1_squared = Symbolics.simplify(sum(block_11[i, 1]^2 for i in 1:half))
    s1_squared = Symbolics.simplify(sum(block_21[i, 1]^2 for i in 1:half))
    
    # Try to extract c₁ and s₁ from c₁² and s₁² using pattern matching
    trig_result = _try_extract_trig_from_squared(c1_squared, s1_squared)
    
    if !isnothing(trig_result)
        c1, s1 = trig_result
        # Successfully extracted c₁ = cos(θ₁), s₁ = sin(θ₁)
        # Now get M = block_11 / c₁ and recursively extract its pairs
        M = Symbolics.simplify.(block_11 ./ c1)
        
        inner_pairs = _extract_so2_kron_pairs_from_first_column(M, k - 1)
        if !isnothing(inner_pairs)
            return vcat([(c1, s1)], inner_pairs)
        end
    end
    
    # Fallback: try direct extraction using first column ratios
    # This works when the trig functions appear directly in matrix entries
    
    # For each factor j, compute (cⱼ, sⱼ) from the first column
    # mat[1,1] = c₁c₂...cₖ
    # mat[2^(k-1)+1, 1] = s₁c₂...cₖ (for first factor)
    # etc.
    
    # Try using the product structure directly
    # If mat[1,1] = cos(α)cos(β)... we can try to factor it
    
    # Last resort: compute scaled pairs and hope they work
    # c₁ * M[1,1] = mat[1,1], s₁ * M[1,1] = mat[half+1, 1]
    # If we set "virtual c₁" = mat[1,1] / inner_product and 
    # "virtual s₁" = mat[half+1,1] / inner_product, we get valid eigenvalues
    
    # Get M's first entry from recursive structure
    inner_pairs = _extract_so2_kron_pairs_from_first_column(block_11, k - 1)
    if isnothing(inner_pairs)
        return nothing
    end
    
    # M[1,1] = product of first elements of inner pairs (all the "c" values)
    M_11 = prod(p[1] for p in inner_pairs)
    
    # Now extract c₁ and s₁
    c1 = Symbolics.simplify(mat[1, 1] / M_11)
    s1 = Symbolics.simplify(mat[half + 1, 1] / M_11)
    
    # Verify c₁² + s₁² = 1
    check = Symbolics.simplify(c1^2 + s1^2)
    if !_issymzero_trig(Symbolics.simplify(check - 1))
        return nothing
    end
    
    return vcat([(c1, s1)], inner_pairs)
end

"""
    _so2_kron_eigenvalues_from_pairs_v2(pairs::Vector) -> Vector

Compute eigenvalues from (cos, sin) pairs.
If pairs are directly (cos(θ), sin(θ)), extracts θ and uses the angle-sum formula
for cleaner output: cos(±θ₁±θ₂±...) + i·sin(±θ₁±θ₂±...).
"""
function _so2_kron_eigenvalues_from_pairs_v2(pairs)
    k = length(pairs)
    
    # Try to extract angles directly from the pairs
    # If c = cos(θ) and s = sin(θ), extract θ
    angles = Any[]
    for (c, s) in pairs
        θ = _try_extract_angle_from_cos_sin(c, s)
        if isnothing(θ)
            # Fallback to multiplication approach
            return _so2_kron_eigenvalues_from_pairs_multiply(pairs)
        end
        push!(angles, θ)
    end
    
    # Use the clean angle-sum formula
    eigenvalues = []
    for signs in Iterators.product(fill([-1, 1], k)...)
        angle_sum = sum(s * θ for (s, θ) in zip(signs, angles))
        push!(eigenvalues, cos(angle_sum) + im * sin(angle_sum))
    end
    
    return eigenvalues
end

"""
    _try_extract_angle_from_cos_sin(c, s) -> Union{Num, Nothing}

Try to extract θ from (cos(θ), sin(θ)) pair by pattern matching.
Returns θ if both c and s are trig functions of the same angle, nothing otherwise.
"""
function _try_extract_angle_from_cos_sin(c, s)
    unwrapped_c = Symbolics.unwrap(c)
    unwrapped_s = Symbolics.unwrap(s)
    
    # Check if c = cos(θ) for some θ
    if !Symbolics.iscall(unwrapped_c) || Symbolics.operation(unwrapped_c) !== cos
        return nothing
    end
    θ_c = Num(Symbolics.arguments(unwrapped_c)[1])
    
    # Check if s = sin(θ) for the same θ
    if !Symbolics.iscall(unwrapped_s) || Symbolics.operation(unwrapped_s) !== sin
        return nothing
    end
    θ_s = Num(Symbolics.arguments(unwrapped_s)[1])
    
    # Verify same angle
    if !_issymzero(Symbolics.simplify(θ_c - θ_s))
        return nothing
    end
    
    return θ_c
end

"""
    _so2_kron_eigenvalues_from_pairs_multiply(pairs::Vector) -> Vector

Fallback: compute eigenvalues by multiplying (c ± is) factors.
Used when angle extraction fails.
"""
function _so2_kron_eigenvalues_from_pairs_multiply(pairs)
    k = length(pairs)
    
    eigenvalues = [1 + 0im]
    
    for (c, s) in pairs
        new_eigenvalues = []
        for λ in eigenvalues
            push!(new_eigenvalues, Symbolics.simplify(λ * (c + im*s)))
            push!(new_eigenvalues, Symbolics.simplify(λ * (c - im*s)))
        end
        eigenvalues = new_eigenvalues
    end
    
    return [trig_simplify(λ) for λ in eigenvalues]
end

"""
    _is_symbolic_orthogonal(mat)

Check if mat is orthogonal: mat^T * mat = I (works for symbolic matrices).
Uses trigonometric simplification to handle sin²(θ) + cos²(θ) = 1.
"""
function _is_symbolic_orthogonal(mat)
    n = size(mat, 1)
    product = mat' * mat
    
    for i in 1:n, j in 1:n
        expected = i == j ? 1 : 0
        diff = Symbolics.simplify(product[i, j] - expected)
        if !_issymzero_trig(diff)
            return false
        end
    end
    return true
end

"""
    _try_so2_kronecker_decomposition(mat, k)

Try to decompose a 2^k × 2^k orthogonal matrix as R₁ ⊗ R₂ ⊗ ... ⊗ Rₖ
where each Rᵢ is a 2×2 rotation matrix.

Uses a recursive approach:
1. Extract the top-left 2×2 block structure to find R₁
2. Check if the remaining structure matches R₁ ⊗ (something)
3. Recursively decompose the "something"

Returns a tuple (eigenvalues, cos_sin_pairs) where cos_sin_pairs are the (c, s) 
for each rotation factor, allowing clean eigenvalue computation.
"""
function _try_so2_kronecker_decomposition(mat, k)
    # Use the helper that returns both eigenvalues and angle info
    result = _try_so2_kronecker_decomposition_impl(mat, k)
    if isnothing(result)
        return nothing
    end
    
    eigenvalues, _ = result
    return eigenvalues
end

"""
    _try_so2_kronecker_decomposition_impl(mat, k)

Implementation that returns (eigenvalues, angle_info) where angle_info
contains the (cos, sin) pairs for each rotation factor.

This allows extracting cos/sin directly from matrix entries instead of sqrt.
"""
function _try_so2_kronecker_decomposition_impl(mat, k)
    N = size(mat, 1)
    
    if k == 1
        # Base case: 2×2 matrix, should be SO(2)
        lie_result = _detect_lie_group(mat)
        if !isnothing(lie_result) && lie_result[1] == :SO2
            c, s = lie_result[2]
            eigenvalues = [c + im*s, c - im*s]
            return (eigenvalues, [(c, s)])
        end
        return nothing
    end
    
    # For k > 1: Try to factor as R ⊗ M where R is 2×2 and M is 2^(k-1) × 2^(k-1)
    half_N = N ÷ 2
    
    # In a Kronecker product R ⊗ M:
    # mat = [R[1,1]*M  R[1,2]*M]
    #       [R[2,1]*M  R[2,2]*M]
    # 
    # For R = [c -s; s c]:
    # block_11 = c*M, block_12 = -s*M, block_21 = s*M, block_22 = c*M
    
    block_11 = mat[1:half_N, 1:half_N]
    block_12 = mat[1:half_N, half_N+1:N]
    block_21 = mat[half_N+1:N, 1:half_N]
    block_22 = mat[half_N+1:N, half_N+1:N]
    
    # Check: block_11 should equal block_22 (both are c*M)
    for i in 1:half_N, j in 1:half_N
        diff = Symbolics.simplify(block_11[i, j] - block_22[i, j])
        if !_issymzero(diff)
            return nothing
        end
    end
    
    # Check: block_12 should equal -block_21 (both are ∓s*M)
    for i in 1:half_N, j in 1:half_N
        diff = Symbolics.simplify(block_12[i, j] + block_21[i, j])
        if !_issymzero(diff)
            return nothing
        end
    end
    
    # Key insight: Extract M first (recursively), then get c and s from matrix entries.
    # 
    # For mat = R ⊗ M:
    #   mat[1,1] = c * M[1,1]
    #   mat[half_N+1, 1] = s * M[1,1]
    #
    # So if we can find M[1,1] from the recursive decomposition, we get:
    #   c = mat[1,1] / M[1,1]
    #   s = mat[half_N+1, 1] / M[1,1]
    #
    # This avoids sqrt entirely!
    
    # To recursively decompose M, we need M = block_11 / c, but we don't know c yet.
    # However, we can use the ratio approach: M is proportional to block_11.
    # 
    # For the recursive call, we need to verify M is a valid rotation Kronecker product.
    # We can do this by checking if block_11 / block_11[1,1] has rotation structure.
    
    # Get the (1,1) element of block_11
    b11_ref = block_11[1, 1]
    b21_ref = block_21[1, 1]
    
    if _issymzero(b11_ref) && _issymzero(b21_ref)
        return nothing  # Can't normalize
    end
    
    # Normalize block_11 to get M_normalized = block_11 / block_11[1,1]
    # This makes M_normalized[1,1] = 1
    # But we want M itself, not normalized...
    
    # Alternative approach: Use the structure of the inner product.
    # For mat = R ⊗ M where both are orthogonal:
    #   block_11^T block_11 + block_21^T block_21 = I
    #   (c*M)^T(c*M) + (s*M)^T(s*M) = (c² + s²) M^T M = M^T M = I
    #   So M is orthogonal and block_11^T block_11 = c² I
    #
    # Therefore: c² = (block_11^T block_11)[1,1] = sum of squares of first column
    
    c_squared = Symbolics.simplify(sum(block_11[i, 1]^2 for i in 1:half_N))
    s_squared = Symbolics.simplify(sum(block_21[i, 1]^2 for i in 1:half_N))
    
    # Verify c² + s² = 1 (use trig simplification for same-angle case)
    sum_check = Symbolics.simplify(c_squared + s_squared - 1)
    if !_issymzero_trig(sum_check)
        return nothing
    end
    
    # Now the key: extract c and s directly without sqrt!
    # 
    # For SO(2) rotation, c = cos(θ), s = sin(θ) are the actual trig functions.
    # We have: block_11 = c*M, block_21 = s*M
    # 
    # If we can identify M, we get c = block_11[i,j]/M[i,j] for any non-zero M[i,j].
    #
    # Strategy: Recursively decompose M (scaled by c) to get M's structure,
    # then use that to extract c and s.
    #
    # But block_11 = c*M, so to get M we need to divide by c...
    #
    # Alternative: Express M in terms of its recursive decomposition.
    # If M = R2 ⊗ ... ⊗ Rk, then M[1,1] = cos(θ₂) * cos(θ₃) * ... * cos(θₖ)
    # 
    # And block_11[1,1] = c * M[1,1] = cos(θ₁) * cos(θ₂) * ... * cos(θₖ)
    #     block_21[1,1] = s * M[1,1] = sin(θ₁) * cos(θ₂) * ... * cos(θₖ)
    #
    # So if we first decompose the normalized matrix, we can extract the angles.
    
    # Normalize: M_prop = block_11 / c where we pick c such that M_prop is orthogonal
    # Since block_11 = c*M with M orthogonal, we have block_11/c is orthogonal.
    # We don't know c, but we know c² = c_squared.
    # 
    # Try: if c_squared simplifies to something like cos(θ)², we can extract cos(θ).
    # This happens when the inner M is the identity (k=2 case).
    
    if k == 2
        # Special case: 4×4 = R₁ ⊗ R₂
        # block_11 = c₁ * R₂ = c₁ * [c₂ -s₂; s₂ c₂]
        # block_21 = s₁ * R₂ = s₁ * [c₂ -s₂; s₂ c₂]
        #
        # mat[1,1] = c₁*c₂, mat[2,1] = c₁*s₂
        # mat[3,1] = s₁*c₂, mat[4,1] = s₁*s₂
        #
        # From mat orthogonal: column 1 has norm 1
        # c₁²c₂² + c₁²s₂² + s₁²c₂² + s₁²s₂² = c₁²(c₂²+s₂²) + s₁²(c₂²+s₂²) = c₁² + s₁² = 1 ✓
        #
        # c_squared = c₁²c₂² + c₁²s₂² = c₁² (for this column)
        # Wait, no. c_squared = sum of block_11[:,1]² = (c₁c₂)² + (c₁s₂)² = c₁²
        # s_squared = sum of block_21[:,1]² = (s₁c₂)² + (s₁s₂)² = s₁²
        # 
        # So c_squared = c₁² and s_squared = s₁²!
        # And block_11 = c₁ * R₂
        #
        # We can get R₂ by checking if block_11 / sqrt(c_squared) is a rotation.
        # But sqrt(c₁²) = |c₁| which may differ from c₁ if c₁ < 0.
        #
        # Better: Recognize that block_11[1,1] = c₁*c₂ and block_11[2,1] = c₁*s₂
        # So c₂ = block_11[1,1] / c₁ and s₂ = block_11[2,1] / c₁
        # And c₁ can be extracted if we know c₂.
        #
        # The matrix entries tell us everything:
        # mat[1,1] = c₁c₂, mat[2,1] = c₁s₂, mat[3,1] = s₁c₂, mat[4,1] = s₁s₂
        #
        # Ratios: mat[2,1]/mat[1,1] = s₂/c₂ = tan(θ₂) if c₁ ≠ 0, c₂ ≠ 0
        #         mat[3,1]/mat[1,1] = s₁/c₁ = tan(θ₁) if c₂ ≠ 0
        #         mat[4,1]/mat[1,1] = (s₁/c₁)(s₂/c₂) = tan(θ₁)tan(θ₂)
        #
        # Cross-multiply: mat[1,1]*mat[4,1] = c₁c₂s₁s₂
        #                 mat[2,1]*mat[3,1] = c₁s₂s₁c₂ = same!
        # 
        # We can get c₁,s₁,c₂,s₂ from the four corner elements:
        # Let a = mat[1,1], b = mat[2,1], c = mat[3,1], d = mat[4,1]
        # a = c₁c₂, b = c₁s₂, c = s₁c₂, d = s₁s₂
        #
        # c₁² = a² + b² (from column of R₁⊗R₂ restricted to block_11)
        # c₂² = a² + c² (? no, this mixes c₁ and s₁)
        #
        # Actually from the inner structure:
        # For R₂ = block_11/c₁, we have R₂ = [a/c₁, b₁₂/c₁; b/c₁, a/c₁]
        # where we need R₂ to be SO(2): [c₂ -s₂; s₂ c₂]
        # So: a/c₁ = c₂, b/c₁ = s₂, and we need the second column to match.
        #
        # The constraint c₂² + s₂² = 1 gives: (a/c₁)² + (b/c₁)² = 1
        # → a² + b² = c₁² which matches our earlier c_squared!
        #
        # So c₁ = sqrt(a² + b²) = sqrt(c_squared) -- back to sqrt :-(
        #
        # BUT: We can use the ORIGINAL matrix entries directly!
        # mat = kron(R₁, R₂) where R₁, R₂ are symbolic SO(2) matrices.
        # R₁ = [cos(θ), -sin(θ); sin(θ), cos(θ)]
        # 
        # If we look at mat[1,1] = cos(θ)cos(φ), this is already a product of trig functions.
        # The issue is we computed c_squared = cos²(θ) but Symbolics doesn't simplify
        # sqrt(cos²(θ)) to cos(θ).
        #
        # Solution: Try to FACTOR c_squared to recognize it as a perfect square.
        # If c_squared = cos(θ)², we can recognize this pattern and extract cos(θ).
        
        # Check if block_11 is SO(2) scaled by some factor
        lie_result_11 = _detect_lie_group(block_11)
        if !isnothing(lie_result_11) && lie_result_11[1] == :SO2
            # block_11 detected as SO(2) means it equals some R₂!
            # But wait, block_11 = c₁ * R₂, which is NOT SO(2) unless c₁ = 1.
            # Actually _detect_lie_group checks for the structure, which might 
            # allow scaled rotations. Let's check what it returns.
            c2, s2 = lie_result_11[2]
            # The detection returns the cos/sin from the actual matrix entries,
            # so block_11[1,1] = c2, block_11[2,1] = s2 (up to scaling).
            # We need to check if these are the original trig functions.
        end
        
        # Better approach: extract (c₁, s₁) and (c₂, s₂) directly from ratios
        # that are invariant to the coupling.
        #
        # From: a = c₁c₂, b = c₁s₂, c = s₁c₂, d = s₁s₂
        # We have: a*d = c₁c₂s₁s₂ = b*c (cross-ratio identity)
        #
        # And: a/c = c₁/s₁ = cot(θ₁), b/d = c₁/s₁ = cot(θ₁)
        # Similarly: a/b = c₂/s₂ = cot(θ₂), c/d = c₂/s₂ = cot(θ₂)
        #
        # So we can verify the structure, but to get c₁ directly we still need sqrt.
        #
        # UNLESS: we recognize that mat[1,1] = cos(θ)cos(φ) is already the 
        # symbolic expression we want! 
        #
        # The issue was that we computed c_squared = cos²(θ) and then took sqrt.
        # But actually, in the 4×4 case:
        #   c_squared = mat[1,1]² + mat[2,1]² = (c₁c₂)² + (c₁s₂)² = c₁²
        # 
        # If we look at mat[1,1] directly, it might already be cos(θ)*M[1,1].
        # And for k=2, M = R₂, so M[1,1] = cos(φ).
        # Thus mat[1,1] = cos(θ)cos(φ).
        #
        # The solution: RECOGNIZE cos²(x) + sin²(x) = 1 patterns in c_squared.
        # Or: use the fact that if c_squared = cos²(θ), then:
        #   mat[1,1]² + mat[2,1]² should factor as (cos(θ))² × (something)
        #
        # Actually, let's try: look at the expressions in c_squared and s_squared.
        # If they're of the form cos²(θ) and sin²(θ), we should detect this.
        
        # Try to recognize c_squared as cos(something)^2
        trig_result = _try_extract_trig_from_squared(c_squared, s_squared)
        
        if !isnothing(trig_result)
            c1, s1 = trig_result
            # Successfully extracted c₁ = cos(θ₁), s₁ = sin(θ₁)
            # Now extract R₂ = block_11 / c₁
            M = Symbolics.simplify.(block_11 ./ c1)
            
            # Get eigenvalues of M (should be 2×2 SO(2))
            lie_result_M = _detect_lie_group(M)
            if !isnothing(lie_result_M) && lie_result_M[1] == :SO2
                c2, s2 = lie_result_M[2]
                
                # Build eigenvalues: (c₁ + is₁)(c₂ + is₂), (c₁ + is₁)(c₂ - is₂), etc.
                eig_R1 = [c1 + im*s1, c1 - im*s1]
                eig_R2 = [c2 + im*s2, c2 - im*s2]
                
                eigenvalues = []
                for λ1 in eig_R1, λ2 in eig_R2
                    push!(eigenvalues, Symbolics.simplify(λ1 * λ2))
                end
                
                return (eigenvalues, [(c1, s1), (c2, s2)])
            end
        end
        
        # Fallback: Use the original matrix entries directly
        # mat[1,1], mat[2,1], mat[3,1], mat[4,1] are c₁c₂, c₁s₂, s₁c₂, s₁s₂
        # These are exactly the products we need for eigenvalues!
        #
        # Eigenvalues of R₁⊗R₂ are (c₁±is₁)(c₂±is₂) = c₁c₂-s₁s₂ ± i(c₁s₂+s₁c₂), etc.
        # Which equals cos(θ₁±θ₂) + i·sin(θ₁±θ₂)
        #
        # From mat entries: 
        #   c₁c₂ = mat[1,1], c₁s₂ = mat[2,1], s₁c₂ = mat[3,1], s₁s₂ = mat[4,1]
        #
        # Wait, need to verify the block structure.
        # mat[1,:] = [c₁c₂, -c₁s₂, -s₁c₂, s₁s₂] for kron([c₁,-s₁;s₁,c₁], [c₂,-s₂;s₂,c₂])
        # Let me check: kron(R₁, R₂)[1,1] = R₁[1,1]*R₂[1,1] = c₁c₂ ✓
        # kron(R₁, R₂)[2,1] = R₁[1,1]*R₂[2,1] = c₁*s₂ ✓
        # kron(R₁, R₂)[3,1] = R₁[2,1]*R₂[1,1] = s₁*c₂ ✓  
        # kron(R₁, R₂)[4,1] = R₁[2,1]*R₂[2,1] = s₁*s₂ ✓
        #
        # So from first column: a=c₁c₂, b=c₁s₂, c=s₁c₂, d=s₁s₂
        #
        # Eigenvalues:
        # λ₁ = (c₁+is₁)(c₂+is₂) = c₁c₂-s₁s₂ + i(c₁s₂+s₁c₂) = a-d + i(b+c)
        # λ₂ = (c₁+is₁)(c₂-is₂) = c₁c₂+s₁s₂ + i(-c₁s₂+s₁c₂) = a+d + i(c-b)
        # λ₃ = (c₁-is₁)(c₂+is₂) = c₁c₂+s₁s₂ + i(c₁s₂-s₁c₂) = a+d + i(b-c)
        # λ₄ = (c₁-is₁)(c₂-is₂) = c₁c₂-s₁s₂ + i(-c₁s₂-s₁c₂) = a-d + i(-b-c)
        #
        # These are exactly cos(θ₁±θ₂) + i·sin(θ₁±θ₂)!
        
        a = mat[1, 1]
        b = mat[2, 1]
        c_entry = mat[3, 1]
        d = mat[4, 1]
        
        # Verify the cross-ratio: a*d should equal b*c
        cross_check = Symbolics.simplify(a*d - b*c_entry)
        if !_issymzero(cross_check)
            return nothing
        end
        
        # Build eigenvalues and apply trig simplification for clean output
        raw_eigenvalues = [
            (a - d) + im*(b + c_entry),     # e^{i(θ+φ)}
            (a + d) + im*(c_entry - b),     # e^{i(θ-φ)}
            (a + d) + im*(b - c_entry),     # e^{i(-θ+φ)} = e^{i(φ-θ)}
            (a - d) + im*(-b - c_entry)     # e^{-i(θ+φ)}
        ]
        
        # Apply trig simplification to get clean cos(θ±φ) + i·sin(θ±φ) form
        eigenvalues = [trig_simplify(λ) for λ in raw_eigenvalues]
        
        # We don't have the individual (c,s) pairs, but eigenvalues are correct
        return (eigenvalues, nothing)
    end
    
    # For k > 2: Use recursive decomposition
    # Extract c and s using sqrt (fallback, less clean)
    c = Symbolics.simplify(sqrt(c_squared))
    s = Symbolics.simplify(sqrt(s_squared))
    
    # Eigenvalues of R = [c -s; s c] are c ± is
    eig_R = [c + im*s, c - im*s]
    
    # Extract M = block_11 / c
    if _issymzero(c)
        M = block_21
    else
        M = Symbolics.simplify.(block_11 ./ c)
    end
    
    # Recursively decompose M
    result_M = _try_so2_kronecker_decomposition_impl(M, k - 1)
    if isnothing(result_M)
        return nothing
    end
    
    eig_M, angle_info_M = result_M
    
    eigenvalues = []
    for λ in eig_R, μ in eig_M
        push!(eigenvalues, Symbolics.simplify(λ * μ))
    end
    
    angle_info = isnothing(angle_info_M) ? nothing : vcat([(c, s)], angle_info_M)
    return (eigenvalues, angle_info)
end

"""
    _try_extract_trig_from_squared(c_squared, s_squared)

Try to recognize c_squared and s_squared as cos²(θ) and sin²(θ) for some θ.
Returns (cos(θ), sin(θ)) if successful, nothing otherwise.

This avoids the sqrt(cos²(θ)) problem by pattern matching.
"""
function _try_extract_trig_from_squared(c_squared, s_squared)
    # Check if c_squared is of the form cos(θ)^2
    # by looking at its structure
    
    unwrapped_c = Symbolics.unwrap(c_squared)
    unwrapped_s = Symbolics.unwrap(s_squared)
    
    # Check if c_squared = f^2 for some f
    try
        if Symbolics.iscall(unwrapped_c) && Symbolics.operation(unwrapped_c) === (^)
            args = Symbolics.arguments(unwrapped_c)
            # Note: args[2] may be a symbolic 2, so use isequal with value extraction
            if length(args) == 2 && isequal(Symbolics.value(args[2]), 2)
                base = args[1]
                # Check if base is cos(something)
                if Symbolics.iscall(base)
                    op = Symbolics.operation(base)
                    if op === cos
                        θ = Symbolics.arguments(base)[1]
                        c = Num(base)  # cos(θ)
                        s = sin(Num(θ))  # sin(θ)
                        
                        # Verify s_squared = sin²(θ)
                        expected_s_squared = Symbolics.simplify(s^2)
                        if _issymzero(Symbolics.simplify(s_squared - expected_s_squared))
                            return (c, s)
                        end
                    end
                end
            end
        end
    catch
        # Pattern matching failed
    end
    
    return nothing
end

# ============================================================================
# SO(3) Kronecker Product Detection and Eigenvalue Computation
# ============================================================================

"""
    _detect_so3_kronecker_product(mat)

Detect if mat is a Kronecker product of SO(3) rotation matrices.

For a 9×9 matrix, check if it's R₁ ⊗ R₂ where R₁, R₂ ∈ SO(3).
For a 27×27 matrix, check if it's R₁ ⊗ R₂ ⊗ R₃ where Rᵢ ∈ SO(3).

The eigenvalues of SO(3) ⊗ SO(3) are products of individual SO(3) eigenvalues:
- SO(3) eigenvalues: {1, e^{iθ}, e^{-iθ}}
- SO(3) ⊗ SO(3) eigenvalues: all 9 products

Returns a vector of eigenvalues if detected, nothing otherwise.
"""
function _detect_so3_kronecker_product(mat)
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
    
    # For 9×9, try to extract two SO(3) factors
    if k == 2
        return _detect_so3_kron_2fold(mat)
    elseif k == 3
        return _detect_so3_kron_3fold(mat)
    end
    
    # For higher powers, could implement recursively but skip for now
    return nothing
end

"""
    _detect_so3_kron_2fold(mat)

Detect SO(3) ⊗ SO(3) structure in a 9×9 matrix.

Uses the first column pattern to extract rotation information:
- The 9×9 matrix has block structure with 3×3 blocks
- Block (i,j) = R1[i,j] * R2

Extracts R1 and R2 by analyzing the block structure, then computes
eigenvalues as products of SO(3) eigenvalues.
"""
function _detect_so3_kron_2fold(mat)
    size(mat) == (9, 9) || return nothing
    
    # NOTE: We do NOT check orthogonality of the 9×9 matrix directly, as this
    # can be extremely slow for complex symbolic expressions (e.g., Euler ⊗ Euler).
    # Instead, we try to factor the matrix and check orthogonality of the 3×3 factors,
    # which is much faster.
    
    # Extract 3×3 blocks: blocks[i][j] = R1[i,j] * R2 for R1 ⊗ R2
    blocks = [[mat[3(i-1)+1:3i, 3(j-1)+1:3j] for j in 1:3] for i in 1:3]
    
    # Strategy 1: Check diagonal blocks for SO(3) structure
    # For R1 ⊗ R2, blocks[k][k] = R1[k,k] * R2
    # If R1[k,k] = 1 for some k (axis of rotation), then blocks[k][k] = R2
    # Rz has R1[3,3] = 1, Ry has R1[2,2] = 1, Rx has R1[1,1] = 1
    
    for diag_idx in [3, 2, 1]  # Try in order: Rz, Ry, Rx
        Bkk = blocks[diag_idx][diag_idx]
        if _is_so3_trig(Bkk)
            R2_candidate = Bkk
            
            # Extract R1 from block ratios
            # R1[i,j] = blocks[i][j][p,q] / R2[p,q] for some non-zero R2[p,q]
            # Choose p,q based on which diagonal block we're using
            ref_idx = diag_idx  # Use the same index for the reference element
            R2_ref = R2_candidate[ref_idx, ref_idx]
            
            # If R2_ref is zero or symbolic-zero, try (1,1)
            if _issymzero(R2_ref)
                R2_ref = R2_candidate[1, 1]
                ref_idx = 1
            end
            
            if !_issymzero(R2_ref)
                R1_entries = [Symbolics.simplify(blocks[i][j][ref_idx, ref_idx] / R2_ref) for i in 1:3, j in 1:3]
                R1_candidate = R1_entries
                
                if _is_so3_trig(R1_candidate)
                    # Success! Compute eigenvalues using SO(3) formula
                    return _compute_so3_kron_eigenvalues(R1_candidate, R2_candidate)
                end
            end
        end
    end
    
    # Strategy 2: Try the general Kronecker factorization with SO(3) eigenvalue formula
    # This handles cases like Euler ⊗ Euler where no diagonal block is SO(3)
    kron_info = _try_kronecker_factorization(mat, 3, 3)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        # The extracted factors A and B may be scaled versions of SO(3) matrices.
        # Check if they are orthogonal (up to a scalar) by verifying A^T A is diagonal
        # For a scaled rotation c*R: (c*R)^T (c*R) = c² R^T R = c² I
        
        # Simpler approach: check if A and B are SO(3) using _is_so3_trig
        A_simplified = Symbolics.simplify.(A)
        B_simplified = Symbolics.simplify.(B)
        
        if _is_so3_trig(A_simplified) && _is_so3_trig(B_simplified)
            return _compute_so3_kron_eigenvalues(A_simplified, B_simplified)
        end
        
        # If A and B aren't directly SO(3), they might be scaled.
        # For Euler ⊗ Euler, the factorization normalizes A[ref_i, ref_j] = 1,
        # but the true rotation has no entry = 1.
        # We can still compute eigenvalues using the trace formula.
        
        # Check if A is a scaled SO(3): A^T A should equal c² I for some scalar c²
        AtA = Symbolics.simplify.(A' * A)
        BtB = Symbolics.simplify.(B' * B)
        
        # Check if AtA and BtB are scalar multiples of identity
        A_scale_sq = _check_scalar_identity(AtA)
        B_scale_sq = _check_scalar_identity(BtB)
        
        if !isnothing(A_scale_sq) && !isnothing(B_scale_sq)
            # A = c_A * R_A, B = c_B * R_B where R_A, R_B ∈ SO(3)
            # The Kronecker product A ⊗ B = c_A * c_B * (R_A ⊗ R_B)
            # Eigenvalues: c_A * c_B * (products of R_A, R_B eigenvalues)
            
            # For SO(3), eigenvalues are {1, e^{iθ}, e^{-iθ}} where cos(θ) = (tr(R) - 1) / 2
            # tr(R_A) = tr(A) / c_A, where c_A = sqrt(A_scale_sq)
            # But sqrt is messy symbolically...
            
            # Alternative: compute eigenvalues directly from the unscaled trace
            # tr(A) = c_A * tr(R_A), so tr(R_A) = tr(A) / c_A = tr(A) / sqrt(A_scale_sq)
            # cos(θ_A) = (tr(R_A) - 1) / 2 = (tr(A) / sqrt(A_scale_sq) - 1) / 2
            
            # This still needs sqrt. Better approach: verify via determinant
            # det(A) = c_A³ * det(R_A) = c_A³ * 1 = c_A³
            # So c_A = det(A)^(1/3) = (det(AtA))^(1/6)
            
            # For now, skip this complex case and use Strategy 3
        end
    end
    
    # Strategy 3: Direct trace-based computation for Kronecker products of rotations
    # For K = R₁ ⊗ R₂, we can compute tr(R₁) and tr(R₂) from the block structure
    # without needing to extract clean factors.
    #
    # Key insight: tr(K) = tr(R₁) * tr(R₂)
    # And from the diagonal blocks: sum_i tr(blocks[i][i]) = tr(R₁) * tr(R₂)
    # (since blocks[i][i] = R₁[i,i] * R₂, so tr(blocks[i][i]) = R₁[i,i] * tr(R₂))
    # → sum_i R₁[i,i] * tr(R₂) = tr(R₁) * tr(R₂)
    #
    # If we can verify the Kronecker structure another way, we can compute eigenvalues
    # from the trace information.
    
    result = _try_so3_kron_from_traces(mat, blocks)
    if !isnothing(result)
        return result
    end
    
    return nothing
end

"""
    _check_scalar_identity(M)

Check if M is a scalar multiple of the identity matrix.
Returns the scalar if so, nothing otherwise.

Uses trig-aware zero checking to handle expressions like cos²θ + sin²θ - 1 = 0.
"""
function _check_scalar_identity(M)
    n = size(M, 1)
    size(M, 2) == n || return nothing
    
    # Get the diagonal value (should all be the same)
    diag_val = Symbolics.simplify(M[1, 1])
    
    # Check all diagonal entries are equal (using trig-aware comparison)
    for i in 2:n
        if !_issymzero_trig(M[i, i] - diag_val)
            return nothing
        end
    end
    
    # Check all off-diagonal entries are zero (using trig-aware comparison)
    for i in 1:n, j in 1:n
        if i != j && !_issymzero_trig(M[i, j])
            return nothing
        end
    end
    
    return diag_val
end

"""
    _try_so3_kron_from_traces(mat, blocks)

Try to compute SO(3) ⊗ SO(3) eigenvalues using trace information.
Works even when the factors can't be cleanly extracted.

For K = R₁ ⊗ R₂:
- tr(K) = tr(R₁) * tr(R₂)
- Eigenvalues of R₁: {1, cos(θ₁) ± i*sin(θ₁)} where cos(θ₁) = (tr(R₁) - 1) / 2
- Eigenvalues of K: all 9 products

This approach extracts tr(R₂) from block ratios and tr(R₁) from the sum of diagonal block traces.
"""
function _try_so3_kron_from_traces(mat, blocks)
    # For K = R₁ ⊗ R₂, blocks[i][j] = R₁[i,j] * R₂
    # tr(blocks[i][i]) = R₁[i,i] * tr(R₂)
    # sum_i tr(blocks[i][i]) = tr(R₁) * tr(R₂) = tr(K)
    
    # First, compute tr(K) from the matrix
    tr_K = Symbolics.simplify(sum(mat[i, i] for i in 1:9))
    
    # Strategy A: Try to find a diagonal block that is directly SO(3)
    for k in 1:3
        Bkk = blocks[k][k]
        if _is_so3_trig(Bkk)
            # blocks[k][k] = R₁[k,k] * R₂ = R₂ (since R₁[k,k] = 1 for axis of rotation)
            tr_R2 = Symbolics.simplify(sum(Bkk[i, i] for i in 1:3))
            tr_R1 = Symbolics.simplify(tr_K / tr_R2)
            return _compute_so3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
        end
    end
    
    # Strategy B: Use determinant to extract scale factor
    # For blocks[k][k] = c_k * R₂ where R₂ is SO(3):
    # det(blocks[k][k]) = c_k³ * det(R₂) = c_k³ (since det(R₂) = 1)
    # tr(blocks[k][k]) = c_k * tr(R₂)
    # So: tr(R₂) = tr(blocks[k][k]) / cbrt(det(blocks[k][k]))
    #
    # This avoids the expensive orthogonality check via _issymzero_trig
    for k in 1:3
        Bkk = blocks[k][k]
        det_Bkk = trig_simplify(det(Bkk))
        
        # Skip if determinant is zero
        if _issymzero(det_Bkk)
            continue
        end
        
        # Check if det_Bkk is a perfect cube (like cos(β)³)
        c_k = _try_symbolic_cbrt(det_Bkk)
        if !isnothing(c_k)
            tr_Bkk = Symbolics.simplify(sum(Bkk[i, i] for i in 1:3))
            tr_R2 = trig_simplify(tr_Bkk / c_k)
            tr_R1 = trig_simplify(tr_K / tr_R2)
            return _compute_so3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
        end
    end
    
    # Strategy C: Use orthogonality to extract tr(R₂)² from block traces
    # For K = R₁ ⊗ R₂, blocks[i][j] = R₁[i,j] * R₂
    # So tr(blocks[i][j]) = R₁[i,j] * tr(R₂)
    #
    # Using orthogonality of R₁ (sum of squares of any row = 1):
    # R₁[1,1]² + R₁[1,2]² + R₁[1,3]² = 1
    # Therefore:
    # tr(B₁₁)² + tr(B₁₂)² + tr(B₁₃)² = (R₁[1,1]² + R₁[1,2]² + R₁[1,3]²) * tr(R₂)² = tr(R₂)²
    #
    # This gives us tr(R₂)² directly, then tr(R₂) = sqrt(tr(R₂)²)
    
    # Compute block traces for first row
    tr_row1 = [Symbolics.simplify(sum(blocks[1][j][i, i] for i in 1:3)) for j in 1:3]
    
    # tr(R₂)² = sum of squared traces
    tr_R2_sq = trig_simplify(sum(t^2 for t in tr_row1))
    
    # Try to get tr(R₂) as sqrt(tr_R2_sq)
    # For symbolic expressions, we use sqrt directly (it will simplify if possible)
    tr_R2 = sqrt(tr_R2_sq)
    
    # If tr_R2_sq is a perfect square, sqrt will simplify it
    # Otherwise, tr_R2 will contain sqrt() which is fine for eigenvalue expressions
    
    # Compute tr(R₁) = tr(K) / tr(R₂)
    tr_R1 = trig_simplify(tr_K / tr_R2)
    
    # Return eigenvalues computed from traces
    return _compute_so3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
end

"""
    _check_scaled_orthogonal(M)

Check if M = c * R where R is orthogonal and c is a scalar.
Returns (c², R) if so, where c² is the squared scaling factor and R is the normalized matrix.
Returns nothing if M is not a scaled orthogonal matrix.

Uses trig-aware simplification for expressions like cos²θ + sin²θ = 1.
"""
function _check_scaled_orthogonal(M)
    n = size(M, 1)
    size(M, 2) == n || return nothing
    
    # For M = c * R where R^T R = I:
    # M^T M = c² R^T R = c² I
    # So M^T M should be a scalar multiple of identity
    
    # Compute M^T M with trig-aware simplification
    # The product can generate expressions like cos²θ₁*cos²θ₂ + sin²θ₁ + cos²θ₁*sin²θ₂
    # which simplify to cos²θ₁ (by pulling out cos²θ₁ from the sum with sin²θ₂ + cos²θ₂ = 1)
    MtM = M' * M
    
    # Apply trig simplification to each element
    MtM_simplified = similar(MtM)
    for i in 1:n, j in 1:n
        MtM_simplified[i, j] = trig_simplify(MtM[i, j])
    end
    
    c_sq = _check_scalar_identity(MtM_simplified)
    
    isnothing(c_sq) && return nothing
    _issymzero_trig(c_sq) && return nothing  # Degenerate case
    
    # R = M / c where c = sqrt(c_sq)
    # But sqrt is messy symbolically. Instead, normalize by M[1,1] / R[1,1]
    # We know R is orthogonal, so |R[i,j]| ≤ 1
    
    # Alternative: return c_sq and let caller use it
    # For eigenvalue purposes, we might not need R explicitly
    
    # Try to find a clean sqrt by checking if c_sq is a perfect square
    c = _try_symbolic_sqrt(c_sq)
    if !isnothing(c)
        R = Symbolics.simplify.(M ./ c)
        return (c_sq, R)
    end
    
    # Can't extract clean sqrt, but we know the structure
    # Return c_sq and the original matrix
    return (c_sq, M)
end

"""
    _try_symbolic_sqrt(expr)

Try to compute sqrt(expr) symbolically if expr is a perfect square.
Returns the square root if successful, nothing otherwise.
"""
function _try_symbolic_sqrt(expr)
    unwrapped = Symbolics.unwrap(expr)
    
    # Check if it's already a power
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (^)
        args = Symbolics.arguments(unwrapped)
        if length(args) == 2
            base, exp = args
            # Check if exp is 2
            if isequal(Symbolics.value(exp), 2)
                return Num(base)
            end
        end
    end
    
    # Check if it's a product of squares
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (*)
        args = Symbolics.arguments(unwrapped)
        sqrt_factors = []
        for arg in args
            if Symbolics.iscall(arg) && Symbolics.operation(arg) === (^)
                sub_args = Symbolics.arguments(arg)
                if length(sub_args) == 2 && isequal(Symbolics.value(sub_args[2]), 2)
                    push!(sqrt_factors, sub_args[1])
                else
                    return nothing  # Not a perfect square
                end
            else
                return nothing  # Not a perfect square
            end
        end
        return Num(prod(sqrt_factors))
    end
    
    # Check if it's a simple number (not a symbolic Num)
    if expr isa Real && !(expr isa Num)
        s = sqrt(expr)
        if s^2 == expr
            return s
        end
    end
    
    # Check common trig patterns
    # cos²(θ) → cos(θ), sin²(θ) → sin(θ) (with sign ambiguity)
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (^)
        args = Symbolics.arguments(unwrapped)
        if length(args) == 2 && isequal(Symbolics.value(args[2]), 2)
            base = args[1]
            if Symbolics.iscall(base)
                op = Symbolics.operation(base)
                if op === cos || op === sin
                    # cos²(θ) or sin²(θ) - return abs value (positive sqrt)
                    return Num(base)  # Assuming positive branch
                end
            end
        end
    end
    
    return nothing
end

"""
    _try_symbolic_cbrt(expr)

Try to compute cbrt(expr) (cube root) symbolically if expr is a perfect cube.
Returns the cube root if successful, nothing otherwise.
"""
function _try_symbolic_cbrt(expr)
    unwrapped = Symbolics.unwrap(expr)
    
    # Check if it's already a power x^3
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (^)
        args = Symbolics.arguments(unwrapped)
        if length(args) == 2
            base, exp = args
            # Check if exp is 3
            if isequal(Symbolics.value(exp), 3)
                return Num(base)
            end
        end
    end
    
    # Check if it's a product of cubes
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (*)
        args = Symbolics.arguments(unwrapped)
        cbrt_factors = []
        for arg in args
            if Symbolics.iscall(arg) && Symbolics.operation(arg) === (^)
                sub_args = Symbolics.arguments(arg)
                if length(sub_args) == 2 && isequal(Symbolics.value(sub_args[2]), 3)
                    push!(cbrt_factors, sub_args[1])
                else
                    return nothing  # Not a perfect cube
                end
            else
                return nothing  # Not a perfect cube
            end
        end
        return Num(prod(cbrt_factors))
    end
    
    # Check if it's a simple number (not a symbolic Num)
    if expr isa Real && !(expr isa Num)
        c = cbrt(expr)
        if c^3 ≈ expr
            return c
        end
    end
    
    # Check common trig patterns
    # cos³(θ) → cos(θ), sin³(θ) → sin(θ)
    if Symbolics.iscall(unwrapped) && Symbolics.operation(unwrapped) === (^)
        args = Symbolics.arguments(unwrapped)
        if length(args) == 2 && isequal(Symbolics.value(args[2]), 3)
            base = args[1]
            if Symbolics.iscall(base)
                op = Symbolics.operation(base)
                if op === cos || op === sin
                    return Num(base)
                end
            end
        end
    end
    
    return nothing
end

"""
    _compute_so3_kron_eigenvalues_from_traces(tr_R1, tr_R2)

Compute eigenvalues of SO(3) ⊗ SO(3) from the traces of the factors.

For SO(3) with trace t: cos(θ) = (t - 1) / 2, so eigenvalues are {1, e^{iθ}, e^{-iθ}}.
"""
function _compute_so3_kron_eigenvalues_from_traces(tr_R1, tr_R2)
    # cos(θ₁) = (tr_R1 - 1) / 2
    cos_theta1 = Symbolics.simplify((tr_R1 - 1) / 2)
    # sin(θ₁) = sqrt(1 - cos²(θ₁)) - this introduces sqrt, which is messy
    # Better: use cos(θ₁) and sin(θ₁) = sqrt(1 - cos²(θ₁)) symbolically
    
    cos_theta2 = Symbolics.simplify((tr_R2 - 1) / 2)
    
    # Eigenvalues of R₁: 1, cos(θ₁) + i*sin(θ₁), cos(θ₁) - i*sin(θ₁)
    # Eigenvalues of R₂: 1, cos(θ₂) + i*sin(θ₂), cos(θ₂) - i*sin(θ₂)
    
    # For clean output, use sqrt(1 - cos²) for sin
    sin_theta1_sq = Symbolics.simplify(1 - cos_theta1^2)
    sin_theta2_sq = Symbolics.simplify(1 - cos_theta2^2)
    
    # sin_theta1 = sqrt(sin_theta1_sq) - this may not simplify nicely
    # For Euler angles, cos_theta is a complex expression of multiple angles
    # sqrt(1 - cos²(complex_expr)) won't simplify to sin(complex_expr)
    
    # Alternative: compute the products directly
    # λ₁ * μ₁ where λ = 1, cos(θ₁) ± i*sin(θ₁) and μ = 1, cos(θ₂) ± i*sin(θ₂)
    
    # Let's use the symbolic sqrt and hope it simplifies in some cases
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
    _compute_so3_kron_eigenvalues(R1, R2)

Compute eigenvalues of R1 ⊗ R2 where R1, R2 ∈ SO(3).

Uses the trace formula: SO(3) eigenvalues are {1, e^{iθ}, e^{-iθ}}
where cos(θ) = (tr(R) - 1) / 2.

For Kronecker product, eigenvalues are all 9 products of individual eigenvalues.
"""
function _compute_so3_kron_eigenvalues(R1, R2)
    vals_R1 = _so3_eigenvalues(R1)
    vals_R2 = _so3_eigenvalues(R2)
    
    if isnothing(vals_R1) || isnothing(vals_R2)
        return nothing
    end
    
    eigenvalues = []
    for λ in vals_R1, μ in vals_R2
        product = λ * μ
        if product isa Complex
            # Use light simplification: just Symbolics.simplify, not aggressive_simplify
            # This avoids expensive trig simplification that can time out
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
    _detect_so3_kron_3fold(mat)

Detect SO(3) ⊗ SO(3) ⊗ SO(3) structure in a 27×27 matrix.
"""
function _detect_so3_kron_3fold(mat)
    size(mat) == (27, 27) || return nothing
    
    # First try to factor as 3 ⊗ 9
    kron_info = _try_kronecker_factorization(mat, 3, 9)
    if !isnothing(kron_info)
        A, B, m, n = kron_info
        
        # A should be SO(3), B should be SO(3) ⊗ SO(3)
        A_simplified = Symbolics.simplify.(A)
        if _is_so3_trig(A_simplified)
            vals_A = _so3_eigenvalues(A_simplified)
            
            # Recursively detect B as SO(3) ⊗ SO(3)
            vals_B = _detect_so3_kron_2fold(B)
            
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
        
        # A should be SO(3) ⊗ SO(3), B should be SO(3)
        B_simplified = Symbolics.simplify.(B)
        if _is_so3_trig(B_simplified)
            vals_B = _so3_eigenvalues(B_simplified)
            
            # Recursively detect A as SO(3) ⊗ SO(3)
            vals_A = _detect_so3_kron_2fold(A)
            
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

"""
    _is_so3_trig(A)

Check if A is an SO(3) matrix using trig-aware simplification.
"""
function _is_so3_trig(A)
    size(A) == (3, 3) || return false
    
    # Check orthogonality with trig simplification
    _is_orthogonal(A) || return false
    
    # Check det = 1 with trig simplification
    d = det(A)
    _issymzero_trig(d - 1) || return false
    
    return true
end
