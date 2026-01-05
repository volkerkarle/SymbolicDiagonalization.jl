# ============================================================================
# Tensor Products: Kronecker Product Detection and Eigenvalue Computation
# ============================================================================
#
# Kronecker products A ⊗ B arise from tensor product representations in group
# theory. If V and W are representations of groups G and H, then V ⊗ W is a
# representation of G × H.
#
# Key eigenvalue property:
#   λ(A ⊗ B) = {λᵢ(A) · λⱼ(B)} for all pairs (i,j)
#
# This reduces an (mn) × (mn) eigenvalue problem to solving m × m and n × n
# problems separately, then forming products.
#
# For Lie group Kronecker products (SO(2)⊗k, SU(2)⊗k), see patterns/lie/.
# Those have additional structure giving clean trigonometric forms.
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

