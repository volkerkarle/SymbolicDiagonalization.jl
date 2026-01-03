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

This detection exploits the structure of rotation Kronecker products:
1. The matrix is orthogonal: R^T R = I
2. The dimension is a power of 2
3. The trace and other invariants match the expected form

Returns a vector of eigenvalues if detected, nothing otherwise.
"""
function _detect_rotation_kronecker_product(mat)
    N = size(mat, 1)
    size(mat, 2) == N || return nothing
    
    # Only for symbolic matrices - numeric matrices go through generic Kronecker detection
    _is_numeric_matrix(mat) && return nothing
    
    # Check if dimension is a power of 2
    N < 4 && return nothing  # Need at least 4×4 for this
    log2_N = log2(N)
    isinteger(log2_N) || return nothing
    k = Int(log2_N)
    k >= 2 || return nothing  # Need at least 2 factors for this to be useful
    
    # Check if orthogonal
    if !_is_symbolic_orthogonal(mat)
        return nothing
    end
    
    # Try to recursively decompose as SO(2) ⊗ (something)
    result = _try_so2_kronecker_decomposition(mat, k)
    return result
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
            if length(args) == 2 && args[2] == 2
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
