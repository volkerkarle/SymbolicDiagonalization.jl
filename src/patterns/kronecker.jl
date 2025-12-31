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
