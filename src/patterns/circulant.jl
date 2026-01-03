# ============================================================================
# Circulant Matrix Pattern Detection and Eigenvalue Computation
# Includes standard circulant and block circulant matrices
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
    # Num is a subtype of Number, so we need to explicitly check for it
    return all(x -> x isa Number && !(x isa Num), mat)
end

"""
    _block_circulant_eigenvalues(mat, n_blocks, block_size, blocks; var=nothing, timeout=300, max_terms=10000)

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
