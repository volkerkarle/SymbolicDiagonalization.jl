# ============================================================================
# Finite Group: Cyclic Group Zₙ - Circulant Matrices
# ============================================================================
#
# Circulant matrices commute with the cyclic permutation operator, making them
# invariant under the cyclic group Zₙ. This symmetry forces the characteristic
# polynomial to factor completely (abelian group → 1-dimensional irreps).
#
# The DFT matrix diagonalizes all circulants because it is the character table
# of Zₙ. Eigenvalues are given by evaluating the polynomial p(z) = Σ cⱼzʲ at
# the n-th roots of unity.
#
# Also includes block circulant matrices (Zₙ acting on blocks).
# ============================================================================

"""
    _exact_root_of_unity(n, k)

Return the exact k-th n-th root of unity ω^k where ω = exp(2πi/n).

For small n, returns exact rational/algebraic values:
- n=1: 1
- n=2: {1, -1}
- n=3: {1, (-1+i√3)/2, (-1-i√3)/2}
- n=4: {1, i, -1, -i}
- n=6: {1, (1+i√3)/2, (-1+i√3)/2, -1, (-1-i√3)/2, (1-i√3)/2}

For other n, returns cos(2πk/n) + i·sin(2πk/n) using numeric complex exponential.
"""
function _exact_root_of_unity(n::Int, k::Int)
    # Normalize k to [0, n-1]
    k = mod(k, n)
    
    # Exact values for small n
    if n == 1
        return 1
    elseif n == 2
        return k == 0 ? 1 : -1
    elseif n == 3
        if k == 0
            return 1
        elseif k == 1
            return Complex(-0.5, sqrt(3)/2)
        else  # k == 2
            return Complex(-0.5, -sqrt(3)/2)
        end
    elseif n == 4
        return [1, im, -1, -im][k+1]
    elseif n == 6
        sqrt3_2 = sqrt(3)/2
        roots = [
            Complex(1.0, 0.0),           # k=0
            Complex(0.5, sqrt3_2),       # k=1
            Complex(-0.5, sqrt3_2),      # k=2
            Complex(-1.0, 0.0),          # k=3
            Complex(-0.5, -sqrt3_2),     # k=4
            Complex(0.5, -sqrt3_2)       # k=5
        ]
        return roots[k+1]
    elseif n == 8
        sqrt2_2 = sqrt(2)/2
        roots = [
            Complex(1.0, 0.0),           # k=0
            Complex(sqrt2_2, sqrt2_2),   # k=1
            Complex(0.0, 1.0),           # k=2
            Complex(-sqrt2_2, sqrt2_2),  # k=3
            Complex(-1.0, 0.0),          # k=4
            Complex(-sqrt2_2, -sqrt2_2), # k=5
            Complex(0.0, -1.0),          # k=6
            Complex(sqrt2_2, -sqrt2_2)   # k=7
        ]
        return roots[k+1]
    else
        # General case: use numeric complex exponential
        θ = 2 * π * k / n
        return Complex(cos(θ), sin(θ))
    end
end

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

For symbolic matrices, uses exact roots of unity when possible to avoid
floating-point contamination in the output.
"""
function _circulant_eigenvalues(mat)
    n = size(mat, 1)
    first_row = mat[1, :]
    
    # Check if matrix has symbolic entries
    has_symbolic = any(x -> x isa Num, first_row)
    
    # Compute eigenvalues using DFT formula
    eigenvalues = Vector{Any}(undef, n)
    
    for j in 0:(n-1)
        # λⱼ = Σₖ cₖ ω^(jk)
        λ = first_row[1]  # c₀ term (ω^0 = 1)
        
        for k in 1:(n-1)
            if has_symbolic
                # Use exact roots of unity for symbolic matrices
                ω_power = _exact_root_of_unity(n, j * k)
            else
                # Use floating point for numeric matrices
                θ = 2 * π * j * k / n
                ω_power = cos(θ) + im * sin(θ)
            end
            λ = λ + first_row[k+1] * ω_power
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

# ============================================================================
# BCCB: Block Circulant with Circulant Blocks (Zₙ × Zₘ)
# ============================================================================
#
# A BCCB matrix is a block circulant matrix where each block is itself
# circulant. This corresponds to the direct product group Zₙ × Zₘ.
#
# Diagonalized by the 2D DFT (Kronecker product of two DFT matrices).
# Common in: image processing, 2D convolution, PDEs on rectangular grids.
#
# For an (nm) × (nm) BCCB with n blocks of size m×m:
#   λⱼₖ = Σₚ₌₀ⁿ⁻¹ Σᵧ₌₀ᵐ⁻¹ c[p,q] · ωₙʲᵖ · ωₘᵏᵧ
# where ωₙ = e^{2πi/n}, ωₘ = e^{2πi/m}
# ============================================================================

"""
    _is_bccb(mat)

Check if a matrix is BCCB (Block Circulant with Circulant Blocks).

A BCCB matrix is:
1. Block circulant (blocks arranged in circulant pattern)
2. Each block is itself circulant

This corresponds to invariance under the group Zₙ × Zₘ.

Returns (n, m, first_rows) if BCCB, where:
- n: number of blocks
- m: size of each block  
- first_rows: n×m matrix where row p is the first row of block p

Returns `nothing` if not BCCB.
"""
function _is_bccb(mat)
    # First check if block circulant
    result = _is_block_circulant(mat)
    result === nothing && return nothing
    
    n_blocks, block_size, blocks = result
    
    # Now check if each block is circulant
    first_rows = Matrix{Any}(undef, n_blocks, block_size)
    
    for p in 1:n_blocks
        block = blocks[p]
        
        # Check if this block is circulant
        m = block_size
        block_first_row = block[1, :]
        
        for i in 2:m
            for j in 1:m
                expected_idx = mod1(j - (i - 1), m)
                if !_issymzero(block[i, j] - block_first_row[expected_idx])
                    return nothing  # Block is not circulant
                end
            end
        end
        
        first_rows[p, :] = block_first_row
    end
    
    return (n_blocks, block_size, first_rows)
end

"""
    _bccb_eigenvalues(n, m, first_rows)

Compute eigenvalues of a BCCB matrix using 2D DFT.

For a BCCB matrix with n blocks of size m, the nm eigenvalues are:

    λⱼₖ = Σₚ₌₀ⁿ⁻¹ Σᵧ₌₀ᵐ⁻¹ c[p,q] · ωₙʲᵖ · ωₘᵏᵧ

where:
- c[p,q] is entry q of the first row of block p
- ωₙ = exp(2πi/n), ωₘ = exp(2πi/m)
- j = 0, ..., n-1 and k = 0, ..., m-1

This is the 2D DFT of the "generating matrix" c[p,q].
"""
function _bccb_eigenvalues(n, m, first_rows)
    eigenvalues = Vector{Any}(undef, n * m)
    idx = 1
    
    for j in 0:(n-1)
        for k in 0:(m-1)
            # λⱼₖ = Σₚ Σᵧ c[p,q] · ωₙʲᵖ · ωₘᵏᵧ
            λ = zero(first_rows[1, 1])
            
            for p in 0:(n-1)
                for q in 0:(m-1)
                    c_pq = first_rows[p + 1, q + 1]
                    
                    # ωₙʲᵖ · ωₘᵏᵧ = exp(2πi·jp/n) · exp(2πi·kq/m)
                    θ1 = 2 * π * j * p / n
                    θ2 = 2 * π * k * q / m
                    
                    ω_power = (cos(θ1) + im * sin(θ1)) * (cos(θ2) + im * sin(θ2))
                    λ = λ + c_pq * ω_power
                end
            end
            
            eigenvalues[idx] = Symbolics.simplify(λ)
            idx += 1
        end
    end
    
    return eigenvalues
end

# ============================================================================
# Circulant Eigenvectors (DFT basis)
# ============================================================================
#
# The eigenvectors of any n×n circulant matrix are the columns of the DFT matrix.
# This is because circulant matrices commute with the cyclic shift operator,
# and the DFT diagonalizes all such operators simultaneously.
#
# The k-th eigenvector (0-indexed) is:
#   vₖ = [1, ωᵏ, ω²ᵏ, ..., ω⁽ⁿ⁻¹⁾ᵏ]ᵀ / √n  (normalized)
#
# where ω = exp(2πi/n) is the primitive n-th root of unity.
# ============================================================================

"""
    _dft_column(n, k; normalized=false)

Return the k-th column of the n×n DFT matrix (0-indexed).

This is the eigenvector common to all n×n circulant matrices:
    vₖ = [1, ωᵏ, ω²ᵏ, ..., ω⁽ⁿ⁻¹⁾ᵏ]ᵀ

where ω = exp(2πi/n).

If `normalized=true`, divides by √n for unit norm.

For small n, uses exact algebraic roots of unity to avoid floating-point errors.
"""
function _dft_column(n::Int, k::Int; normalized::Bool=false)
    # Normalize k to [0, n-1]
    k = mod(k, n)
    
    # Check if any symbolic entries would benefit from exact roots
    # For now, always use exact roots for cleaner output
    vec = Vector{Any}(undef, n)
    
    for j in 0:(n-1)
        # Entry j is ω^(jk)
        vec[j+1] = _exact_root_of_unity(n, j * k)
    end
    
    if normalized
        scale = 1 / sqrt(Symbolics.Num(n))
        vec = [v * scale for v in vec]
    end
    
    return vec
end

"""
    _circulant_eigenvectors(n; normalized=false)

Return all n eigenvectors of an n×n circulant matrix.

The eigenvectors are the columns of the DFT matrix, which diagonalize
all circulant matrices. They are returned in the same order as the
eigenvalues from `_circulant_eigenvalues`.

Returns a vector of vectors: [v₀, v₁, ..., vₙ₋₁] where vₖ corresponds
to eigenvalue λₖ.

If `normalized=true`, eigenvectors have unit norm (divided by √n).
"""
function _circulant_eigenvectors(n::Int; normalized::Bool=false)
    return [_dft_column(n, k; normalized=normalized) for k in 0:(n-1)]
end

"""
    _circulant_eigenpairs(mat; normalized=false)

Compute eigenpairs (eigenvalue, eigenvector) for a circulant matrix.

For a circulant matrix, the eigenvectors are independent of the matrix entries
(they are always the DFT basis vectors). The eigenvalues depend on the
first row of the matrix.

Returns a vector of tuples [(λ₀, [v₀]), (λ₁, [v₁]), ..., (λₙ₋₁, [vₙ₋₁])].
Each eigenvector is wrapped in a single-element array to match the general
eigenpair format (which supports eigenspaces of dimension > 1).

If `normalized=true`, eigenvectors have unit norm.
"""
function _circulant_eigenpairs(mat; normalized::Bool=false)
    n = size(mat, 1)
    
    # Get eigenvalues
    eigenvalues = _circulant_eigenvalues(mat)
    
    # Get eigenvectors (same for all circulant matrices of this size)
    eigenvectors = _circulant_eigenvectors(n; normalized=normalized)
    
    # Package as eigenpairs
    pairs = Vector{Tuple{Any, Vector}}(undef, n)
    for k in 1:n
        pairs[k] = (eigenvalues[k], [eigenvectors[k]])
    end
    
    return pairs
end
