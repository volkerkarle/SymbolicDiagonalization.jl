# ============================================================================
# Permutation Matrix Pattern Detection and Eigenvalue Computation
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
