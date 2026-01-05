# ============================================================================
# Finite Group: Graph Patterns with Algebraic Eigenvalues
# ============================================================================
#
# Hypercube Graphs Qₙ:
#   The n-dimensional hypercube is the Cayley graph of (Z₂)ⁿ. Since (Z₂)ⁿ
#   is abelian, the Walsh-Hadamard matrix diagonalizes Qₙ. Eigenvalues are
#   λₖ = n - 2k with multiplicity C(n,k), giving closed form for any dimension.
#
# Strongly Regular Graphs srg(n, k, λ, μ):
#   The regularity constraints force exactly 3 distinct eigenvalues, which
#   can be computed from the parameters via a quadratic formula. This
#   algebraic structure avoids solving the characteristic polynomial.
#
# Both patterns exploit group-theoretic or algebraic constraints to bypass
# the Abel-Ruffini limitation on polynomial root-finding.
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
    eigenvalues = [Float64(k)]  # Multiplicity 1
    for _ in 1:f
        push!(eigenvalues, r)
    end
    for _ in 1:g
        push!(eigenvalues, s)
    end
    
    return eigenvalues
end
