function _issymzero(x)
    try
        v = Base.iszero(x)
        if v isa Bool
            return v
        end
    catch
    end
    try
        sx = Symbolics.simplify(x)
        v = Symbolics.iszero(sx)
        return v === true
    catch
    end
    try
        v = Symbolics.iszero(x)
        return v === true
    catch
    end
    # Fall back to "not proven zero" to keep elimination moving.
    return false
end

function _rref(M)
    A = copy(Matrix(M))
    m, n = size(A)
    pivots = Int[]
    row = 1
    for col in 1:n
        pivot_row = _find_pivot(A, row, col)
        isnothing(pivot_row) && continue
        if pivot_row != row
            A[row, :], A[pivot_row, :] = A[pivot_row, :], A[row, :]
        end
        pivot = Symbolics.simplify(A[row, col])
        A[row, :] .= Symbolics.simplify.(A[row, :] ./ pivot)
        # Eliminate the current column in every other row.
        for r in 1:m
            r == row && continue
            factor = A[r, col]
            _issymzero(factor) && continue
            A[r, :] .= Symbolics.simplify.(A[r, :] .- factor .* A[row, :])
        end
        push!(pivots, col)
        row += 1
        row > m && break
    end
    return A, pivots
end

function _find_pivot(A, start_row, col)
    m = size(A, 1)
    for r in start_row:m
        !_issymzero(A[r, col]) && return r
    end
    return nothing
end

function _nullspace(M)
    simplified = Symbolics.simplify.(M)
    R, pivots = _rref(Symbolics.expand.(simplified))
    m, n = size(R)
    free = setdiff(1:n, pivots)
    vectors = Vector{Vector{eltype(R)}}()
    for f in free
        v = fill(zero(eltype(R)), n)
        v[f] = one(eltype(R))
        for (row, pivot_col) in enumerate(pivots)
            v[pivot_col] = Symbolics.simplify(-R[row, f])
        end
        push!(vectors, Symbolics.simplify.(v))
    end
    return vectors
end
