# ============================================================================
# Structure Detection and Matrix Property Analysis
# Functions to detect matrix properties and exploitable structure
# ============================================================================

"""
    _prepare_diagonal_shift(mat)

Prepare a diagonal shift transformation to reduce symbolic complexity.

For a matrix A with n symbolic variables, if A[n,n] is a simple symbolic variable `f`
that doesn't appear elsewhere, we can:
1. Shift: A' = A - f*I (makes A'[n,n] = 0)
2. Substitute: Replace compound diagonal terms (a-f) with fresh variables (ã)
3. Solve: Compute eigenvalues with fewer effective variables
4. Back-substitute: Replace fresh variables with original expressions

This reduces a 6-variable 3×3 problem to a 5-variable problem.

Returns `(shifted_mat, shift_val, back_substitution_dict)` if beneficial, otherwise `nothing`.
The back_substitution_dict maps fresh variables back to their original expressions.
"""
function _prepare_diagonal_shift(mat)
    n = size(mat, 1)
    n >= 2 || return nothing
    
    # Get the [n,n] element
    corner = mat[n, n]
    
    # Check if corner is a simple symbolic variable (not compound)
    if !(corner isa Num)
        return nothing
    end
    
    corner_vars = Symbolics.get_variables(corner)
    corner_vars_list = collect(corner_vars)
    if !(length(corner_vars_list) == 1 && isequal(corner, corner_vars_list[1]))
        return nothing  # Corner is not a simple variable
    end
    
    # Check if corner variable appears elsewhere in the matrix
    for i in 1:n, j in 1:n
        (i == n && j == n) && continue
        elem_vars = Symbolics.get_variables(mat[i, j])
        if corner in elem_vars
            return nothing  # Corner variable appears elsewhere, can't eliminate it
        end
    end
    
    # Shift the matrix: A' = A - corner*I
    shifted_mat = mat - corner * I
    
    # Now create fresh variables for diagonal elements that became compound (x - corner)
    # and build the back-substitution dictionary
    back_subs = Dict{Num, Any}()
    
    for i in 1:(n-1)
        orig_diag = mat[i, i]
        shifted_diag = shifted_mat[i, i]  # This is (orig_diag - corner)
        
        # Check if the original diagonal element has symbolic variables
        orig_vars = Symbolics.get_variables(orig_diag)
        if !isempty(orig_vars)
            # Create a fresh variable for this shifted diagonal entry
            # Use a name without underscore prefix to avoid Julia code generation issues
            fresh_var = Symbolics.variable(Symbol("ξ", i))  # Greek xi for shifted variables
            back_subs[fresh_var] = shifted_diag  # fresh_var → (orig - corner)
            
            # Replace in the shifted matrix
            shifted_mat[i, i] = fresh_var
        end
    end
    
    # If no substitutions were made (all diagonal elements were numeric), 
    # we still benefit from eliminating the corner variable
    return (shifted_mat, corner, back_subs)
end

"""
    _apply_back_substitution(exprs, back_subs)

Apply back-substitution to replace fresh variables with original expressions.
"""
function _apply_back_substitution(exprs, back_subs)
    if isempty(back_subs)
        return exprs
    end
    
    # Build substitution rules
    rules = [fresh => orig for (fresh, orig) in back_subs]
    
    return [Symbolics.substitute(expr, Dict(rules)) for expr in exprs]
end

"""
    _detect_structure(mat)

Detect the structure of a matrix for eigenvalue computation optimization.
Returns one of: :diagonal, :triangular, :hermitian, :symmetric, :unitary, :none
"""
function _detect_structure(mat)
    _is_diagonal(mat) && return :diagonal
    _is_triangular(mat) && return :triangular
    _is_hermitian(mat) && return :hermitian
    _is_symmetric(mat) && return :symmetric
    _is_unitary(mat) && return :unitary
    return :none
end

# ============================================================================
# Matrix Property Checks
# Boolean functions to check matrix properties
# ============================================================================

"""
    _is_diagonal(mat)

Check if a matrix is diagonal (all off-diagonal elements are zero).
"""
function _is_diagonal(mat)
    m, n = size(mat)
    for i in 1:m, j in 1:n
        i == j && continue
        !_issymzero(mat[i, j]) && return false
    end
    return true
end

"""
    _is_triangular(mat)

Check if a matrix is triangular (upper or lower).
"""
function _is_triangular(mat)
    m, n = size(mat)
    upper = true
    lower = true
    for i in 1:m, j in 1:n
        if i < j && !_issymzero(mat[i, j])
            upper = false
        elseif i > j && !_issymzero(mat[i, j])
            lower = false
        end
    end
    return upper || lower
end

"""
    _is_symmetric(mat)

Check if a matrix is symmetric (A = A^T).
"""
function _is_symmetric(mat)
    m, n = size(mat)
    m == n || return false
    for i in 1:m, j in i+1:n
        !_issymzero(mat[i, j] - mat[j, i]) && return false
    end
    return true
end

"""
    _is_hermitian(mat)

Check if a matrix is Hermitian (A = A^H).
"""
function _is_hermitian(mat)
    m, n = size(mat)
    m == n || return false
    for i in 1:m, j in i+1:n
        !_issymzero(mat[i, j] - conj(mat[j, i])) && return false
    end
    return true
end

"""
    _is_unitary(mat)

Check if a matrix is unitary (A * A^H = I).
"""
function _is_unitary(mat)
    m, n = size(mat)
    m == n || return false
    U = mat * adjoint(mat)
    for i in 1:m, j in 1:n
        target = i == j ? one(eltype(U)) : zero(eltype(U))
        !_issymzero(U[i, j] - target) && return false
    end
    return true
end

# ============================================================================
# Block Structure Detection
# ============================================================================

"""
    _block_split(mat)

Detect a block-diagonal split; returns the split index or nothing.
"""
function _block_split(mat)
    n = size(mat, 1)
    n <= 1 && return nothing
    for k in 1:n-1
        if all(_issymzero, mat[1:k, k+1:n]) && all(_issymzero, mat[k+1:n, 1:k])
            return k
        end
    end
    return nothing
end

"""
    _detect_multiple_blocks(mat)

Detect all block-diagonal structure in a matrix. Returns a vector of block ranges
[(start1, end1), (start2, end2), ...] or nothing if not block-diagonal.

This generalizes _block_split to handle multiple blocks. The algorithm finds maximal
contiguous blocks where each block has no interaction (zero off-diagonal blocks) with
other blocks.

Example:
```
[A 0 0]
[0 B 0]  →  [(1, 2), (3, 3), (4, 5)]
[0 0 C]
```
"""
function _detect_multiple_blocks(mat)
    n = size(mat, 1)
    n <= 1 && return nothing
    
    # Use a greedy algorithm: find the smallest independent block starting from position 1,
    # then recursively find blocks in the remaining matrix
    blocks = Tuple{Int,Int}[]
    current_start = 1
    
    while current_start <= n
        # Find the smallest k such that mat[current_start:k, *] and mat[*, current_start:k]
        # form an independent block (no interaction with indices k+1:n)
        block_end = n  # Default: extend to the end if no boundary found
        
        for k in current_start:(n-1)
            # Check if there's a block boundary after position k
            # i.e., mat[current_start:k, k+1:n] and mat[k+1:n, current_start:k] are zero
            upper_right_zero = all(_issymzero, mat[i, j] for i in current_start:k, j in (k+1):n)
            lower_left_zero = all(_issymzero, mat[i, j] for i in (k+1):n, j in current_start:k)
            
            if upper_right_zero && lower_left_zero
                # Found a boundary! This is the end of the current block
                block_end = k
                break
            end
        end
        
        push!(blocks, (current_start, block_end))
        current_start = block_end + 1
    end
    
    # Only return if we found actual blocks (more than just the whole matrix)
    return length(blocks) > 1 ? blocks : nothing
end

# ============================================================================
# Persymmetric Structure Detection
# ============================================================================

"""
    _is_persymmetric(mat)

Check if a matrix is persymmetric: Q[i,j] == Q[n+1-j, n+1-i]
"""
function _is_persymmetric(mat)
    n = size(mat, 1)
    for i in 1:n, j in 1:n
        if !isequal(Symbolics.simplify(mat[i, j] - mat[n+1-j, n+1-i]), 0)
            return false
        end
    end
    return true
end

"""
    _persymmetric_split(mat)

For symmetric persymmetric matrices, compute a transformation that reduces
the eigenvalue problem. Returns `(P, S, D)` where:
- `P` is the transformation matrix
- `S` and `D` are the symmetric and skew-symmetric parts
- Eigenvalues of `mat` = eigenvalues of `S` ∪ eigenvalues of `D`

If the matrix is not symmetric persymmetric, returns `nothing`.
"""
function _persymmetric_split(mat)
    n = size(mat, 1)
    n % 2 != 0 && return nothing  # Only works for even dimensions
    # Check if symmetric or Hermitian (both work with this decomposition)
    if !(_is_symmetric(mat) || _is_hermitian(mat))
        return nothing
    end
    !_is_persymmetric(mat) && return nothing
    
    # For symmetric persymmetric matrices, we can use the fact that they
    # commute with the exchange matrix J. The eigenspaces of J give us
    # a block decomposition.
    
    # Exchange matrix J (anti-identity)
    J = zeros(eltype(mat), n, n)
    for i in 1:n
        J[i, n+1-i] = one(eltype(mat))
    end
    
    # Build transformation matrix P = [(I+J)/√2, (I-J)/√2]
    # But for symbolic computation, avoid √2 and just use [I+J, I-J]
    # The eigenvalues are the same, eigenvectors are just scaled
    I_n = Matrix{eltype(mat)}(I, n, n)
    
    # For a symmetric persymmetric matrix Q:
    # Q*J = J*Q (they commute)
    # So Q can be decomposed using eigenvectors of J
    
    # J has eigenvalues +1 and -1 (for even n)
    # Eigenvectors of eigenvalue +1: (e_i + e_{n+1-i})/√2
    # Eigenvectors of eigenvalue -1: (e_i - e_{n+1-i})/√2
    
    # Compute Q in the transformed basis
    # This gives us two smaller blocks
    
    # For n=4: P = [e1+e4, e2+e3, e1-e4, e2-e3]
    # In this basis, Q becomes block diagonal
    
    half = div(n, 2)
    P = zeros(eltype(mat), n, n)
    
    # First half of columns: e_i + e_{n+1-i} for i=1:half
    for i in 1:half
        P[i, i] = one(eltype(mat))
        P[n+1-i, i] = one(eltype(mat))
    end
    
    # Second half of columns: e_i - e_{n+1-i} for i=1:half  
    for i in 1:half
        P[i, half+i] = one(eltype(mat))
        P[n+1-i, half+i] = -one(eltype(mat))
    end
    
    # Transform: Q_new = P^T * Q * P
    # Note: P^T * P = 2I (since we didn't normalize by 1/√2), so the eigenvalues
    # of Q_transformed are 2× the eigenvalues of Q. We need to divide by 2.
    Q_transformed = transpose(P) * mat * P
    
    # Extract the two blocks
    block1 = Q_transformed[1:half, 1:half]
    block2 = Q_transformed[half+1:end, half+1:end]
    
    # Divide by 2 to account for P not being normalized
    # (P^T * P = 2I instead of I, so eigenvalues are scaled by 2)
    block1 = Symbolics.simplify.(block1 ./ 2)
    block2 = Symbolics.simplify.(block2 ./ 2)
    
    # For Hermitian matrices, the blocks are real but may have Complex{Num} type
    # Convert to real Num if all imaginary parts are zero
    if eltype(mat) <: Complex
        # Check if blocks are actually real and convert if so
        block1_real = all(x -> isequal(Symbolics.simplify(imag(x)), 0), block1)
        block2_real = all(x -> isequal(Symbolics.simplify(imag(x)), 0), block2)
        
        if block1_real
            block1 = real.(block1)
        end
        if block2_real
            block2 = real.(block2)
        end
    end
    
    return (block1, block2, P)
end

# ============================================================================
# Complexity Analysis
# ============================================================================

"""
    _count_symbolic_vars(A)

Count distinct symbolic variables appearing in matrix entries.
"""
function _count_symbolic_vars(A)
    vars = Set{Any}()
    for elem in A
        if elem isa Num
            try
                syms = Symbolics.get_variables(elem)
                for s in syms
                    push!(vars, s)
                end
            catch
            end
        end
    end
    return length(vars)
end

"""
    _check_complexity(A; threshold = 5, quiet = false)

Warn about complexity based on number of distinct symbolic variables.
"""
function _check_complexity(A; threshold = 5, quiet = false)
    n_vars = _count_symbolic_vars(A)
    if !quiet && n_vars > threshold
        @warn "Matrix contains $n_vars symbolic variables (threshold: $threshold). Operations may be slow or timeout with many parameters. Consider structured matrices, partial substitution, or expand=false."
    end
    return n_vars
end

# ============================================================================
# Adjugate Matrix Computation
# Helper functions for computing adjugate matrices (used in eigenvector computation)
# ============================================================================

"""
    _adjugate_vectors(M)

Compute eigenvectors using the adjugate matrix method.
Only used for small matrices (n ≤ MAX_ADJUGATE_SIZE) due to O(n!) complexity.
"""
function _adjugate_vectors(M)
    n = size(M, 1)
    # Adjugate computation has O(n!) complexity, only use for small matrices
    n <= MAX_ADJUGATE_SIZE || return Vector{Vector{eltype(M)}}()
    adj = _adjugate(M)
    vecs = Vector{Vector{eltype(M)}}()
    for j in 1:n
        col = Symbolics.simplify.(adj[:, j])
        all(_issymzero, col) && continue
        push!(vecs, col)
        # Return immediately after finding the first non-zero column to avoid
        # collecting multiple linearly dependent eigenvectors from the adjugate.
        return vecs
    end
    return vecs
end

"""
    _adjugate(M)

Compute the adjugate (classical adjoint) of matrix M.
"""
function _adjugate(M)
    n = size(M, 1)
    adj = Matrix{eltype(M)}(undef, n, n)
    if n == 1
        adj[1, 1] = one(eltype(M))
        return adj
    end
    for i in 1:n, j in 1:n
        minor_det = _minor_det(M, i, j)
        adj[j, i] = (-one(eltype(M)))^(i + j) * minor_det
    end
    return Symbolics.simplify.(adj)
end

"""
    _minor_det(M, row, col)

Compute the (row, col) minor determinant of matrix M.
"""
function _minor_det(M, row, col)
    n = size(M, 1)
    if n == 1
        return one(eltype(M))
    elseif n == 2
        r = setdiff(1:2, row)
        c = setdiff(1:2, col)
        return M[r[1], c[1]]
    elseif n == 3
        rows = setdiff(1:3, row)
        cols = setdiff(1:3, col)
        a, b = rows
        c, d = cols
        return Symbolics.simplify(M[a, c] * M[b, d] - M[a, d] * M[b, c])
    else
        error("minor_det only implemented for n ≤ 3")
    end
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    _diag_error_msg(struct_hint, found, n)

Generate an error message for non-diagonalizable matrices.
"""
function _diag_error_msg(struct_hint, found, n)
    base = "matrix is not diagonalizable: found $found independent eigenvectors for size $n"
    if struct_hint in (:hermitian, :symmetric, :unitary)
        return "$base (note: $struct_hint matrices should be diagonalizable; check symbolic simplification or assumptions)"
    end
    return base
end

"""
    _build_eigenvalue_result(vals, λ, expand)

Build the standard return tuple for eigenvalue computations.
Returns `(vals, poly, λ)` where `poly` is the characteristic polynomial.
"""
function _build_eigenvalue_result(vals, λ, expand)
    poly_raw = prod(λ .- vals)
    poly = expand ? Symbolics.expand(poly_raw) : poly_raw
    return vals, poly, λ
end

# ============================================================================
# Numerical Evaluation Utilities
# ============================================================================

"""
    _safe_complex_sqrt(x)

Compute sqrt that handles tiny negative numbers due to floating-point precision.
Returns a complex result when the argument is negative.
"""
function _safe_complex_sqrt(x::Real)
    if x >= 0
        return Complex(sqrt(x), 0.0)
    else
        # Negative argument - return purely imaginary result
        return Complex(0.0, sqrt(-x))
    end
end
_safe_complex_sqrt(x::Complex) = sqrt(x)

"""
    _evaluate_symbolic_expr(expr, var_values::Dict)

Evaluate a symbolic expression with numerical values, handling complex arithmetic
robustly to avoid domain errors from sqrt of tiny negative numbers.

# Arguments
- `expr`: A symbolic expression (Num or Complex{Num})
- `var_values`: Dict mapping symbolic variables to their numerical values

# Returns
A Complex{Float64} value

# Example
```julia
@variables a b c
expr = sqrt(a^2 + b^2) + c
result = _evaluate_symbolic_expr(expr, Dict(a => 3.0, b => 4.0, c => 1.0))
# result == 6.0 + 0.0im
```

# Notes
This function is designed to handle expressions from cubic/quartic formulas where
floating-point precision can cause discriminants that should be zero to become
tiny negative numbers. It uses complex arithmetic throughout to avoid NaN results.
"""
function _evaluate_symbolic_expr(expr, var_values::Dict)
    # Handle Complex{Num}
    re_expr = real(expr)
    im_expr = imag(expr)
    
    # Substitute values
    re_subbed = Symbolics.substitute(re_expr, var_values)
    im_subbed = Symbolics.substitute(im_expr, var_values)
    
    # Convert to Julia expression and evaluate with complex-safe functions
    re_val = _eval_with_complex_sqrt(re_subbed)
    im_val = _eval_with_complex_sqrt(im_subbed)
    
    return Complex(re_val, im_val)
end

"""
    _eval_with_complex_sqrt(expr)

Evaluate a substituted symbolic expression using complex-safe sqrt.
"""
function _eval_with_complex_sqrt(expr)
    # Handle Num type first (since Num <: Number in Symbolics)
    if expr isa Num
        # Try to extract the underlying value
        val = Symbolics.value(expr)
        # Check if it's a plain number (not symbolic)
        if val isa AbstractFloat || val isa Integer || val isa Rational
            return Float64(val)
        end
        # Otherwise, convert to Julia expression and evaluate
        julia_expr = Symbolics.toexpr(expr)
        result = _eval_expr_safe(julia_expr)
        return real(result)
    end
    
    # If it's already a plain Julia number, return it
    if expr isa Real
        return Float64(expr)
    end
    
    # Shouldn't reach here, but handle gracefully
    julia_expr = Symbolics.toexpr(expr)
    result = _eval_expr_safe(julia_expr)
    return real(result)
end

"""
    _eval_expr_safe(expr)

Recursively evaluate an expression tree, using complex-safe sqrt.
"""
function _eval_expr_safe(expr)
    if expr isa Number
        return Complex{Float64}(expr)
    elseif expr isa Symbol
        # Should not happen after substitution, but handle gracefully
        error("Unsubstituted variable: $expr")
    elseif expr isa Expr
        if expr.head == :call
            fn = expr.args[1]
            args = [_eval_expr_safe(a) for a in expr.args[2:end]]
            
            # Handle both Symbol (:sqrt) and function object (sqrt) cases
            fn_name = fn isa Symbol ? fn : nameof(fn)
            
            if fn_name == :sqrt || fn_name == :√
                return sqrt(Complex{Float64}(args[1]))
            elseif fn_name == :cbrt || fn_name == :∛
                # cbrt of negative real is real, but we use complex for consistency
                return cbrt(Complex{Float64}(real(args[1])))
            elseif fn_name == :^
                base, exp = args
                if exp isa Rational && exp.den != 1
                    # Fractional power - use complex
                    return Complex{Float64}(base)^exp
                else
                    return base^exp
                end
            elseif fn_name == :+
                return sum(args)
            elseif fn_name == :-
                if length(args) == 1
                    return -args[1]
                else
                    return args[1] - args[2]
                end
            elseif fn_name == :*
                return prod(args)
            elseif fn_name == :/
                return args[1] / args[2]
            elseif fn_name == :sin
                return sin(args[1])
            elseif fn_name == :cos
                return cos(args[1])
            elseif fn_name == :tan
                return tan(args[1])
            elseif fn_name == :atan
                if length(args) == 1
                    return atan(args[1])
                else
                    # Two-argument atan(y, x) - use real parts since this is for angle computation
                    return Complex{Float64}(atan(real(args[1]), real(args[2])))
                end
            elseif fn_name == :acos
                return acos(args[1])
            elseif fn_name == :asin
                return asin(args[1])
            elseif fn_name == :exp
                return exp(args[1])
            elseif fn_name == :log
                return log(Complex{Float64}(args[1]))
            elseif fn_name == :abs
                return abs(args[1])
            elseif fn_name == :real
                return real(args[1])
            elseif fn_name == :imag
                return imag(args[1])
            elseif fn_name == :conj
                return conj(args[1])
            else
                # Try to call the function directly
                if fn isa Function
                    return fn(args...)
                else
                    f = getfield(Base, fn_name)
                    return f(args...)
                end
            end
        elseif expr.head == :if || expr.head == :elseif
            # Handle conditional expressions
            cond = _eval_expr_safe(expr.args[1])
            if real(cond) != 0
                return _eval_expr_safe(expr.args[2])
            elseif length(expr.args) >= 3
                return _eval_expr_safe(expr.args[3])
            else
                return Complex{Float64}(0)
            end
        else
            error("Unsupported expression head: $(expr.head)")
        end
    elseif expr isa Function
        # Function objects can appear in expressions, just return them for later use
        # This shouldn't happen in numerical evaluation context
        error("Unexpected function object in expression: $expr")
    else
        return Complex{Float64}(expr)
    end
end
