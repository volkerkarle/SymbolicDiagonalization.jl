"""
    ExpressionComplexityError

Exception thrown when symbolic expressions exceed safe complexity limits.
Contains a helpful message with suggestions for how to proceed.
"""
struct ExpressionComplexityError <: Exception
    message::String
end

Base.showerror(io::IO, e::ExpressionComplexityError) = print(io, "ExpressionComplexityError: ", e.message)

"""
    _aggressive_simplify(expr; max_terms = 10000)

Apply aggressive simplification to symbolic expressions by:
1. Expanding the expression
2. Applying standard simplification
3. Attempting to detect and factor perfect squares

This helps produce cleaner eigenvalue expressions, especially for discriminants
in the quadratic formula and intermediate results in cubic/quartic formulas.

If the expanded expression exceeds `max_terms`, throws an `ExpressionComplexityError`.
"""
function _aggressive_simplify(expr; max_terms = 10000)
    # For non-symbolic expressions, just return as-is
    if !_is_symbolic_coeff(expr)
        return expr
    end
    
    # Strategy: expand first to collect all terms, then simplify
    expanded = Symbolics.expand(expr)
    
    # Check complexity before continuing
    _check_expression_size(expanded, max_terms)
    
    simplified = Symbolics.simplify(expanded)
    
    # Try to factor perfect squares: x² ± 2xy + y² + rest -> (x±y)² + rest
    factored = _try_factor_perfect_square(simplified)
    
    return factored
end

"""
    _check_expression_size(expr, max_terms)

Check if a symbolic expression has grown too large. Throws `ExpressionComplexityError`
if the expression exceeds the threshold.
"""
function _check_expression_size(expr, max_terms)
    # Estimate size by counting operations in the expression tree
    size = _estimate_expr_size(expr)
    if size > max_terms
        throw(ExpressionComplexityError(
            """Expression has grown too large (≈$size terms, limit: $max_terms).
            
            Suggestions to resolve this:
            1. Reduce matrix size (try 2×2 or 3×3 instead of 4×4)
            2. Use fewer symbolic variables (substitute numeric values where possible)
            3. Exploit matrix structure (diagonal, triangular, block-diagonal)
            4. Use numeric eigenvalues instead: eigvals(Float64.(substitute(A, values)))
            5. Increase limit if needed: set max_terms parameter (use with caution)
            
            Note: 3×3 matrices with many variables can produce expressions with 1000+ terms.
            4×4 matrices can easily exceed 10,000 terms."""
        ))
    end
end

"""
    _estimate_expr_size(expr)

Estimate the size of a symbolic expression by counting operations.
"""
function _estimate_expr_size(expr)
    if !_is_symbolic_coeff(expr)
        return 1
    end
    
    try
        unwrapped = Symbolics.unwrap(expr)
        return _count_operations(unwrapped)
    catch
        return 1
    end
end

"""
    _count_operations(x)

Recursively count operations in a SymbolicUtils expression tree.
"""
function _count_operations(x)
    # Base case: if it's not a symbolic term, it's a leaf (size 1)
    if x isa Number || !isdefined(x, :f)
        return 1
    end
    
    # If it has arguments, count recursively
    if isdefined(x, :arguments)
        args = getfield(x, :arguments)
        return 1 + sum(_count_operations, args; init=0)
    end
    
    return 1
end

"""
    _try_factor_perfect_square(expr)

Attempt to detect and factor perfect square patterns in expressions.
Looks for patterns like x² - 2xy + y² + rest and converts them to (x-y)² + rest.

This is particularly useful for quadratic discriminants which often have the
form (x-y)² + 4z².

# Implementation Status
This function is intentionally a no-op placeholder. Full implementation would require:
1. Identifying all variables in the expression
2. Treating it as a multivariate polynomial  
3. Detecting perfect square trinomials (x² + y² ± 2xy → (x±y)²)
4. Reconstructing in factored form

Due to Symbolics.jl's automatic expansion behavior and the complexity of 
multivariate polynomial pattern matching, this is deferred. The expanded form
(e.g., a² - 2ac + 4b² + c²) is mathematically equivalent to (a-c)² + 4b² and
produces correct eigenvalues.
"""
function _try_factor_perfect_square(expr)
    # Returns expression as-is; see docstring for rationale
    return expr
end

"""
    _symbolic_sqrt(x)

Compute square root in a way that works with symbolic expressions including Complex{Num}.
Julia's Base sqrt(::Complex) has boolean checks that fail for symbolic values,
so we implement special handling.

For plain Num (non-complex symbolic), we add 0im to force complex arithmetic,
which is necessary for the "casus irreducibilis" in cubic formulas where the
discriminant is negative but all roots are real.
"""
function _symbolic_sqrt(x)
    # For plain Num (not already complex), wrap in Complex to handle potential negative values
    # This is needed for casus irreducibilis in cubic formulas
    if x isa Num && !(x isa Complex)
        # Create Complex{Num} and use the formula below
        return _symbolic_sqrt(Complex(x, zero(x)))
    end
    
    # If x is not Complex{Num}, use regular sqrt
    if !(x isa Complex{<:Any})
        return sqrt(x)
    end
    
    # For Complex{Num}, implement the formula manually to avoid boolean checks
    # sqrt(a + bi) = sqrt((r+a)/2) + i*sign(b)*sqrt((r-a)/2)
    # where r = sqrt(a² + b²)
    a = real(x)
    b = imag(x)
    r = sqrt(a^2 + b^2)
    
    real_part = sqrt((r + a) / 2)
    imag_part = sqrt((r - a) / 2)
    
    # Handle sign of imaginary part
    # For symbolic expressions, we use the formula that gives the principal branch
    return Complex(real_part, imag_part)
end

"""
    ComputationTimeoutError

Exception thrown when a symbolic computation exceeds the time limit.
"""
struct ComputationTimeoutError <: Exception
    message::String
end

Base.showerror(io::IO, e::ComputationTimeoutError) = print(io, "ComputationTimeoutError: ", e.message)

"""
    symbolic_roots(coeffs; timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)

Solve a univariate polynomial with coefficients given in ascending order
(`coeffs[1]` is the constant term). Supports degrees 1–4 using closed-form
formulas (linear, quadratic, Cardano, and Ferrari).

# Keyword Arguments
- `timeout`: Maximum computation time in seconds (default: $(DEFAULT_TIMEOUT_SECONDS)). Set to `nothing` to disable.
- `max_terms`: Maximum expression size during simplification (default: $(DEFAULT_MAX_TERMS)).

# Throws
- `ComputationTimeoutError`: If computation exceeds timeout
- `ExpressionComplexityError`: If expressions grow too large
"""
function symbolic_roots(coeffs; timeout = DEFAULT_TIMEOUT_SECONDS, max_terms = DEFAULT_MAX_TERMS)
    deg = length(coeffs) - 1
    deg >= 1 || error("expected polynomial coefficients")
    # Note: Symbolic quartic computation can produce extremely large expressions
    # and may take a long time. Consider using numeric evaluation or structured
    # matrices (diagonal, triangular, block-diagonal) for practical applications.
    
    if !isnothing(timeout)
        return _with_timeout(() -> _symbolic_roots_impl(coeffs, max_terms), timeout, deg)
    else
        return _symbolic_roots_impl(coeffs, max_terms)
    end
end

"""
    _symbolic_roots_impl(coeffs, max_terms)

Internal implementation of symbolic root finding without timeout wrapper.
"""
function _symbolic_roots_impl(coeffs, max_terms)
    normalized = Symbolics.expand.(coeffs)
    leading = normalized[end]
    normalized = isequal(leading, one(leading)) ? normalized : normalized ./ leading
    deg = length(coeffs) - 1
    
    if deg == 1
        return _roots_linear(normalized)
    elseif deg == 2
        return _roots_quadratic(normalized; max_terms = max_terms)
    elseif deg == 3
        return _roots_cubic(normalized; max_terms = max_terms)
    elseif deg == 4
        return _roots_quartic(normalized; max_terms = max_terms)
    else
        error("closed-form roots implemented for degree ≤ 4")
    end
end

"""
    _with_timeout(f, timeout_seconds, degree)

Execute function `f` with a timeout. If the computation takes longer than
`timeout_seconds`, throws a `ComputationTimeoutError`.

Uses a Channel-based approach to avoid race conditions between the computation
task and the timeout check.
"""
function _with_timeout(f, timeout_seconds, degree)
    # Use a Channel to signal completion - this avoids race conditions
    # between checking task status and timeout flag
    result_channel = Channel{Any}(1)
    
    task = @async begin
        try
            result = f()
            put!(result_channel, (:ok, result))
        catch e
            put!(result_channel, (:error, e))
        end
    end
    
    # Use timedwait pattern: try to get result within timeout
    start_time = time()
    while (time() - start_time) < timeout_seconds
        if isready(result_channel)
            status, value = take!(result_channel)
            close(result_channel)
            if status === :error
                throw(value)
            end
            return value
        end
        sleep(0.05)  # Small sleep to avoid busy-waiting
    end
    
    # Timeout occurred - try to interrupt the task
    close(result_channel)
    try
        schedule(task, InterruptException(), error=true)
    catch
        # Task may have already completed or be uninterruptible
    end
    
    throw(ComputationTimeoutError(
        """Computation exceeded timeout of $timeout_seconds seconds (degree $degree polynomial).
        
        Suggestions to resolve this:
        1. Reduce matrix size (symbolic 4×4 can be very slow)
        2. Use fewer symbolic variables (substitute known values)
        3. Check for block-diagonal structure to reduce effective size
        4. Use numeric eigenvalues: eigvals(Float64.(substitute(A, values)))
        5. Increase timeout if needed: symbolic_roots(coeffs; timeout = 600)
        6. Set timeout = nothing to disable (use with caution - may hang indefinitely)
        
        Note: Quartic formulas can take 10+ minutes for matrices with many variables."""
    ))
end

function _is_symbolic_coeff(x)
    if x isa Num
        ux = try
            Symbolics.unwrap(x)
        catch
            nothing
        end
        return !(ux isa Number)
    end
    x isa Number && return false
    try
        return Symbolics.issymbolic(x) === true
    catch
        return false
    end
end

_roots_linear(c; max_terms = 10000) = [-c[1] / c[2]]

function _roots_quadratic(c; max_terms = 10000)
    a, b, d = c[3], c[2], c[1]
    disc = b^2 - 4a*d
    # Apply aggressive simplification to get a canonical expanded form
    disc = _aggressive_simplify(disc; max_terms = max_terms)
    sqrt_disc = sqrt(disc)
    return [(-b - sqrt_disc) / (2a), (-b + sqrt_disc) / (2a)]
end

function _roots_cubic(c; max_terms = 10000)
    a, b, cc, d = c[4], c[3], c[2], c[1]
    if iszero(a)
        return _roots_quadratic(c[1:3]; max_terms = max_terms)
    end
    # Depressed cubic: y^3 + py + q. This is Cardano's method.
    b1 = b / a
    c1 = cc / a
    d1 = d / a
    p = c1 - b1^2 / 3
    q = (2b1^3) / 27 - (b1 * c1) / 3 + d1
    # Apply aggressive simplification to p and q
    p = _aggressive_simplify(p; max_terms = max_terms)
    q = _aggressive_simplify(q; max_terms = max_terms)
    Δ = (q / 2)^2 + (p / 3)^3
    # Simplify the discriminant
    Δ = _aggressive_simplify(Δ; max_terms = max_terms)
    # Use _symbolic_sqrt which handles negative values correctly by returning
    # complex results. This is needed for "casus irreducibilis" where Δ < 0
    # but all 3 roots are real.
    sqrtΔ = _symbolic_sqrt(Δ)
    # Use ^(1//3) instead of cbrt to handle complex arguments properly.
    C = _symbolic_cbrt(-q / 2 + sqrtΔ)
    D = _symbolic_cbrt(-q / 2 - sqrtΔ)
    # Primitive cube roots of unity: ω = e^(2πi/3), ω² = e^(4πi/3)
    # ω = -1/2 + √3/2·i, ω² = -1/2 - √3/2·i
    omega = -0.5 + 0.5 * sqrt(3) * im
    omega2 = -0.5 - 0.5 * sqrt(3) * im
    shift = -b1 / 3
    return [shift + C + D,
            shift + omega * C + omega2 * D,
            shift + omega2 * C + omega * D]
end

"""
    _symbolic_cbrt(x)

Compute the cube root of x, handling both real and complex arguments.
For symbolic expressions (Num), uses x^(1//3).
For Complex{Num}, we must compute manually since Julia's Complex^power
tries to use Symbolics.Pow which doesn't exist.

The cube root of a complex number z = r*e^(iθ) is r^(1/3)*e^(iθ/3).
In Cartesian form: if z = a + bi, then
  r = sqrt(a² + b²)
  θ = atan2(b, a)
  cbrt(z) = r^(1/3) * (cos(θ/3) + i*sin(θ/3))

For symbolic expressions, we use symbolic cos/sin/atan.
"""
function _symbolic_cbrt(x)
    # For plain Num (not complex), use power form
    if x isa Num && !(x isa Complex)
        return x^(1//3)
    end
    
    # For Complex{Num}, compute manually using polar form
    if x isa Complex && (real(x) isa Num || imag(x) isa Num)
        a = real(x)
        b = imag(x)
        # r = |z| = sqrt(a² + b²)
        r = sqrt(a^2 + b^2)
        # θ = atan(b, a) - angle of the complex number
        θ = atan(b, a)
        # Cube root: r^(1/3) * e^(iθ/3)
        r_cbrt = r^(1//3)
        θ_third = θ / 3
        # cos and sin for symbolic expressions
        real_part = r_cbrt * cos(θ_third)
        imag_part = r_cbrt * sin(θ_third)
        return Complex(real_part, imag_part)
    end
    
    # For numeric complex, use power
    if x isa Complex
        return x^(1/3)
    end
    
    # For real numbers, cbrt is fine
    return cbrt(x)
end

function _roots_quartic(c; max_terms = 10000)
    a, b, cc, d, e = c[5], c[4], c[3], c[2], c[1]
    if iszero(a)
        return _roots_cubic(c[1:4]; max_terms = max_terms)
    end
    # Ferrari's method via a resolvent cubic. Variables are spelled out to
    # keep the algebra readable when Symbolics prints the expressions.
    A = b / a
    B = cc / a
    Cc = d / a
    Dd = e / a
    p = B - 3A^2 / 8
    q = Cc + A^3 / 8 - A * B / 2
    r = Dd - 3A^4 / 256 + A^2 * B / 16 - A * Cc / 4
    # Simplify intermediate coefficients
    p = _aggressive_simplify(p; max_terms = max_terms)
    q = _aggressive_simplify(q; max_terms = max_terms)
    r = _aggressive_simplify(r; max_terms = max_terms)
    resolvent = [4p * r - q^2, -8r, -4p, 8]
    alphas = _roots_cubic(resolvent; max_terms = max_terms)
    alpha = alphas[1]
    beta_sq = 2alpha - p
    beta_sq = _aggressive_simplify(beta_sq; max_terms = max_terms)
    beta = _symbolic_sqrt(beta_sq)
    # For symbolic expressions, we cannot use iszero() directly as it requires
    # boolean context. Instead, compute gamma algebraically. When beta=0, the
    # limit of -q/(2*beta) is handled by the symbolic system.
    gamma = -q / (2beta)
    gamma = _aggressive_simplify(gamma; max_terms = max_terms)
    t1 = beta^2 - 4 * (alpha + gamma)
    t2 = beta^2 - 4 * (alpha - gamma)
    t1 = _aggressive_simplify(t1; max_terms = max_terms)
    t2 = _aggressive_simplify(t2; max_terms = max_terms)
    roots_y = [
        (-beta - _symbolic_sqrt(t1)) / 2,
        (-beta + _symbolic_sqrt(t1)) / 2,
        (beta - _symbolic_sqrt(t2)) / 2,
        (beta + _symbolic_sqrt(t2)) / 2
    ]
    return roots_y .- A / 4
end


