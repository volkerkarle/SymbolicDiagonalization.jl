# ============================================================================
# Symbolic Simplification Utilities
# Aggressive rule-based rewriting for symbolic eigenvalue expressions
# ============================================================================

using SymbolicUtils
using SymbolicUtils.Rewriters

# ============================================================================
# Helper for Safe Rule Application
# ============================================================================

"""
    _safe_rule(rule)

Wrap a rewrite rule to catch BoundsError exceptions that occur when an @acrule
with N terms is applied to an expression with fewer than N terms.
Returns nothing instead of throwing an error.
"""
function _safe_rule(rule)
    function (x)
        try
            return rule(x)
        catch e
            if e isa BoundsError
                return nothing
            else
                rethrow(e)
            end
        end
    end
end

# ============================================================================
# Trigonometric Simplification Rules
# ============================================================================

"""
    TRIG_RULES

Core trigonometric identities for simplification.
"""
const TRIG_RULES = [
    # Pythagorean identity
    @acrule sin(~x)^2 + cos(~x)^2 => 1
    
    # Squared Pythagorean identity
    _safe_rule(@acrule sin(~x)^4 + 2*(sin(~x)^2)*(cos(~x)^2) + cos(~x)^4 => 1)
    _safe_rule(@acrule sin(~x)^4 + cos(~x)^4 + 2*(sin(~x)^2)*(cos(~x)^2) => 1)
    
    # Angle addition formulas
    @acrule cos(~a)*cos(~b) + (-1)*sin(~a)*sin(~b) => cos(~a + ~b)
    @acrule sin(~a)*cos(~b) + cos(~a)*sin(~b) => sin(~a + ~b)
    
    # Angle subtraction formulas  
    @acrule cos(~a)*cos(~b) + sin(~a)*sin(~b) => cos(~a - ~b)
    @acrule sin(~a)*cos(~b) + (-1)*cos(~a)*sin(~b) => sin(~a - ~b)
    
    # Negative sum
    @acrule (-1)*sin(~a)*cos(~b) + (-1)*cos(~a)*sin(~b) => -sin(~a + ~b)
    
    # Negative angle identities
    @rule sin(-1 * ~x) => -sin(~x)
    @rule cos(-1 * ~x) => cos(~x)
    @rule sin(-(~x)) => -sin(~x)
    @rule cos(-(~x)) => cos(~x)
    
    # Double angle formulas
    @rule 2*sin(~x)*cos(~x) => sin(2*(~x))
    @rule (-2)*sin(~x)*cos(~x) => -sin(2*(~x))
    @acrule cos(~x)^2 + (-1)*sin(~x)^2 => cos(2*(~x))
]

# ============================================================================
# Algebraic Simplification Rules
# ============================================================================

"""
    ALGEBRAIC_RULES

Algebraic identities useful for eigenvalue expressions.
"""
const ALGEBRAIC_RULES = [
    # sqrt(x^2) -> |x| (we assume positive for symbolic)
    @rule sqrt((~x)^2) => ~x
    
    # sqrt(a) * sqrt(b) -> sqrt(a*b)
    @rule sqrt(~a) * sqrt(~b) => sqrt((~a) * (~b))
    
    # (sqrt(x))^2 -> x
    @rule sqrt(~x)^2 => ~x
    
    # x / sqrt(x) -> sqrt(x) (when x > 0)
    @rule (~x) / sqrt(~x) => sqrt(~x)
    
    # sqrt(x) / x -> 1/sqrt(x)
    @rule sqrt(~x) / (~x) => 1 / sqrt(~x)
]

# ============================================================================
# Sqrt-Trig Simplification Rules (The Key Innovation)
# ============================================================================

"""
    SQRT_TRIG_RULES

Rules for simplifying sqrt expressions involving trig functions.
These are critical for getting clean eigenvalue expressions.

Key insight: sqrt(1 - cos^2(x)) = |sin(x)| = sin(x) for x in [0, pi]
"""
const SQRT_TRIG_RULES = [
    # sqrt(1 - cos^2(x)) -> sin(x) (assuming x in [0, pi])
    @rule sqrt(1 - cos(~x)^2) => sin(~x)
    @rule sqrt(1 + (-1)*cos(~x)^2) => sin(~x)
    @rule sqrt((-1)*cos(~x)^2 + 1) => sin(~x)
    
    # sqrt(1 - sin^2(x)) -> cos(x) (assuming x in [-pi/2, pi/2])
    @rule sqrt(1 - sin(~x)^2) => cos(~x)
    @rule sqrt(1 + (-1)*sin(~x)^2) => cos(~x)
    @rule sqrt((-1)*sin(~x)^2 + 1) => cos(~x)
    
    # More complex forms that appear in SO(4) eigenvalues
    # sqrt((1-cos(x))/2) -> |sin(x/2)| = sin(x/2) for x in [0, 2pi]
    @rule sqrt((1 - cos(~x)) / 2) => sin((~x) / 2)
    @rule sqrt((1 + cos(~x)) / 2) => cos((~x) / 2)
    
    # sqrt(2 - 2*cos(x)) -> 2*|sin(x/2)|
    @rule sqrt(2 - 2*cos(~x)) => 2*sin((~x) / 2)
    @rule sqrt(2 + 2*cos(~x)) => 2*cos((~x) / 2)
]

# ============================================================================
# Complex Number Simplification Rules
# ============================================================================

"""
    COMPLEX_RULES

Rules for simplifying complex expressions, especially eigenvalues.
"""
const COMPLEX_RULES = [
    # i^2 -> -1 (handled by Symbolics, but just in case)
    # |e^(ix)| -> 1
]

# ============================================================================
# Combined Rewriters
# ============================================================================

"""Internal rewriter for trigonometric simplification."""
const _trig_rewriter = Fixpoint(Prewalk(PassThrough(Chain(TRIG_RULES))))

"""Internal rewriter for algebraic simplification."""
const _algebraic_rewriter = Fixpoint(Prewalk(PassThrough(Chain(ALGEBRAIC_RULES))))

"""Internal rewriter for sqrt-trig simplification."""
const _sqrt_trig_rewriter = Fixpoint(Prewalk(PassThrough(Chain(SQRT_TRIG_RULES))))

"""Combined aggressive rewriter."""
const _aggressive_rewriter = Fixpoint(Prewalk(PassThrough(Chain([
    TRIG_RULES...,
    ALGEBRAIC_RULES...,
    SQRT_TRIG_RULES...
]))))

# ============================================================================
# Public Simplification Functions
# ============================================================================

"""
    trig_simplify(expr)

Apply trigonometric simplification rules to a symbolic expression.

# Examples
```julia
@variables θ φ
trig_simplify(sin(θ)^2 + cos(θ)^2)  # → 1
trig_simplify(cos(θ)*cos(φ) - sin(θ)*sin(φ))  # → cos(θ + φ)
```
"""
function trig_simplify(expr)
    if expr isa Num
        unwrapped = Symbolics.unwrap(expr)
        result = _trig_rewriter(unwrapped)
        return isnothing(result) ? expr : Symbolics.wrap(result)
    elseif expr isa Complex
        re_simp = trig_simplify(real(expr))
        im_simp = trig_simplify(imag(expr))
        return re_simp + im * im_simp
    else
        result = _trig_rewriter(expr)
        return isnothing(result) ? expr : result
    end
end

"""
    aggressive_simplify(expr)

Apply aggressive simplification combining:
- Trigonometric identities (sin²+cos²=1, angle addition, etc.)
- Algebraic simplifications (sqrt(x)^2 = x, etc.)
- Sqrt-trig rules (sqrt(1-cos²x) = sin(x))

This is designed for eigenvalue expressions where we want the cleanest
possible symbolic form.

# Examples
```julia
@variables θ
aggressive_simplify(sqrt(1 - cos(θ)^2))  # → sin(θ)
aggressive_simplify(cos(θ)^2 + sin(θ)^2)  # → 1
```
"""
function aggressive_simplify(expr)
    if expr isa Num
        unwrapped = Symbolics.unwrap(expr)
        # First expand, then apply rules, then simplify
        expanded = Symbolics.expand(Symbolics.wrap(unwrapped))
        unwrapped_exp = Symbolics.unwrap(expanded)
        result = _aggressive_rewriter(unwrapped_exp)
        if isnothing(result)
            return Symbolics.simplify(expanded)
        end
        simplified = Symbolics.wrap(result)
        return Symbolics.simplify(simplified)
    elseif expr isa Complex
        re_simp = aggressive_simplify(real(expr))
        im_simp = aggressive_simplify(imag(expr))
        return re_simp + im * im_simp
    else
        result = _aggressive_rewriter(expr)
        return isnothing(result) ? expr : result
    end
end

"""
    simplify_eigenvalue(expr)

Specialized simplification for eigenvalue expressions.

Applies multiple rounds of simplification to get the cleanest form:
1. Expand the expression
2. Apply trigonometric identities
3. Simplify sqrt expressions involving trig
4. Final cleanup with Symbolics.simplify

For rotation matrices, this should convert:
  sqrt((tr(A)-1)/2)^2 → (tr(A)-1)/2
  sqrt(1 - cos(θ)^2) → sin(θ)
  cos(θ) + i*sin(θ) → e^(iθ) (conceptually, keeps trig form)
"""
function simplify_eigenvalue(expr)
    # Handle complex eigenvalues
    if expr isa Complex
        re = simplify_eigenvalue(real(expr))
        im_part = simplify_eigenvalue(imag(expr))
        return re + im * im_part
    end
    
    if !(expr isa Num)
        return expr
    end
    
    # Multiple passes for thorough simplification
    result = expr
    for _ in 1:3
        # Expand first to expose structure
        result = Symbolics.expand(result)
        
        # Apply our custom rules
        unwrapped = Symbolics.unwrap(result)
        rewritten = _aggressive_rewriter(unwrapped)
        if !isnothing(rewritten)
            result = Symbolics.wrap(rewritten)
        end
        
        # Standard simplification
        result = Symbolics.simplify(result)
    end
    
    return result
end

"""
    simplify_eigenvalues(vals)

Apply simplification to a vector of eigenvalues.
"""
function simplify_eigenvalues(vals::AbstractVector)
    return [simplify_eigenvalue(v) for v in vals]
end

"""
    trig_simplify_matrix(mat)

Apply trigonometric simplification to all elements of a matrix.
"""
function trig_simplify_matrix(mat)
    return trig_simplify.(mat)
end

# ============================================================================
# Enhanced Zero Check with Trig Simplification
# ============================================================================

"""
    _issymzero_trig(expr)

Check if an expression is symbolically zero, using trigonometric simplification
to handle cases like sin²(θ) + cos²(θ) - 1.
"""
function _issymzero_trig(expr)
    if _issymzero(expr)
        return true
    end
    
    simplified = trig_simplify(expr)
    return _issymzero(simplified)
end

# ============================================================================
# Euler Form Conversion (Future Enhancement)
# ============================================================================

"""
    to_euler_form(expr)

Attempt to convert cos(θ) + i·sin(θ) to e^(iθ) form.

Note: Full Euler conversion requires pattern matching on the structure
of complex expressions. Currently returns simplified trig form.
"""
function to_euler_form(expr)
    return trig_simplify(expr)
end

# ============================================================================
# Utility: Detect if expression involves trig functions
# ============================================================================

"""
    _has_trig(expr)

Check if an expression contains trigonometric functions.
"""
function _has_trig(expr)
    if expr isa Num
        str = string(expr)
        return occursin("sin", str) || occursin("cos", str) || 
               occursin("tan", str) || occursin("cot", str)
    end
    return false
end
