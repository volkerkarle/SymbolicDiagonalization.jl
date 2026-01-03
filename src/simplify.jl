# ============================================================================
# Symbolic Simplification Utilities
# Rule-based rewriting for trigonometric and other expressions
# ============================================================================

using SymbolicUtils
using SymbolicUtils.Rewriters

# ============================================================================
# Trigonometric Simplification Rules
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

"""
    TRIG_RULES

A collection of rewrite rules for trigonometric simplification:
- Pythagorean identity: sin²(x) + cos²(x) → 1
- Angle addition: cos(a)cos(b) - sin(a)sin(b) → cos(a+b)
- Angle subtraction: cos(a)cos(b) + sin(a)sin(b) → cos(a-b)
- Double angle formulas
- Negative angle handling
"""
const TRIG_RULES = [
    # Pythagorean identity
    @acrule sin(~x)^2 + cos(~x)^2 => 1
    
    # Squared Pythagorean identity: (sin² + cos²)² = sin⁴ + 2sin²cos² + cos⁴ = 1
    # These rules match 3-term sums, so they need the safe wrapper to avoid
    # BoundsError when applied to 2-term expressions
    _safe_rule(@acrule sin(~x)^4 + 2*(sin(~x)^2)*(cos(~x)^2) + cos(~x)^4 => 1)
    _safe_rule(@acrule sin(~x)^4 + cos(~x)^4 + 2*(sin(~x)^2)*(cos(~x)^2) => 1)
    
    # Angle addition formulas
    # cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
    @acrule cos(~a)*cos(~b) + (-1)*sin(~a)*sin(~b) => cos(~a + ~b)
    # sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
    @acrule sin(~a)*cos(~b) + cos(~a)*sin(~b) => sin(~a + ~b)
    
    # Angle subtraction formulas  
    # cos(a-b) = cos(a)cos(b) + sin(a)sin(b)
    @acrule cos(~a)*cos(~b) + sin(~a)*sin(~b) => cos(~a - ~b)
    # sin(a-b) = sin(a)cos(b) - cos(a)sin(b)
    @acrule sin(~a)*cos(~b) + (-1)*cos(~a)*sin(~b) => sin(~a - ~b)
    
    # Negative sum: -sin(a)cos(b) - cos(a)sin(b) = -sin(a+b)
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

"""
    _trig_rewriter

Internal rewriter that applies trigonometric simplification rules.
Uses Fixpoint to apply rules until no more changes occur.
"""
const _trig_rewriter = Fixpoint(Prewalk(PassThrough(Chain(TRIG_RULES))))

"""
    trig_simplify(expr)

Apply trigonometric simplification rules to a symbolic expression.

Simplifies expressions using identities like:
- sin²(x) + cos²(x) = 1
- cos(a)cos(b) - sin(a)sin(b) = cos(a+b)
- sin(a)cos(b) + cos(a)sin(b) = sin(a+b)

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
        # Try to apply directly
        result = _trig_rewriter(expr)
        return isnothing(result) ? expr : result
    end
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
    # First try standard simplification
    if _issymzero(expr)
        return true
    end
    
    # Try with trig simplification
    simplified = trig_simplify(expr)
    return _issymzero(simplified)
end

# ============================================================================
# Complex/Euler Form Simplification
# ============================================================================

"""
    EULER_RULES

Rules for converting between trigonometric and exponential forms.
"""
const EULER_RULES = [
    # e^(ix) = cos(x) + i*sin(x) - can't easily represent this
    # But we can simplify products of complex exponentials
]

"""
    to_euler_form(expr)

Attempt to convert cos(θ) + i·sin(θ) to e^(iθ) form.
Currently just returns the simplified trig form.
"""
function to_euler_form(expr)
    # For now, just apply trig simplification
    # Full Euler conversion would require detecting cos + i*sin patterns
    return trig_simplify(expr)
end
