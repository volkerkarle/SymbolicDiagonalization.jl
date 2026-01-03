# Upstream Fix for Symbolics.jl: Complex{Num}^Real Power Bug

## Summary

`Complex{Num}^Real` operations fail with `UndefVarError: Pow not defined in Symbolics`.

```julia
using Symbolics
@variables x y
(x + y*im)^(1/3)  # ERROR: UndefVarError: `Pow` not defined in `Symbolics`
```

This prevents using Cardano's cubic formula (and similar) for symbolic matrices where complex intermediate values are required.

## Root Cause

In `src/solver/main.jl`, line 1:

```julia
Base.:^(a::Complex{<:Real}, b::Num) = Symbolics.Pow(a, b)
```

Two problems:

1. **`Symbolics.Pow` doesn't exist** - causes the `UndefVarError`
2. **Overly broad type signature** - `Complex{<:Real}` matches `Complex{Num}` since `Num <: Real`

### Dispatch Chain

For `Complex{Num}^Float64`:
1. `^(z::Complex{T}, p::S) where {T<:Real, S<:Real}` at `Base/complex.jl:885` matches
2. Calls `Complex{P}(z) ^ P(p)` with type promotion
3. Eventually hits the buggy `^(a::Complex{<:Real}, b::Num)` method
4. `Symbolics.Pow(...)` fails

## Proposed Fix (Option C)

Add an explicit `Complex{Num}^Real` method using polar form in `src/num.jl`:

```julia
# Add near line 196, alongside existing Complex{Num} operations
"""
    ^(z::Complex{Num}, p::Real)

Compute complex symbolic power using polar form:
z^p = |z|^p * exp(i*p*arg(z))

This is necessary because Base's complex power implementation uses
boolean conditionals that fail with symbolic values.
"""
function Base.:^(z::Complex{Num}, p::Real)
    a, b = reim(z)
    # r = |z| = sqrt(a² + b²)
    r = sqrt(a^2 + b^2)
    # θ = arg(z) = atan(b, a)
    θ = atan(b, a)
    # z^p = r^p * (cos(pθ) + i*sin(pθ))
    r_pow = r^p
    θ_pow = p * θ
    return Complex(r_pow * cos(θ_pow), r_pow * sin(θ_pow))
end

# Also add Rational power for completeness
Base.:^(z::Complex{Num}, p::Rational) = z^float(p)
```

### Why Polar Form?

- Avoids boolean conditionals (`isreal`, `iszero`, etc.) that fail on symbolic values
- Mathematically correct for all complex numbers
- Produces clean symbolic expressions with `atan`, `cos`, `sin`
- Already works because `atan(Num, Num)`, `cos(Num)`, `sin(Num)` are defined

## Test Cases

```julia
@testset "Complex{Num} power" begin
    @variables x y
    
    # Basic power operations
    z = Complex(x, y)
    @test z^2 == x^2 - y^2 + 2im*x*y  # Integer power (already works)
    
    # Non-integer powers should not error
    result = z^(1/3)
    @test result isa Complex{Num}
    
    # Numerical verification
    z_num = Complex(8.0, 0.0)
    r = real(z^(1/3))
    i = imag(z^(1/3))
    r_fn = Symbolics.build_function(r, [x, y]; expression=Val{false})
    i_fn = Symbolics.build_function(i, [x, y]; expression=Val{false})
    @test r_fn([8.0, 0.0]) ≈ 2.0
    @test i_fn([8.0, 0.0]) ≈ 0.0
    
    # Complex cube root
    @test r_fn([0.0, 8.0]) ≈ sqrt(3)  # real part of cbrt(8i)
    @test i_fn([0.0, 8.0]) ≈ 1.0      # imag part of cbrt(8i)
end
```

## Secondary Issue

The line in `solver/main.jl`:
```julia
Base.:^(a::Complex{<:Real}, b::Num) = Symbolics.Pow(a, b)
```

Should also be fixed. Assuming the intent was to create a symbolic power term:
```julia
Base.:^(a::Complex{T}, b::Num) where {T<:Union{AbstractFloat,Integer,Rational,AbstractIrrational}} = 
    Symbolics.wrap(SymbolicUtils.term(^, a, Symbolics.unwrap(b)))
```

This restricts the signature to concrete `Real` types, excluding `Num`.

## Related Issues

- #1016: `build_function` and `expand` don't work correctly on `Complex{Num}`
- #1487: `\` operator not working with symbolic complex variables
- #1661: Cannot convert BasicSymbolic{Number} to Complex{Num}

## Files to Modify

1. `src/num.jl` - Add `^(::Complex{Num}, ::Real)` method (primary fix)
2. `src/solver/main.jl` - Fix or restrict `^(::Complex{<:Real}, ::Num)` (secondary fix)
3. `test/complex.jl` or `test/num.jl` - Add test cases

## Verification

After the fix, this should work:

```julia
using Symbolics
@variables x y

# Direct power
z = Complex(x, y)
z^(1/3)  # Should return Complex{Num} with symbolic expressions

# Use case: Cardano's formula for cubic roots
# When discriminant Δ < 0 (casus irreducibilis), need Complex{Num}^(1/3)
Δ = x  # Could be negative symbolically
sqrtΔ = sqrt(Complex(Δ, zero(Δ)))  # Forces complex path
cbrt_term = (-y/2 + sqrtΔ)^(1//3)  # This should work
```
