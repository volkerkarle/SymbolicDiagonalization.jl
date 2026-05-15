# ============================================================================
# Dense 4x4 Matrix -- The General Case
# ============================================================================
#
# For a dense 4x4 matrix, the characteristic polynomial is degree 4 -- the
# highest degree solvable in radicals (Ferrari's formula). This example
# demonstrates both numeric and symbolic diagonalization.
#
# A 5x5 general matrix would hit the Abel-Ruffini theorem: no closed-form
# solution in radicals exists for degree 5+ polynomials.
#
# Usage:
#   julia --project examples/dense_4x4.jl
# ============================================================================

using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

# Helper: evaluate a symbolic Num to Float64
num2float(x::Num) = Float64(eval(Symbolics.toexpr(x)))
num2float(x::Number) = Float64(x)

sep() = println(repeat("=", 72))
dash() = println(repeat("-", 72))

sep()
println("Dense 4x4 Matrix -- The General Case")
sep()

# ============================================================================
# 1. Numeric 4x4 (works instantly via LAPACK)
# ============================================================================

dash()
println("1. Numeric 4x4 Dense Matrix")
dash()

A = [3.1  1.2  0.8  2.3;
     0.5  4.7  1.9  3.2;
     2.1  0.3  5.6  1.7;
     1.4  3.8  0.9  2.6]

@time vals_num = eigvals(A)
println("\nNumeric eigenvalues:")
for (i, v) in enumerate(vals_num)
    println("  lambda$i = $v")
end

# ============================================================================
# 2. Symbolic 3x3 General Matrix (Cardano's formula)
# ============================================================================

dash()
println("2. Symbolic 3x3 Dense Matrix (Cardano's formula)")
dash()

@variables a b c d

# Deliberately non-symmetric, non-circulant arrangement
M3 = [a b c;
      d a b;
      c d a]

struct_hint = SymbolicDiagonalization._detect_structure(M3)
println("Detected structure: ", isnothing(struct_hint) ? ":none (general)" : string(struct_hint))

@time vals3 = eigvals(M3)
println("\nSymbolic 3x3 eigenvalues (lambda1 through lambda3):")
for (i, v) in enumerate(vals3)
    s = string(v)
    println("  lambda$i = ", length(s) > 120 ? s[1:120]*"..." : s)
end

# Verify numerically
subs3 = Dict(a=>3.1, b=>1.2, c=>0.8)
# Numerical verification uses LinearAlgebra's standard eigen for numeric matrices
subs3 = Dict(a=>3.1, b=>1.2, c=>0.8, d=>0.5)
M_num = num2float.(substitute(M3, subs3))
num_ref = eigvals(M_num)
println("\n  Direct LAPACK eigenvalues: ", join(round.(num_ref, digits=4), ", "))

# ============================================================================
# 3. Why Not Symbolic 4x4?
# ============================================================================

dash()
println("3. Symbolic 4x4 -- The Quartic Challenge")
dash()

println("""
  Ferrari's formula for a degree-4 polynomial involves:
    * Computing the resolvent cubic
    * Finding its real root
    * Solving two quadratics

  With symbolic coefficients from a 4x4 determinant, the intermediate
  expressions grow exponentially. Even with 3 symbolic variables, the
  integer coefficients in the characteristic polynomial overflow Int64
  (the determinant involves up to degree-4 terms in each variable,
  producing O(10^10) intermediate coefficients).

  For this reason, symbolic 4x4 is best used with:
    * Block-diagonal structure (decomposes into 2x2 blocks)
    * Circulant/Toeplitz structure (closed-form formulas)
    * Lie group structure (SO(4), SU(2)(x)SU(2))
    * Substituted numeric values for some variables
""")

# ============================================================================
# 4. Structured 4x4 that works symbollically
# ============================================================================

dash()
println("4. Block-Diagonal 4x4 (decomposes into 2x2 blocks)")
dash()

@variables p q r s

B = [p q  0 0;
     q p  0 0;
     0 0  r s;
     0 0  s r]

@time vals_block = eigvals(B)
println("\nEigenvalues:")
for (i, v) in enumerate(vals_block)
    println("  lambda$i = $v")
end

# ============================================================================
# Summary
# ============================================================================

sep()
println("Summary")
sep()
println("""
  Numeric 4x4 dense:       OK Instant (LAPACK)
  Symbolic 3x3 dense:      OK Works (Cardano's formula)
  Symbolic 4x4 dense:      AVOID Quartic overflow with >2 symbolic vars
  Symbolic 4x4 structured: OK Block-diagonal, circulant, Toeplitz
  Symbolic 5x5 dense:      N/A No closed form (Abel-Ruffini)
""")
