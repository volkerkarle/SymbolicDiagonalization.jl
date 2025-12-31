# ============================================================================
# Interactive Pattern Exploration for SymbolicDiagonalization.jl
# ============================================================================
#
# This script demonstrates how to explore matrix patterns and their eigenvalues.
# Run interactively in the Julia REPL or as a script.
#
# Usage:
#   julia --project examples/explore_patterns.jl
#
# Or in the REPL:
#   include("examples/explore_patterns.jl")
# ============================================================================

using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

println("=" ^ 70)
println("SymbolicDiagonalization.jl - Pattern Exploration Examples")
println("=" ^ 70)

# ============================================================================
# Example 1: Simple 2×2 Symbolic Matrix
# ============================================================================

println("\n📊 Example 1: General 2×2 Symbolic Matrix")
println("-" ^ 50)

@variables a b c d
A = [a b; c d]
println("Matrix A = ")
display(A)

vals, poly, λ = symbolic_eigenvalues(A)
println("\nCharacteristic polynomial: ", poly)
println("Eigenvalues:")
for (i, v) in enumerate(vals)
    println("  λ$i = ", v)
end

# ============================================================================
# Example 2: Symmetric Matrix (Real Eigenvalues)
# ============================================================================

println("\n\n📊 Example 2: Symmetric 2×2 Matrix")
println("-" ^ 50)

@variables x y
S = [x y; y x]
println("Symmetric matrix S = ")
display(S)

vals_s, _, _ = symbolic_eigenvalues(S)
println("\nEigenvalues (always real for symmetric matrices):")
for (i, v) in enumerate(vals_s)
    println("  λ$i = ", simplify(v))
end

# ============================================================================
# Example 3: Circulant Matrix (Uses DFT-based formula)
# ============================================================================

println("\n\n📊 Example 3: Circulant Matrix")
println("-" ^ 50)

# Circulant matrix: each row is a cyclic shift of the first row
C = [1 2 3; 3 1 2; 2 3 1]
println("Circulant matrix C = ")
display(C)

vals_c, _, _ = symbolic_eigenvalues(C)
println("\nEigenvalues (computed via DFT formula):")
for (i, v) in enumerate(vals_c)
    println("  λ$i = ", v)
end

# ============================================================================
# Example 4: Permutation Matrix (Roots of Unity)
# ============================================================================

println("\n\n📊 Example 4: Permutation Matrix (3-cycle)")
println("-" ^ 50)

# Permutation matrix for cycle (1 2 3)
P = [0 1 0; 0 0 1; 1 0 0]
println("Permutation matrix P (cycle 1→2→3→1) = ")
display(P)

vals_p, _, _ = symbolic_eigenvalues(P)
println("\nEigenvalues (3rd roots of unity):")
for (i, v) in enumerate(vals_p)
    println("  λ$i = ", v)
end

# ============================================================================
# Example 5: Toeplitz Tridiagonal (Closed-form for any size)
# ============================================================================

println("\n\n📊 Example 5: Toeplitz Tridiagonal Matrix")
println("-" ^ 50)

# Toeplitz tridiagonal: constant diagonals
n = 5
T = diagm(-1 => fill(1.0, n-1), 0 => fill(2.0, n), 1 => fill(1.0, n-1))
println("Toeplitz tridiagonal T ($(n)×$(n)) = ")
display(T)

vals_t, _, _ = symbolic_eigenvalues(T)
println("\nEigenvalues (closed-form formula):")
for (i, v) in enumerate(vals_t)
    println("  λ$i = ", round(v, digits=6))
end

# ============================================================================
# Example 6: Block-Diagonal Matrix (Decomposed)
# ============================================================================

println("\n\n📊 Example 6: Block-Diagonal Matrix")
println("-" ^ 50)

@variables p q r s
B = [p 0 0 0;
     0 q 0 0;
     0 0 r s;
     0 0 s r]
println("Block-diagonal matrix B = ")
display(B)

vals_b, _, _ = symbolic_eigenvalues(B)
println("\nEigenvalues (each block solved independently):")
for (i, v) in enumerate(vals_b)
    println("  λ$i = ", simplify(v))
end

# ============================================================================
# Example 7: Hypercube Graph (Q₃ = 3-dimensional cube)
# ============================================================================

println("\n\n📊 Example 7: Hypercube Graph Q₃ (8 vertices)")
println("-" ^ 50)

# Adjacency matrix of 3-dimensional hypercube
Q3 = [0 1 1 0 1 0 0 0;
      1 0 0 1 0 1 0 0;
      1 0 0 1 0 0 1 0;
      0 1 1 0 0 0 0 1;
      1 0 0 0 0 1 1 0;
      0 1 0 0 1 0 0 1;
      0 0 1 0 1 0 0 1;
      0 0 0 1 0 1 1 0]
println("Hypercube Q₃ adjacency matrix (8×8):")
println("Eigenvalues are n-2k for k=0,1,...,n with multiplicity C(n,k)")

vals_q, _, _ = symbolic_eigenvalues(Q3)
println("\nEigenvalues:")
unique_vals = unique(vals_q)
for v in sort(unique_vals, rev=true)
    mult = count(x -> x == v, vals_q)
    println("  λ = $v  (multiplicity $mult)")
end

# ============================================================================
# Example 8: Full Diagonalization with Eigenvectors
# ============================================================================

println("\n\n📊 Example 8: Full Diagonalization (P, D, λ)")
println("-" ^ 50)

@variables α β
M = [α β; β α]
println("Matrix M = ")
display(M)

P_mat, D_mat, λ = symbolic_diagonalize(M)
println("\nDiagonal matrix D = ")
display(D_mat)
println("\nEigenvector matrix P = ")
display(P_mat)

# Verify: M * P = P * D
println("\nVerification: M*P - P*D should be zero:")
diff = simplify.(M * P_mat - P_mat * D_mat)
display(diff)

# ============================================================================
# Summary
# ============================================================================

println("\n\n" * "=" ^ 70)
println("Pattern Detection Summary")
println("=" ^ 70)
println("""
SymbolicDiagonalization.jl automatically detects these patterns:

📐 Structure-based:
   • Diagonal/Triangular - eigenvalues from diagonal entries
   • Block-diagonal - recursive decomposition
   • Persymmetric - half-size reduction for even dimensions

🔢 Special matrices:
   • Circulant - DFT-based closed form for any size
   • Block circulant - eigenvalues from block DFT
   • Toeplitz tridiagonal - trigonometric formula for any size
   • Permutation - roots of unity from cycle decomposition

📈 Graph patterns:
   • Hypercube Qₙ - eigenvalues n-2k with binomial multiplicities
   • Strongly regular graphs - three distinct eigenvalues from parameters

For matrices up to 4×4 without special structure, the package uses
Cardano's and Ferrari's formulas for exact closed-form solutions.
""")

println("Done! Explore more by modifying the examples above.")
