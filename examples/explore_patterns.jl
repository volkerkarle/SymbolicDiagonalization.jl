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

println("\n\n" * "=" ^ 70)
println("BONUS: Matrices That Are 'Impossible' to Diagonalize Symbolically")
println("=" ^ 70)
println("""
The Abel-Ruffini theorem proves that polynomials of degree 5+ have no 
general closed-form solution. Yet these matrices work because of structure!
""")

# ============================================================================
# Example 9: 10×10 Circulant Matrix
# ============================================================================

println("\n\n📊 Example 9: 10×10 Circulant Matrix (Degree-10 Polynomial)")
println("-" ^ 50)

# Build a 10×10 circulant matrix
n_circ = 10
first_row = collect(1:n_circ)
C10 = [first_row[mod(j - i, n_circ) + 1] for i in 1:n_circ, j in 1:n_circ]
println("Circulant matrix C (10×10) with first row [1, 2, 3, ..., 10]")
println("This requires solving a degree-10 polynomial - no closed form in general!")

vals_c10, _, _ = symbolic_eigenvalues(C10)
println("\nEigenvalues computed via DFT formula (instant, closed-form):")
println("  λ₁ = $(round(real(vals_c10[1]), digits=4))")
println("  λ₅ = $(round(real(vals_c10[5]), digits=4) + round(imag(vals_c10[5]), digits=4))im")
println("  Total eigenvalues: $(length(vals_c10))")
println("\nNote: Works for ANY size - 100×100, 1000×1000, etc!")

# ============================================================================
# Example 10: 16×16 Hypercube Graph Q₄
# ============================================================================

println("\n\n📊 Example 10: 16×16 Hypercube Graph Q₄")
println("-" ^ 50)

# Build Q₄ adjacency matrix (4-dimensional hypercube = 16 vertices)
function build_hypercube(n)
    N = 2^n
    Q = zeros(Int, N, N)
    for i in 0:N-1
        for j in 0:N-1
            # Adjacent if differ by exactly one bit
            if count_ones(xor(i, j)) == 1
                Q[i+1, j+1] = 1
            end
        end
    end
    return Q
end

Q4 = build_hypercube(4)
println("Hypercube Q₄: 16 vertices, each connected to 4 neighbors")
println("General 16×16 matrix requires degree-16 polynomial - impossible!")

vals_q4, _, _ = symbolic_eigenvalues(Q4)
unique_vals_q4 = unique(vals_q4)
println("\nEigenvalues (closed-form via Walsh-Hadamard theory):")
for v in sort(unique_vals_q4, rev=true)
    mult = count(x -> x == v, vals_q4)
    println("  λ = $v  with multiplicity $mult = C(4, $(div(4-v,2)))")
end
println("\nNote: Works for ANY dimension - Q₆ (64×64), Q₁₀ (1024×1024), etc!")

# ============================================================================
# Example 11: 20×20 Path Graph Laplacian
# ============================================================================

println("\n\n📊 Example 11: 20×20 Path Graph Laplacian")
println("-" ^ 50)

P20 = path_laplacian(20)
println("Path graph P₂₀ Laplacian (tridiagonal, degree-20 characteristic polynomial)")
println("Structure: [1,-1,0,...], [-1,2,-1,0,...], ..., [...,0,-1,1]")

vals_p20, _, _ = symbolic_eigenvalues(P20)
println("\nEigenvalues (closed-form: λₖ = 2 - 2cos(πk/20)):")
println("  λ₁ = $(round(vals_p20[1], digits=6))  (algebraic connectivity)")
println("  λ₁₀ = $(round(vals_p20[10], digits=6))")
println("  λ₂₀ = $(round(vals_p20[20], digits=6))  (spectral radius)")
println("\nNote: Works for ANY size - 50×50, 100×100, 1000×1000!")

# ============================================================================
# Example 12: 16×16 Hadamard Matrix
# ============================================================================

println("\n\n📊 Example 12: 16×16 Hadamard Matrix")
println("-" ^ 50)

H16 = hadamard_matrix(4)  # 2^4 = 16
println("Sylvester-Hadamard H₁₆: 16×16 matrix of ±1 entries")
println("Degree-16 characteristic polynomial, but only 2 distinct eigenvalues!")

vals_h16, _, _ = symbolic_eigenvalues(H16)
unique_h16 = unique(vals_h16)
println("\nEigenvalues (closed-form: λ = ±√16 = ±4):")
for v in sort(unique_h16, rev=true)
    mult = count(x -> x == v, vals_h16)
    println("  λ = $(round(v, digits=4)) with multiplicity $mult")
end
println("\nNote: Works for ANY 2ⁿ - H₃₂, H₆₄, H₁₀₂₄, etc!")

# ============================================================================
# Example 13: 8×8 DFT (Fourier) Matrix  
# ============================================================================

println("\n\n📊 Example 13: 8×8 DFT Matrix")
println("-" ^ 50)

F8 = dft_matrix(8)
println("Discrete Fourier Transform matrix F₈ (8×8 complex)")
println("Unitary matrix with entries Fⱼₖ = ω^(jk) where ω = e^(2πi/8)")

vals_f8, _, _ = symbolic_eigenvalues(F8)
println("\nEigenvalues (exactly 4 distinct: ±√8, ±i√8):")
# Group by value
sqrt8 = sqrt(8)
counts = Dict{String, Int}()
for v in vals_f8
    if abs(real(v) - sqrt8) < 0.01 && abs(imag(v)) < 0.01
        key = "+√8"
    elseif abs(real(v) + sqrt8) < 0.01 && abs(imag(v)) < 0.01
        key = "-√8"
    elseif abs(imag(v) - sqrt8) < 0.01 && abs(real(v)) < 0.01
        key = "+i√8"
    elseif abs(imag(v) + sqrt8) < 0.01 && abs(real(v)) < 0.01
        key = "-i√8"
    else
        key = "other"
    end
    counts[key] = get(counts, key, 0) + 1
end
for (k, v) in sort(collect(counts), by=x->x[1])
    println("  λ = $k with multiplicity $v")
end

# ============================================================================
# Example 14: Cartan Matrix Type A₇ (7×7)
# ============================================================================

println("\n\n📊 Example 14: 7×7 Cartan Matrix (Type A₇)")
println("-" ^ 50)

A7 = cartan_matrix_A(7)
println("Cartan matrix of type A₇ (root system of SU(8))")
println("Tridiagonal with 2 on diagonal, -1 on off-diagonals")
println("display(A7):")
display(A7)

vals_a7, _, _ = symbolic_eigenvalues(A7)
println("\nEigenvalues (closed-form: λₖ = 4sin²(πk/16)):")
for (k, v) in enumerate(vals_a7)
    println("  λ$k = $(round(v, digits=6))")
end
println("\nNote: Works for ANY rank - A₁₀, A₅₀, A₁₀₀!")

# ============================================================================
# Example 15: 64×64 Nested Kronecker (6-fold)
# ============================================================================

println("\n\n📊 Example 15: 64×64 Nested Kronecker Product (6-fold)")
println("-" ^ 50)

# Create 6 diagonal 2×2 matrices
diag_matrices = [[i 0; 0 i+1] for i in 1:6]
K64 = reduce(kron, diag_matrices)  # 64×64
println("K = D₁ ⊗ D₂ ⊗ ... ⊗ D₆ where Dᵢ = diag(i, i+1)")
println("64×64 matrix, but eigenvalues are products of 2×2 eigenvalues!")

vals_k64, _, _ = symbolic_eigenvalues(K64)
println("\nTotal eigenvalues: $(length(vals_k64))")
println("Unique eigenvalues: $(length(unique(vals_k64)))")
println("Min eigenvalue: $(minimum(vals_k64))")
println("Max eigenvalue: $(maximum(vals_k64))")
println("\nNote: Works for 10-fold (1024×1024) in ~33 seconds!")

# ============================================================================
# Example 16: Symbolic Rotation Kronecker (SO(2)^⊗2)
# ============================================================================

println("\n\n📊 Example 16: SO(2)⊗SO(2) (4×4 Symbolic)")
println("-" ^ 50)

@variables θ₁ θ₂
R_kron = kron(SO2_rotation(θ₁), SO2_rotation(θ₂))
println("Double Kronecker product of rotation matrices")
println("4×4 matrix with 2 symbolic angle parameters")

vals_R2, _, _ = symbolic_eigenvalues(R_kron)
println("\nEigenvalues (all combinations e^{i(±θ₁±θ₂)}):")
for (i, v) in enumerate(vals_R2)
    println("  λ$i = ", Symbolics.simplify(v))
end

# ============================================================================
# Summary: Why These Work
# ============================================================================

println("\n\n" * "=" ^ 70)
println("Why These 'Impossible' Matrices Work")
println("=" ^ 70)
println("""
Each example exploits algebraic structure that bypasses polynomial solving:

• Circulant (any size): DFT diagonalizes → eigenvalues from Fourier transform
• Hypercube Qₙ: Walsh-Hadamard basis → eigenvalues n-2k with binomial multiplicities
• Path/Cycle Laplacian: Chebyshev polynomials → trigonometric closed form
• Hadamard 2ⁿ: Self-similar structure → only 2 distinct eigenvalues ±√(2ⁿ)
• DFT matrix: F⁴ = nI → eigenvalues are 4th roots of n
• Cartan Aₙ: Root system theory → sine formula
• Kronecker products: λ(A⊗B) = λ(A)·λ(B) → reduce to smaller problems

The Abel-Ruffini barrier applies to *general* matrices. Structured matrices
have hidden symmetries that make symbolic solutions possible!
""")
