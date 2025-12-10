using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

# Analyze the structure of the user's matrix more carefully
@variables a b c d λ

Q = [a b 0 0; b c d 0; 0 d c b; 0 0 b a]

println("=== Matrix Structure Analysis ===\n")
println("Original matrix Q:")
display(Q)
println("\n")

# Check various structural properties
is_sym = all([isequal(Q[i,j], Q[j,i]) for i in 1:4, j in 1:4])
println("1. Symmetry: Q == Q' ? ", is_sym)

# Check persymmetry: Q[i,j] == Q[n+1-j, n+1-i]
n = 4
J = [0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]  # Exchange matrix
JQJ = J * transpose(Q) * J
is_persym = all([isequal(Q[i,j], JQJ[i,j]) for i in 1:4, j in 1:4])
println("2. Persymmetric (Q == J*Q'*J): ", is_persym)

# Check if it commutes with J
QJ = Q * J
JQ = J * Q
commutes = all([isequal(QJ[i,j], JQ[i,j]) for i in 1:4, j in 1:4])
println("3. Commutes with J (Q*J == J*Q): ", commutes)

# The characteristic polynomial - let's look at it carefully
println("\n=== Characteristic Polynomial ===")
poly, coeffs, _ = characteristic_polynomial(Q; var = λ)
println("Coefficients (constant to highest degree):")
for (i, c) in enumerate(coeffs)
    println("  λ^$(i-1): ", c)
end

println("\n=== Key Observation ===")
println("Notice the block structure:")
println("  [ A    B  ]     where A = [a b]    B = [0 0]")
println("  [ B^T  C  ]           [b c]        [d 0]")
println("                                     ")
println("                        C = [c b]    B^T = [0 d]")
println("                            [b a]          [0 0]")
println()
println("This is NOT block-diagonal, but B is very sparse!")
println("In fact, B has only ONE non-zero entry: B[2,1] = d")

# Let's try a manual simplification approach
println("\n=== Manual Approach ===")
println("For a matrix [A B; B^T C], if B is 'small', we might be able to")
println("use perturbation theory or a specialized decomposition.")
println()
println("Alternatively, we notice:")
println("  - Upper-left 2x2:  [a b; b c]")  
println("  - Lower-right 2x2: [c b; b a]")
println("  - Connection via single off-diagonal entry d")

# Check if the 2x2 blocks would be easy to diagonalize separately
A_block = Q[1:2, 1:2]
C_block = Q[3:4, 3:4]

println("\nUpper-left block eigenvalues (if it were isolated):")
println("  det([a-λ, b; b, c-λ]) = 0")
println("  (a-λ)(c-λ) - b² = 0")
println("  λ² - (a+c)λ + (ac - b²) = 0")
println("  λ = ((a+c) ± sqrt((a+c)² - 4(ac-b²)))/2")
println("    = ((a+c) ± sqrt(a² + 2ac + c² - 4ac + 4b²))/2")
println("    = ((a+c) ± sqrt(a² - 2ac + c² + 4b²))/2")
println("    = ((a+c) ± sqrt((a-c)² + 4b²))/2")

println("\nLower-right block eigenvalues (if it were isolated):")
println("  det([c-λ, b; b, a-λ]) = 0")  
println("  Same as above! λ = ((a+c) ± sqrt((a-c)² + 4b²))/2")

println("\n*** KEY INSIGHT: Both 2x2 blocks have the SAME eigenvalues! ***")
println("This suggests the full 4x4 matrix might have special structure")
println("that makes its eigenvalues related to these 2x2 blocks.")
