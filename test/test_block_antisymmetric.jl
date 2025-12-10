using SymbolicDiagonalization
using Symbolics
using LinearAlgebra
using Test

# Test the problematic 4x4 block-antisymmetric matrix from the user
@testset "User's 4x4 block-antisymmetric matrix" begin
    @variables a b c d λ
    
    # The user's matrix: [a b 0 0; b c d 0; 0 d c b; 0 0 b a]
    Q = [a b 0 0; b c d 0; 0 d c b; 0 0 b a]
    
    println("\n=== Testing User's Matrix ===")
    println("Matrix structure:")
    display(Q)
    println()
    
    # First, let's just try the characteristic polynomial
    println("\n--- Step 1: Computing characteristic polynomial ---")
    try
        @time poly, coeffs, λ_var = characteristic_polynomial(Q; var = λ)
        println("SUCCESS: Characteristic polynomial computed")
        println("Number of coefficients: ", length(coeffs))
        println("Polynomial degree: ", length(coeffs) - 1)
        # Don't print the full polynomial as it might be huge
    catch e
        println("FAILED at characteristic polynomial stage")
        println("Error: ", e)
        rethrow(e)
    end
    
    # Next, try to find eigenvalues with a short timeout
    println("\n--- Step 2: Finding roots (eigenvalues) ---")
    println("Note: This is where it likely fails or takes too long")
    try
        @time vals = eigvals(Q; var = λ, timeout = 30, max_terms = 5000, complexity_threshold = nothing)
        println("SUCCESS: Eigenvalues computed")
        println("Number of eigenvalues: ", length(vals))
        # Don't print the eigenvalues as they might be huge
    catch e
        println("FAILED at root-finding stage")
        println("Error type: ", typeof(e))
        println("Error message: ", e)
        # Don't rethrow - we expect this to fail
    end
    
    println("\n=== Analysis ===")
    println("Matrix properties:")
    println("  - Size: 4x4")
    println("  - Symmetric: ", Q == transpose(Q))
    println("  - Number of distinct variables: 4 (a, b, c, d)")
    println("  - Special structure: block-antisymmetric pattern")
    println("    Upper-left 2x2:  [a b; b c]")
    println("    Lower-right 2x2: [c b; b a]")
    println("    Notice: lower-right has 'flipped' structure")
    println()
    println("The matrix has a special 'persymmetric' or 'anti-diagonal symmetry':")
    println("  Q[i,j] relates to Q[n+1-i, n+1-j] in a specific pattern")
end

# Let's also test if we can identify this pattern
@testset "Detect persymmetric/anti-diagonal structure" begin
    @variables a b c d
    Q = [a b 0 0; b c d 0; 0 d c b; 0 0 b a]
    
    println("\n=== Checking for special structures ===")
    
    # Check if it's symmetric
    is_sym = Q == transpose(Q)
    println("Symmetric: ", is_sym)
    
    # Check if it's persymmetric: Q[i,j] == Q[n+1-j, n+1-i]
    n = size(Q, 1)
    is_persym = true
    for i in 1:n, j in 1:n
        if !isequal(Q[i, j], Q[n+1-j, n+1-i])
            is_persym = false
            break
        end
    end
    println("Persymmetric: ", is_persym)
    
    # Check if upper-right and lower-left blocks are zero
    upper_right_zero = all(iszero, Q[1:2, 3:4])
    lower_left_zero = all(iszero, Q[3:4, 1:2])
    println("Upper-right 2x2 block is zero: ", upper_right_zero)
    println("Lower-left 2x2 block is zero: ", lower_left_zero)
    
    # This suggests an "arrow" or "cross" pattern
    if upper_right_zero && lower_left_zero
        println("\nMatrix has block-cross structure!")
        println("Upper-left block:  ", Q[1:2, 1:2])
        println("Upper-right block: ", Q[1:2, 3:4])
        println("Lower-left block:  ", Q[3:4, 1:2])
        println("Lower-right block: ", Q[3:4, 3:4])
    end
end
