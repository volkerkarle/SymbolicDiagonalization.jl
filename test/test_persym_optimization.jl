using SymbolicDiagonalization
using Symbolics
using LinearAlgebra
using Test

@testset "Persymmetric 4x4 optimization" begin
    @variables a b c d λ
    
    # The user's problematic matrix
    Q = [a b 0 0; b c d 0; 0 d c b; 0 0 b a]
    
    println("\n=== Testing Persymmetric Split ===")
    println("Original matrix Q:")
    display(Q)
    println()
    
    # Test the persymmetric split function directly
    persym_result = SymbolicDiagonalization._persymmetric_split(Q)
    
    if isnothing(persym_result)
        println("ERROR: Persymmetric split returned nothing!")
    else
        block1, block2, P = persym_result
        println("SUCCESS: Matrix split into two 2×2 blocks")
        println("\nBlock 1 (symmetric part):")
        display(block1)
        println("\nBlock 2 (antisymmetric part):")
        display(block2)
        println("\nTransformation matrix P:")
        display(P)
        
        # Verify the transformation
        Q_transformed = transpose(P) * Q * P
        println("\n\nVerifying: P^T * Q * P =")
        display(Symbolics.simplify.(Q_transformed))
        
        # Check if it's block diagonal
        println("\n\nIs P^T*Q*P block diagonal?")
        upper_right = Q_transformed[1:2, 3:4]
        lower_left = Q_transformed[3:4, 1:2]
        println("Upper-right block is zero: ", all(x -> isequal(Symbolics.simplify(x), 0), upper_right))
        println("Lower-left block is zero: ", all(x -> isequal(Symbolics.simplify(x), 0), lower_left))
    end
    
    println("\n=== Testing Full Eigenvalue Computation ===")
    println("Attempting to compute eigenvalues with persymmetric optimization...")
    
    try
        @time vals = eigvals(Q; var = λ, timeout = 30, max_terms = 10000, complexity_threshold = nothing)
        println("\n✓ SUCCESS! Eigenvalues computed:")
        println("Number of eigenvalues: ", length(vals))
        for (i, v) in enumerate(vals)
            println("  λ$i = ", v)
        end
    catch e
        println("\n✗ FAILED")
        println("Error type: ", typeof(e))
        println("Error: ", e)
    end
end
