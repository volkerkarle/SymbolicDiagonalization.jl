using SymbolicDiagonalization
using Symbolics
using LinearAlgebra
using Test

@testset "Error handling and warnings" begin
    # Test 1: Verify exported exception types
    @testset "Exception types exported" begin
        @test ExpressionComplexityError <: Exception
        @test ComputationTimeoutError <: Exception
        
        # Test constructors
        e1 = ExpressionComplexityError("test message")
        @test e1.message == "test message"
        
        e2 = ComputationTimeoutError("timeout test")
        @test e2.message == "timeout test"
    end
    
    # Test 2: Complexity warning system
    @testset "Complexity warnings" begin
        @variables a b c d e f
        A = [a b c; d e f; a b c]
        
        # Should warn with low threshold (6 variables > 3 threshold)
        @test_logs (:warn, r"Matrix contains 6 symbolic variables") eigvals(A; complexity_threshold=3)
        
        # Should not warn when disabled
        @test_logs eigvals([a 0; 0 b]; complexity_threshold=nothing)
    end
    
    # Test 3: Parameter propagation through API
    @testset "Parameter propagation" begin
        @variables a b
        A = [a 0; 0 b]
        
        # All these should work and return the same result
        vals1 = eigvals(A)
        vals2 = eigvals(A; max_terms=5000, timeout=60)
        vals3 = symbolic_eigenvalues(A; max_terms=5000, timeout=60)[1]
        
        @test Set(vals1) == Set(vals2)
        @test Set(vals1) == Set(vals3)
    end
end
