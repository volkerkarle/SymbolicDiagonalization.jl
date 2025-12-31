# ============================================================================
# Edge Cases and Error Handling Tests
# Tests for error conditions, boundary cases, and exception handling
# ============================================================================

@testset "Input Validation" begin
    @testset "Non-square matrix rejection" begin
        # 2x3 matrix should be rejected
        mat = [1 2 3; 4 5 6]
        @test_throws ArgumentError symbolic_eigenvalues(mat)
        @test_throws ArgumentError symbolic_eigenpairs(mat)
        @test_throws ArgumentError symbolic_diagonalize(mat)
        @test_throws ArgumentError characteristic_polynomial(mat)
    end
    
    @testset "1x1 matrix" begin
        @variables a λ
        mat = fill(a, 1, 1)
        vals, poly, _ = symbolic_eigenvalues(mat; var=λ)
        @test length(vals) == 1
        @test isequal(vals[1], a)
        
        # Numeric 1x1
        mat_num = [42.0;;]
        vals_num, _, _ = symbolic_eigenvalues(mat_num)
        @test vals_num == [42.0]
    end
    
    @testset "Invalid structure hint" begin
        mat = [1 0; 0 1]
        @test_throws ArgumentError symbolic_eigenvalues(mat; structure=:invalid)
        @test_throws ArgumentError symbolic_eigenpairs(mat; structure=:bogus)
    end
end

@testset "Exception Types" begin
    @testset "ComputationTimeoutError with very short timeout" begin
        @variables a b c d
        # Create a matrix that requires computation
        mat = [a b; c d]
        
        # With an extremely short timeout (0.001 seconds), this should timeout
        # Note: This test may be flaky on very fast systems
        @test_throws ComputationTimeoutError symbolic_eigenvalues(mat; timeout=0.001)
    end
    
    @testset "ExpressionComplexityError direct test" begin
        # Test exception type directly since triggering it through normal
        # computation is difficult (requires very large expressions)
        err = ExpressionComplexityError("Expression too complex")
        @test err isa Exception
        @test occursin("complex", lowercase(err.message))
    end
    
    @testset "Exception messages are informative" begin
        # Test that exception messages contain helpful information
        try
            throw(ComputationTimeoutError("test timeout message"))
        catch e
            msg = sprint(showerror, e)
            @test occursin("timeout", lowercase(msg))
        end
        
        try
            throw(ExpressionComplexityError("test complexity message"))
        catch e
            msg = sprint(showerror, e)
            @test occursin("complexity", lowercase(msg))
        end
    end
end

@testset "Numeric Edge Cases" begin
    @testset "Zero matrix" begin
        mat = zeros(3, 3)
        vals, _, _ = symbolic_eigenvalues(mat)
        @test all(v -> v == 0 || isequal(v, 0), vals)
    end
    
    @testset "Identity matrix" begin
        for n in 1:4
            mat = Matrix(1.0I, n, n)
            vals, _, _ = symbolic_eigenvalues(mat)
            @test all(v -> v == 1 || isequal(v, 1), vals)
        end
    end
    
    @testset "Matrix with repeated eigenvalues" begin
        # Diagonal with repeated values
        mat = Diagonal([2.0, 2.0, 3.0])
        vals, _, _ = symbolic_eigenvalues(mat)
        @test count(v -> v == 2.0, vals) == 2
        @test count(v -> v == 3.0, vals) == 1
    end
    
    @testset "Complex numeric matrix" begin
        mat = [1.0+im 0; 0 2.0-im]
        vals, _, _ = symbolic_eigenvalues(mat)
        @test Set(vals) == Set([1.0+im, 2.0-im])
    end
end

@testset "Symbolic Edge Cases" begin
    @testset "Matrix with zero symbolic variable" begin
        @variables a
        # Diagonal matrix with a and 0
        mat = Diagonal([a, 0])
        vals, _, _ = symbolic_eigenvalues(mat)
        @test length(vals) == 2
        # Eigenvalues should be a and 0 - check using isequal which works with Num
        @test any(v -> isequal(v, a), vals)
        # For zero check, use a try-catch to handle symbolic comparison
        has_zero = false
        for v in vals
            try
                if v == 0 || isequal(v, 0)
                    has_zero = true
                    break
                end
            catch
                # Symbolic comparison failed, try unwrap
                if isequal(Symbolics.unwrap(v), 0)
                    has_zero = true
                    break
                end
            end
        end
        @test has_zero
    end
    
    @testset "Scalar multiple of identity" begin
        @variables k
        mat = k * Matrix(I, 3, 3)
        vals, _, _ = symbolic_eigenvalues(mat)
        @test length(vals) == 3
        @test all(v -> isequal(simplify(v - k), 0), vals)
    end
end

@testset "Timeout disabled" begin
    @variables a b
    mat = [a b; b a]  # Simple symmetric matrix
    # With timeout=nothing, computation should complete (no timeout check)
    vals, _, _ = symbolic_eigenvalues(mat; timeout=nothing)
    @test length(vals) == 2
end

@testset "Complexity threshold disabled" begin
    @variables a b c d e f
    # 2x2 with many variables (6 total, above default threshold of 5)
    # This would normally warn, but with threshold=nothing it shouldn't
    mat = [a+b+c d; e f]
    # Should complete without warning when threshold is disabled
    vals, _, _ = symbolic_eigenvalues(mat; complexity_threshold=nothing)
    @test length(vals) == 2
end
