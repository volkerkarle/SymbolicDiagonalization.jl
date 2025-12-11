# ============================================================================
# LinearAlgebra Interface Tests
# Tests for eigen() and eigvals() functions matching Julia's standard interface
# ============================================================================

@testset "LinearAlgebra.eigen interface" begin
    @variables a b λ
    
    # Test basic diagonal matrix
    mat = [a 0; 0 b]
    E = eigen(mat; var = λ)
    @test isa(E, LinearAlgebra.Eigen)
    @test length(E.values) == 2
    @test size(E.vectors) == (2, 2)
    @test Set(simplify.(E.values)) == Set([a, b])
    
    # Verify eigenvector property: A*v = λ*v for each column
    for i in 1:2
        v = E.vectors[:, i]
        λ_i = E.values[i]
        result = simplify.(mat * v - λ_i * v)
        @test all(x -> isequal(x, 0), result)
    end
    
    # Test with triangular matrix
    @variables c
    mat2 = [a 1 0; 0 b 1; 0 0 c]
    E2 = eigen(mat2; var = λ)
    @test Set(simplify.(E2.values)) == Set([a, b, c])
    @test size(E2.vectors) == (3, 3)
    
    # Test keyword arguments are passed through
    mat3 = Matrix([a 0; 0 b])
    E3 = eigen(mat3; structure = :diagonal, expand = false, complexity_threshold = nothing)
    @test length(E3.values) == 2
    
    # Test that non-diagonalizable matrix throws error
    @variables x
    jordan = [x 1; 0 x]
    @test_throws ErrorException eigen(jordan; var = λ)
end

@testset "LinearAlgebra.eigvals interface" begin
    @variables a b c λ
    
    # Test basic diagonal matrix
    mat = [a 0; 0 b]
    vals = eigvals(mat; var = λ)
    @test isa(vals, Vector)
    @test length(vals) == 2
    @test Set(simplify.(vals)) == Set([a, b])
    
    # Test with triangular matrix
    mat2 = [a 1 0; 0 b 1; 0 0 c]
    vals2 = eigvals(mat2; var = λ)
    @test length(vals2) == 3
    @test Set(simplify.(vals2)) == Set([a, b, c])
    
    # Test that eigenvalues match symbolic_eigenvalues
    vals_ref, _, _ = symbolic_eigenvalues(mat2; var = λ)
    @test Set(simplify.(vals2)) == Set(simplify.(vals_ref))
    
    # Test with Hermitian matrix
    @variables d
    mat3 = [a b; b' d]
    vals3 = eigvals(mat3; var = λ, structure = :hermitian)
    @test length(vals3) == 2
    
    # Test keyword arguments
    mat4 = Matrix(Diagonal([a, b, c]))
    vals4 = eigvals(mat4; expand = false, complexity_threshold = nothing)
    @test isequal(vals4, [a, b, c])
    
    # Test that non-diagonalizable matrices still return eigenvalues
    @variables x
    jordan = [x 1; 0 x]
    vals5 = eigvals(jordan; var = λ)
    @test length(vals5) == 2
    # Both eigenvalues should be x (repeated root)
    @test all(v -> isequal(simplify(v - x), 0), vals5)
end

@testset "LinearAlgebra interface comparison" begin
    @variables a b λ
    mat = [a 0; 0 b]
    
    # Test that eigen and eigvals are consistent
    E = eigen(mat; var = λ)
    vals = eigvals(mat; var = λ)
    
    @test Set(simplify.(E.values)) == Set(simplify.(vals))
    
    # Test that eigvals is cheaper (no vector computation)
    @test length(vals) == 2
    @test all(v -> v isa Num, vals)
end
