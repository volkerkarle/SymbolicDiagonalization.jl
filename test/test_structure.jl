# ============================================================================
# Structure Detection Tests
# Tests for matrices with specific structure (Hermitian, block-diagonal, persymmetric)
# ============================================================================

@testset "Hermitian symbolic 2x2" begin
    @variables a c b λ
    mat = [a b; b' c]  # Hermitian: off-diagonals are conjugate pairs
    poly, coeffs, λ = characteristic_polynomial(mat; var = λ)
    @test isequal(simplify(coeffs[1] - (a*c - b*b')), 0)
    @test isequal(simplify(coeffs[2] + a + c), 0)
    P, D, _ = symbolic_diagonalize(mat; var = λ, structure = :hermitian)
    diffmat = simplify.(mat * P - P * D)
    @test all(x -> Symbolics.iszero(Symbolics.expand(Symbolics.simplify(x))), diffmat)
    # Test expand=false path does not error and returns an unexpanded polynomial
    vals2, poly2, _ = symbolic_eigenvalues(mat; var = λ, structure = :hermitian, expand = false)
    @test length(vals2) == 2
end

@testset "Hermitian block 3x3" begin
    @variables a c b λ
    mat = [a b 0; b a 0; 0 0 c]  # Hermitian block (real symmetric)
    vals, poly, λ = symbolic_eigenvalues(mat; var = λ, structure = :hermitian)
    # Should get 3 eigenvalues from the 3x3 matrix
    @test length(vals) == 3
    # expand=false path on block split
    vals2, poly2, _ = symbolic_eigenvalues(mat; var = λ, structure = :hermitian, expand = false)
    @test length(vals2) == 3
end

@testset "Block-diagonal Hermitian 4x4" begin
    @variables a b c d λ
    # 2×2 block + 2×2 block structure (should split automatically)
    mat = [
        a  b  0  0;
        b  a  0  0;
        0  0  c  d;
        0  0  d  c
    ]
    # This should succeed because _block_split detects two 2×2 blocks
    vals, poly, λ = symbolic_eigenvalues(
        mat; 
        var = λ, 
        structure = :hermitian,
        expand = false,
        complexity_threshold = nothing
    )
    @test length(vals) == 4
    # Test that diagonalization doesn't error
    P, D, _ = symbolic_diagonalize(
        mat; 
        var = λ, 
        structure = :hermitian,
        expand = false,
        complexity_threshold = nothing
    )
    # Just verify dimensions are correct
    @test size(P) == (4, 4)
    @test size(D) == (4, 4)
end

@testset "Numeric 4x4 eigenvalues" begin
    # 4×4 with numeric entries should work (produces large but valid expressions)
    mat = [
        1.0  2.0  0.0  0.0;
        2.0  1.0  0.0  0.0;
        0.0  0.0  3.0  1.0;
        0.0  0.0  1.0  3.0
    ]
    @variables λ
    vals, poly, _ = symbolic_eigenvalues(
        mat; 
        var = λ,
        complexity_threshold = nothing
    )
    @test length(vals) == 4
    # The eigenvalues should be numeric (no symbolic variables)
    for v in vals
        @test !isa(v, Symbolics.Num) || isempty(Symbolics.get_variables(v))
    end
end

@testset "Persymmetric matrices" begin
    @variables a b c d e f g h i λ
    
    # 3×3 persymmetric (symmetric about anti-diagonal)
    P3 = [a b c; b a b; c b a]
    @test SymbolicDiagonalization._is_persymmetric(P3)
    
    # Non-persymmetric matrix
    M = [a b c; d e f; g h i]
    @test !SymbolicDiagonalization._is_persymmetric(M)
    
    # 4×4 persymmetric
    P4 = [1 2 3 4; 2 5 6 3; 3 6 5 2; 4 3 2 1]
    @test SymbolicDiagonalization._is_persymmetric(P4)
end

@testset "Block splitting detection" begin
    @variables a b c d e f λ
    
    # Block diagonal matrix that should be detected
    mat = [
        a  0  0  0;
        0  b  0  0;
        0  0  c  d;
        0  0  e  f
    ]
    
    # Test that block splitting helps with eigenvalue computation
    vals, _, _ = symbolic_eigenvalues(mat; var = λ)
    @test length(vals) == 4
    
    # First two eigenvalues should be a, b directly
    # Last two come from 2×2 block [c d; e f]
    eigenvals_set = Set(simplify.(vals))
    @test a in eigenvals_set || simplify(vals[1] - a) == 0 || simplify(vals[2] - a) == 0
    @test b in eigenvals_set || simplify(vals[1] - b) == 0 || simplify(vals[2] - b) == 0
end
