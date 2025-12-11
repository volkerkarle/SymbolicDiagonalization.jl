# ============================================================================
# Basic Functionality Tests
# Tests for simple matrices (diagonal, Jordan blocks) and basic functionality
# ============================================================================

@testset "Symbolic diagonal entries" begin
    @variables a b λ
    mat = [a 0; 0 b]
    poly, coeffs, _ = characteristic_polynomial(mat; var = λ)
    @test isequal(simplify(coeffs[1] - a * b), 0)
    @test isequal(simplify(coeffs[2] + a + b), 0)
    pairs, _, _ = symbolic_eigenpairs(mat; var = λ)
    vals = simplify.(first.(pairs))
    @test Set(vals) == Set([a, b])
    P, D, _ = symbolic_diagonalize(mat; var = λ)
    diffmat = simplify.(mat * P - P * D)
    @test all(x -> isequal(x, 0), diffmat)
end

@testset "Jordan-like matrix" begin
    @variables a λ
    mat = [a 1; 0 a]
    poly, coeffs, _ = characteristic_polynomial(mat; var = λ)
    @test isequal(simplify(coeffs[1] - a^2), 0)
    @test isequal(simplify(coeffs[2] + 2a), 0)
    vals = symbolic_roots(coeffs)
    @test isequal(Symbolics.expand(sum(vals) - 2a), 0)
    prod_simplified = Symbolics.expand(prod(vals) - a^2)
    @test isequal(prod_simplified, 0) || isequal(prod_simplified, -0.0)
    @test_throws ErrorException symbolic_diagonalize(mat; var = λ)
end

@testset "Symbolic quartic computation" begin
    @variables c0 c1 c2 c3 c4
    # Quartic roots now work symbolically (though they produce large expressions)
    roots = symbolic_roots([c0, c1, c2, c3, c4])
    @test length(roots) == 4
    # Diagonal matrices should still give direct eigenvalues
    vals, _, _ = symbolic_eigenvalues(Diagonal([c0, c1, c2, c3]); var = c4)
    @test isequal(vals, [c0, c1, c2, c3])
end

@testset "3x3 symbolic triangular" begin
    @variables a b c λ
    mat = [a 1 0; 0 b 1; 0 0 c]
    pairs, poly, _ = symbolic_eigenpairs(mat; var = λ)
    vals = simplify.(first.(pairs))
    @test Set(vals) == Set([a, b, c])
    @test isequal(Symbolics.expand(poly - (λ - a)*(λ - b)*(λ - c)), 0)
    # expand=false path does not change eigenvalues
    vals2, poly2, _ = symbolic_eigenvalues(mat; var = λ, expand = false)
    @test Set(vals2) == Set([a, b, c])
end

@testset "Diagonal 4x4" begin
    @variables a b c d λ
    # Diagonal matrices should work perfectly regardless of size
    mat = Diagonal([a, b, c, d])
    vals, poly, _ = symbolic_eigenvalues(mat; var = λ)
    @test isequal(vals, [a, b, c, d])
    # Full diagonalization should also work
    P, D, _ = symbolic_diagonalize(mat; var = λ, complexity_threshold = nothing)
    @test isequal(D, mat)
end

@testset "Simple 2×2 symbolic matrices" begin
    @variables a b c d λ
    
    # Test 2×2 matrix with 4 symbolic entries
    mat = [a b; c d]
    vals, poly, _ = symbolic_eigenvalues(mat; var = λ)
    @test length(vals) == 2
    
    # Verify using Vieta's formulas
    # Sum of roots = trace = a + d
    sum_vals = simplify(sum(vals))
    @test isequal(simplify(sum_vals - (a + d)), 0)
    
    # Product of roots = determinant = ad - bc
    prod_vals = Symbolics.expand(simplify(prod(vals)))
    det_val = a*d - b*c
    @test isequal(Symbolics.expand(prod_vals - det_val), 0)
end

# NOTE: Skipping general 3×3 symbolic matrix test - too many variables (9) makes it very slow
# For 3×3 symbolic testing, use triangular matrices or matrices with fewer unique variables
