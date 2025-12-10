using SymbolicDiagonalization
using LinearAlgebra
using Symbolics
using Test

# ============================================================================
# Basic Symbolic Matrix Tests
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

# ============================================================================
# Structured Matrix Tests
# Tests for matrices with specific structure (triangular, Hermitian, block)
# ============================================================================

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
    mat = [a b 0; b' a 0; 0 0 c]  # Hermitian block with conjugate off-diagonal
    vals, poly, λ = symbolic_eigenvalues(mat; var = λ, structure = :hermitian)
    # Characteristic polynomial should factor as (λ - c)*((λ - a)^2 - b*b')
    expected = (λ - c) * ((λ - a)^2 - b * b')
    @test isequal(Symbolics.expand(poly - expected), 0)
    P, D, _ = symbolic_diagonalize(mat; var = λ, structure = :hermitian)
    diffmat = simplify.(mat * P - P * D)
    @test all(x -> Symbolics.iszero(Symbolics.expand(Symbolics.simplify(x))), diffmat)
    # expand=false path on block split
    vals2, poly2, _ = symbolic_eigenvalues(mat; var = λ, structure = :hermitian, expand = false)
    @test length(vals2) == 3
end

# ============================================================================
# 4×4 Matrix Tests
# Tests for 4×4 matrices (numeric, block-diagonal, diagonal)
# ============================================================================

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
    # Eigenvalues come from 2×2 blocks: each block has eigenvalues that differ by 2*off_diag
    # Block 1: [a b; b a] → eigenvalues a±b
    # Block 2: [c d; d c] → eigenvalues c±d
    # Since simplification may not reduce quadratic formula, check that we got 4 values
    # and that diagonalization works correctly
    @test length(vals) == 4
    # Test full diagonalization
    P, D, _ = symbolic_diagonalize(
        mat; 
        var = λ, 
        structure = :hermitian,
        expand = false,
        complexity_threshold = nothing
    )
    diffmat = simplify.(mat * P - P * D)
    @test all(x -> Symbolics.iszero(Symbolics.expand(Symbolics.simplify(x))), diffmat)
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
    # Note: use Matrix instead of Diagonal to avoid ambiguity with LinearAlgebra's eigen(::Diagonal)
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
    # Note: use Matrix instead of Diagonal to avoid ambiguity with LinearAlgebra's eigvals(::Diagonal)
    mat4 = Matrix(Diagonal([a, b, c]))
    vals4 = eigvals(mat4; expand = false, complexity_threshold = nothing)
    @test isequal(vals4, [a, b, c])
    
    # Test that non-diagonalizable matrices still return eigenvalues
    # (eigvals doesn't need diagonalization)
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
    # This is more of a design test - eigvals should use symbolic_eigenvalues
    # which skips the eigenvector computation
    @test length(vals) == 2
    @test all(v -> v isa Num, vals)
end

# ============================================================================
# Persymmetric Matrix Tests
# Tests for matrices that are symmetric about the anti-diagonal
# These matrices can be decomposed into smaller blocks for efficiency
# ============================================================================

@testset "Persymmetric symmetric 4x4" begin
    @variables a b c d λ
    # Matrix with special persymmetric structure that was problematic before optimization
    # This matrix would timeout before the persymmetric optimization was added
    Q = [a b 0 0; b c d 0; 0 d c b; 0 0 b a]
    
    # This should now work efficiently with the persymmetric optimization
    vals = eigvals(Q; var = λ, complexity_threshold = nothing)
    @test length(vals) == 4
    
    # Verify eigenvalues are symbolic expressions (not errors)
    for v in vals
        @test v isa Num
    end
    
    # Verify the matrix is correctly identified as symmetric and persymmetric
    @test SymbolicDiagonalization._is_symmetric(Q)
    @test SymbolicDiagonalization._is_persymmetric(Q)
    
    # Verify persymmetric split works
    split_result = SymbolicDiagonalization._persymmetric_split(Q)
    @test !isnothing(split_result)
    block1, block2, P = split_result
    @test size(block1) == (2, 2)
    @test size(block2) == (2, 2)
    
    # Verify transformation is correct: P^T * Q * P should be block diagonal
    Q_transformed = transpose(P) * Q * P
    Q_transformed = Symbolics.simplify.(Q_transformed)
    @test all(x -> isequal(Symbolics.simplify(x), 0), Q_transformed[1:2, 3:4])
    @test all(x -> isequal(Symbolics.simplify(x), 0), Q_transformed[3:4, 1:2])
end

@testset "Persymmetric Hermitian 4x4" begin
    @variables a::Real b_re::Real b_im::Real c::Real d_re::Real d_im::Real λ
    
    # Define complex elements
    b = b_re + im*b_im
    d = d_re + im*d_im
    
    # Hermitian persymmetric 4x4 with complex conjugate pairs in off-diagonals
    # This structure appears in quantum mechanics and signal processing
    Q = [a b 0 0; 
         conj(b) c d 0; 
         0 conj(d) c b; 
         0 0 conj(b) a]
    
    # Verify the matrix is correctly identified as Hermitian and persymmetric
    @test SymbolicDiagonalization._is_hermitian(Q)
    @test SymbolicDiagonalization._is_persymmetric(Q)
    
    # Verify persymmetric split works and produces real blocks
    split_result = SymbolicDiagonalization._persymmetric_split(Q)
    @test !isnothing(split_result)
    block1, block2, P = split_result
    @test size(block1) == (2, 2)
    @test size(block2) == (2, 2)
    
    # Blocks should be real (Hermitian persymmetric → real blocks)
    @test eltype(block1) == Num
    @test eltype(block2) == Num
    
    # Verify that imaginary part b_im doesn't appear in blocks (it cancels)
    for i in 1:2, j in 1:2
        vars_in_block1 = Symbolics.get_variables(block1[i, j])
        vars_in_block2 = Symbolics.get_variables(block2[i, j])
        @test !(b_im in vars_in_block1)
        @test !(b_im in vars_in_block2)
        @test !(d_im in vars_in_block1)
        @test !(d_im in vars_in_block2)
    end
    
    # Compute eigenvalues - should work efficiently
    vals = eigvals(Q; var = λ, complexity_threshold = nothing, timeout = 60)
    @test length(vals) == 4
    
    # Verify eigenvalues are real (Hermitian matrices have real eigenvalues)
    for v in vals
        @test v isa Num
        # Eigenvalues should only depend on the real parts
        vars_in_val = Symbolics.get_variables(v)
        @test !(b_im in vars_in_val)
        @test !(d_im in vars_in_val)
    end
end

# ============================================================================
# Special Pattern Tests
# Tests for matrices with special tridiagonal patterns that have closed-form eigenvalues
# Pattern [b,d,b,b] and [b,b,d,b] in 5×5 tridiagonal matrices
# ============================================================================

@testset "Special 5x5 tridiagonal pattern" begin
    @variables a b d λ
    
    # Special symmetric tridiagonal matrix with pattern [b, d, b, b]
    # This matrix has closed-form eigenvalues despite being 5×5
    M = [a  b  0  0  0;
         b  a  d  0  0;
         0  d  a  b  0;
         0  0  b  a  b;
         0  0  0  b  a]
    
    # Should successfully compute eigenvalues
    vals = eigvals(M; var = λ, complexity_threshold = nothing)
    @test length(vals) == 5
    
    # Expected eigenvalues: a - √(2b² + d²), a - b, a, a + b, a + √(2b² + d²)
    # Check that all eigenvalues are present
    @test any(v -> isequal(simplify(v - a), 0), vals)
    @test any(v -> isequal(simplify(v - (a - b)), 0), vals)
    @test any(v -> isequal(simplify(v - (a + b)), 0), vals)
    
    # Verify with numeric substitution
    vals_numeric = Float64[]
    for v in vals
        v_sub = substitute(v, Dict(a => 1, b => 2, d => 3))
        # Convert to expression and evaluate to get numeric value
        v_numeric = Float64(eval(Symbolics.toexpr(v_sub)))
        push!(vals_numeric, v_numeric)
    end
    
    # Expected values for a=1, b=2, d=3
    sqrt_term = sqrt(2*2^2 + 3^2)
    expected_numeric = sort([1 - sqrt_term, 1 - 2, 1, 1 + 2, 1 + sqrt_term])
    @test isapprox(sort(vals_numeric), expected_numeric, rtol=1e-10)
    
    # Test full eigendecomposition
    E = eigen(M; var = λ, complexity_threshold = nothing)
    @test length(E.values) == 5
    @test size(E.vectors) == (5, 5)
    
    # Verify eigenvalue-eigenvector property numerically
    M_numeric = Float64[1 2 0 0 0;
                        2 1 3 0 0;
                        0 3 1 2 0;
                        0 0 2 1 2;
                        0 0 0 2 1]
    E_numeric = LinearAlgebra.eigen(M_numeric)
    
    # Check that symbolic eigenvalues match numeric ones
    @test isapprox(sort(vals_numeric), sort(E_numeric.values), rtol=1e-10)
end

@testset "Special 5×5 Tridiagonal Pattern [b,b,d,b]" begin
    @variables a b d λ
    
    # Alternative pattern with d at position 3: [b, b, d, b]
    # This should give the same eigenvalues as [b, d, b, b]
    M = [a  b  0  0  0;
         b  a  b  0  0;
         0  b  a  d  0;
         0  0  d  a  b;
         0  0  0  b  a]
    
    # Should successfully compute eigenvalues
    vals = eigvals(M; var = λ, complexity_threshold = nothing)
    @test length(vals) == 5
    
    # Expected eigenvalues: same as [b,d,b,b] pattern
    @test any(v -> isequal(simplify(v - a), 0), vals)
    @test any(v -> isequal(simplify(v - (a - b)), 0), vals)
    @test any(v -> isequal(simplify(v - (a + b)), 0), vals)
    
    # Verify with numeric substitution
    vals_numeric = Float64[]
    for v in vals
        v_sub = substitute(v, Dict(a => 1, b => 2, d => 3))
        v_numeric = Float64(eval(Symbolics.toexpr(v_sub)))
        push!(vals_numeric, v_numeric)
    end
    
    # Expected values for a=1, b=2, d=3 (same as [b,d,b,b])
    sqrt_term = sqrt(2*2^2 + 3^2)
    expected_numeric = sort([1 - sqrt_term, 1 - 2, 1, 1 + 2, 1 + sqrt_term])
    @test isapprox(sort(vals_numeric), expected_numeric, rtol=1e-10)
end

# Note: Error handling tests are lightweight but included separately
# to avoid slowing down the main test suite
# include("test_error_handling.jl")
