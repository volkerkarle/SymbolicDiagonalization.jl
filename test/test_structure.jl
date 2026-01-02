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

# ============================================================================
# Complex Symbolic Hermitian Matrix Tests
# Tests for Hermitian matrices with complex symbolic entries (br + bi*im)
# ============================================================================

# Helper function to evaluate substituted symbolic expression to Float64
function _eval_substituted(v, subs_dict)
    substituted = substitute(v, subs_dict)
    if substituted isa Symbolics.Num
        expr = Symbolics.toexpr(substituted)
        return Float64(real(eval(expr)))
    else
        return Float64(real(substituted))
    end
end

@testset "Complex Symbolic Hermitian 2×2" begin
    @variables d1 d2 br bi
    
    # Hermitian matrix with complex off-diagonal: [d1, br+bi*im; br-bi*im, d2]
    b_complex = br + bi*im
    H = [d1 b_complex; conj(b_complex) d2]
    
    # Compute symbolic eigenvalues
    vals, poly, λ = symbolic_eigenvalues(H)
    @test length(vals) == 2
    
    # Verify symbolic expressions are correct by substituting values
    subs_dict = Dict(d1 => 3.0, d2 => 1.0, br => 1.0, bi => 0.5)
    H_numeric = [3.0 1.0+0.5im; 1.0-0.5im 1.0]
    true_eigs = eigvals(Hermitian(H_numeric))
    
    computed_eigs = [_eval_substituted(v, subs_dict) for v in vals]
    @test sort(true_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Complex Symbolic Hermitian 2×2 - Multiple parameter sets" begin
    @variables d1 d2 br bi
    
    b_complex = br + bi*im
    H = [d1 b_complex; conj(b_complex) d2]
    
    # Get symbolic eigenvalues once
    vals, poly, λ = symbolic_eigenvalues(H)
    @test length(vals) == 2
    
    # Verify the symbolic expressions work for multiple parameter sets
    test_cases = [
        (d1 = 1.0, d2 = 1.0, br = 0.0, bi = 0.0),  # Diagonal
        (d1 = 2.0, d2 = 0.0, br = 1.0, bi = 0.0),  # Real off-diagonal
        (d1 = 0.0, d2 = 0.0, br = 0.0, bi = 1.0),  # Pure imaginary off-diagonal
        (d1 = 5.0, d2 = -3.0, br = 2.0, bi = 1.5), # General case
    ]
    
    for case in test_cases
        subs_dict = Dict(d1 => case.d1, d2 => case.d2, br => case.br, bi => case.bi)
        H_numeric = [case.d1 case.br+case.bi*im; case.br-case.bi*im case.d2]
        true_eigs = eigvals(Hermitian(H_numeric))
        
        computed_eigs = [_eval_substituted(v, subs_dict) for v in vals]
        @test sort(true_eigs) ≈ sort(computed_eigs) atol=1e-10
    end
end

@testset "Symbolic Hermitian Kronecker Product 4×4 = 2×2 ⊗ 2×2" begin
    @variables a b c d
    
    # Two symbolic symmetric 2×2 matrices
    H1 = [a b; b c]
    H2 = [c d; d a]
    K = kron(H1, H2)
    
    # Should detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(K)
    @test !isnothing(kron_info)
    
    # Get symbolic eigenvalues
    vals, poly, λ = symbolic_eigenvalues(K)
    @test length(vals) == 4
    
    # Verify symbolic expressions by substitution
    subs_dict = Dict(a => 2.0, b => 0.5, c => 1.0, d => 0.3)
    K_numeric = kron([2.0 0.5; 0.5 1.0], [1.0 0.3; 0.3 2.0])
    true_eigs = eigvals(K_numeric)
    
    computed_eigs = [_eval_substituted(v, subs_dict) for v in vals]
    @test sort(true_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Symbolic Hermitian Kronecker Product 6×6 = 3×3 ⊗ 2×2" begin
    @variables a b c d e
    
    # 3×3 symmetric tridiagonal (keeps complexity manageable)
    H1 = [a b 0; b c b; 0 b a]
    # 2×2 symmetric
    H2 = [d e; e d]
    K = kron(H1, H2)  # 6×6
    
    # Should detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(K)
    @test !isnothing(kron_info)
    
    # Get symbolic eigenvalues
    vals, poly, λ = symbolic_eigenvalues(K)
    @test length(vals) == 6
    
    # Helper for Complex{Num} evaluation
    function _eval_complex_substituted(v, subs_dict)
        r = real(substitute(v, subs_dict))
        i = imag(substitute(v, subs_dict))
        r_val = r isa Symbolics.Num ? Float64(eval(Symbolics.toexpr(r))) : Float64(r)
        i_val = i isa Symbolics.Num ? Float64(eval(Symbolics.toexpr(i))) : Float64(i)
        return Complex(r_val, i_val)
    end
    
    # Verify symbolic expressions by substitution
    subs_dict = Dict(a => 2.0, b => 0.5, c => 1.0, d => 3.0, e => 0.3)
    H1_num = [2.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 2.0]
    H2_num = [3.0 0.3; 0.3 3.0]
    K_numeric = kron(H1_num, H2_num)
    true_eigs = eigvals(K_numeric)
    
    computed_eigs = [_eval_complex_substituted(v, subs_dict) for v in vals]
    # For Hermitian matrices, imaginary parts should be negligible
    @test all(abs.(imag.(computed_eigs)) .< 1e-10)
    @test sort(true_eigs) ≈ sort(real.(computed_eigs)) atol=1e-10
end

@testset "Symbolic Kronecker 6×6 = 3×3 ⊗ 2×2 (8 variables)" begin
    # Test with more symbolic variables: 5 in 3×3 + 3 in 2×2 = 8 total
    @variables a b c d e   # 3×3 symmetric
    @variables x y z       # 2×2 symmetric
    
    H1 = [a b c; b d e; c e a]  # Full symmetric 3×3 (5 vars)
    H2 = [x y; y z]              # Full symmetric 2×2 (3 vars)
    K = kron(H1, H2)  # 6×6
    
    # Should detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(K)
    @test !isnothing(kron_info)
    
    # Get symbolic eigenvalues
    vals, poly, λ = symbolic_eigenvalues(K)
    @test length(vals) == 6
    
    # Helper for Complex{Num} evaluation
    function _eval_complex_sub(v, subs_dict)
        r = real(substitute(v, subs_dict))
        i = imag(substitute(v, subs_dict))
        r_val = r isa Symbolics.Num ? Float64(eval(Symbolics.toexpr(r))) : Float64(r)
        i_val = i isa Symbolics.Num ? Float64(eval(Symbolics.toexpr(i))) : Float64(i)
        return Complex(r_val, i_val)
    end
    
    # Verify symbolic expressions by substitution
    subs_dict = Dict(a => 2.0, b => 0.3, c => 0.1, d => 1.5, e => 0.2, x => 3.0, y => 0.4, z => 2.0)
    H1_num = [2.0 0.3 0.1; 0.3 1.5 0.2; 0.1 0.2 2.0]
    H2_num = [3.0 0.4; 0.4 2.0]
    K_numeric = kron(H1_num, H2_num)
    true_eigs = eigvals(K_numeric)
    
    computed_eigs = [_eval_complex_sub(v, subs_dict) for v in vals]
    # For symmetric matrices, imaginary parts should be negligible
    @test all(abs.(imag.(computed_eigs)) .< 1e-10)
    @test sort(true_eigs) ≈ sort(real.(computed_eigs)) atol=1e-10
end

@testset "Full 6-variable 3×3 symmetric matrix (diagonal shift optimization)" begin
    # Test that the diagonal shift optimization allows solving a full 6-variable 3×3 symmetric matrix
    # This uses the shift trick: A - f*I has zero corner, reducing to 5 effective variables
    @variables a b c d e f
    A = [a b c; b d e; c e f]
    
    # This should succeed using the diagonal shift optimization
    vals, poly, λ = symbolic_eigenvalues(A; timeout=600, complexity_threshold=nothing)
    @test length(vals) == 3
    
    # Verify that eigenvalues contain all 6 original variables
    all_vars = Set{Symbol}()
    for v in vals
        for var in Symbolics.get_variables(real(v))
            push!(all_vars, Symbolics.tosymbol(var))
        end
    end
    @test all_vars == Set([:a, :b, :c, :d, :e, :f])
    
    # Numerical verification using the robust evaluation helper
    test_vals = Dict(a => 1.0, b => 0.5, c => 0.3, d => 2.0, e => 0.4, f => 3.0)
    A_num = Float64[test_vals[a] test_vals[b] test_vals[c];
                    test_vals[b] test_vals[d] test_vals[e];
                    test_vals[c] test_vals[e] test_vals[f]]
    true_eigs = eigvals(Symmetric(A_num))
    
    # Evaluate symbolic eigenvalues with complex-safe evaluation
    computed = [real(SymbolicDiagonalization._evaluate_symbolic_expr(v, test_vals)) for v in vals]
    
    # Verify eigenvalues match (for symmetric matrices, imaginary parts should be negligible)
    @test isapprox(sort(computed), sort(true_eigs), atol=1e-10)
end

@testset "Kronecker 6×6 = full 6-param 3×3 ⊗ 3-param 2×2 (9 parameters)" begin
    # Test Kronecker product with maximum parameters in both factors
    # 3×3 symmetric with 6 parameters (using diagonal shift optimization)
    @variables a b c d e f
    H1 = [a b c; b d e; c e f]
    
    # 2×2 symmetric with 3 parameters
    @variables x y z
    H2 = [x y; y z]
    
    # Kronecker product: 6×6 with 9 total parameters
    K = kron(H1, H2)
    
    # Should detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(K)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues - uses diagonal shift for the 3×3 factor
    vals, poly, λ = symbolic_eigenvalues(K; timeout=600, complexity_threshold=nothing)
    @test length(vals) == 6
    
    # Verify all 9 variables appear
    all_vars = Set{Symbol}()
    for v in vals
        expr = v isa Complex ? real(v) : v
        for var in Symbolics.get_variables(expr)
            push!(all_vars, Symbolics.tosymbol(var))
        end
    end
    @test all_vars == Set([:a, :b, :c, :d, :e, :f, :x, :y, :z])
    
    # Numerical verification
    test_vals = Dict(
        a => 1.0, b => 0.5, c => 0.3, d => 2.0, e => 0.4, f => 3.0,
        x => 1.5, y => 0.2, z => 2.5
    )
    H1_num = Float64[1.0 0.5 0.3; 0.5 2.0 0.4; 0.3 0.4 3.0]
    H2_num = Float64[1.5 0.2; 0.2 2.5]
    K_num = kron(H1_num, H2_num)
    true_eigs = eigvals(Symmetric(K_num))
    
    computed = [real(SymbolicDiagonalization._evaluate_symbolic_expr(v, test_vals)) for v in vals]
    @test isapprox(sort(computed), sort(true_eigs), atol=1e-8)
end

@testset "Nested Kronecker: 2×2 ⊗ 2×2 ⊗ 2×2 (8×8, 9 parameters)" begin
    # Test triple Kronecker product
    @variables a1 b1 c1  a2 b2 c2  a3 b3 c3
    
    A1 = [a1 b1; b1 c1]
    A2 = [a2 b2; b2 c2]
    A3 = [a3 b3; b3 c3]
    
    K = kron(kron(A1, A2), A3)
    @test size(K) == (8, 8)
    
    vals, _, _ = symbolic_eigenvalues(K; timeout=120, complexity_threshold=nothing)
    @test length(vals) == 8
    
    # Numerical verification
    test_vals = Dict(
        a1 => 1.0, b1 => 0.2, c1 => 2.0,
        a2 => 1.5, b2 => 0.3, c2 => 2.5,
        a3 => 1.2, b3 => 0.1, c3 => 1.8
    )
    A1_num = [1.0 0.2; 0.2 2.0]
    A2_num = [1.5 0.3; 0.3 2.5]
    A3_num = [1.2 0.1; 0.1 1.8]
    K_num = kron(kron(A1_num, A2_num), A3_num)
    true_eigs = eigvals(Symmetric(K_num))
    
    computed = [real(SymbolicDiagonalization._evaluate_symbolic_expr(v, test_vals)) for v in vals]
    @test isapprox(sort(computed), sort(true_eigs), atol=1e-8)
end

@testset "Deep nested Kronecker: 2×2^⊗5 (32×32, 15 parameters)" begin
    # Test 5-fold Kronecker product
    matrices = Matrix{Num}[]
    all_test_vals = Dict{Num, Float64}()
    
    for i in 1:5
        a = Symbolics.variable(Symbol("a$i"))
        b = Symbolics.variable(Symbol("b$i"))
        c = Symbolics.variable(Symbol("c$i"))
        push!(matrices, [a b; b c])
        all_test_vals[a] = 1.0 + 0.1*i
        all_test_vals[b] = 0.1 + 0.02*i
        all_test_vals[c] = 1.5 + 0.15*i
    end
    
    K = reduce(kron, matrices)
    @test size(K) == (32, 32)
    
    vals, _, _ = symbolic_eigenvalues(K; timeout=120, complexity_threshold=nothing)
    @test length(vals) == 32
    
    # Numerical verification
    numeric_matrices = [[1.0 + 0.1*i  0.1 + 0.02*i; 0.1 + 0.02*i  1.5 + 0.15*i] for i in 1:5]
    K_num = reduce(kron, numeric_matrices)
    true_eigs = eigvals(Symmetric(K_num))
    
    computed = [real(SymbolicDiagonalization._evaluate_symbolic_expr(v, all_test_vals)) for v in vals]
    @test isapprox(sort(computed), sort(true_eigs), atol=1e-7)
end

