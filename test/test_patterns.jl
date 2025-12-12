# ============================================================================
# Special Pattern Tests
# Tests for matrices with special patterns (circulant, Kronecker, Toeplitz, etc.)
# ============================================================================

# Helper function to convert symbolic Num to Float64
function num_to_float64(x)
    if x isa Symbolics.Num
        return Float64(eval(Symbolics.toexpr(x)))
    else
        return Float64(x)
    end
end

# ============================================================================
# Circulant Matrices
# ============================================================================

@testset "Circulant Matrices" begin
    # Test 3×3 circulant matrix with symbolic entries
    @variables a b c
    C3 = [a b c;
          c a b;
          b c a]
    
    @test SymbolicDiagonalization._is_circulant(C3)
    
    vals, poly, λ = symbolic_eigenvalues(C3)
    @test length(vals) == 3
    
    # For numeric verification, test with concrete values
    C3_num = [1 2 3; 3 1 2; 2 3 1]
    vals_num, _, _ = symbolic_eigenvalues(C3_num)
    @test length(vals_num) == 3
    
    # Verify eigenvalues numerically
    numeric_eigs = eigvals(float.(C3_num))
    computed_eigs = [complex(float(substitute(v, Dict()))) for v in vals_num]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "Circulant 4×4" begin
    @variables a b c d
    C4 = [a b c d;
          d a b c;
          c d a b;
          b c d a]
    
    @test SymbolicDiagonalization._is_circulant(C4)
    
    vals, poly, λ = symbolic_eigenvalues(C4)
    @test length(vals) == 4
    
    # Numeric test
    C4_num = [1 2 3 4; 4 1 2 3; 3 4 1 2; 2 3 4 1]
    vals_num, _, _ = symbolic_eigenvalues(C4_num)
    @test length(vals_num) == 4
    
    numeric_eigs = eigvals(float.(C4_num))
    computed_eigs = [complex(float(substitute(v, Dict()))) for v in vals_num]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "Circulant 5×5" begin
    # Test that circulant works for n > 4 (beyond quartic formula limit)
    C5_num = [1 2 3 4 5;
              5 1 2 3 4;
              4 5 1 2 3;
              3 4 5 1 2;
              2 3 4 5 1]
    
    @test SymbolicDiagonalization._is_circulant(C5_num)
    
    vals, poly, λ = symbolic_eigenvalues(C5_num)
    @test length(vals) == 5
    
    # Verify against numeric eigenvalues
    numeric_eigs = eigvals(float.(C5_num))
    computed_eigs = [complex(float(substitute(v, Dict()))) for v in vals]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

# ============================================================================
# Toeplitz Tridiagonal Matrices
# ============================================================================

@testset "Toeplitz Tridiagonal Matrices" begin
    # Test 3×3 symmetric tridiagonal (Toeplitz)
    @variables a b
    T3 = [a b 0;
          b a b;
          0 b a]
    
    params = SymbolicDiagonalization._is_toeplitz_tridiagonal(T3)
    @test !isnothing(params)
    @test isequal(params[1], a)  # diagonal
    @test isequal(params[2], b)  # subdiagonal
    @test isequal(params[3], b)  # superdiagonal
    
    vals, poly, λ = symbolic_eigenvalues(T3)
    @test length(vals) == 3
    
    # Numeric verification
    T3_num = [2 1 0; 1 2 1; 0 1 2]
    vals_num, _, _ = symbolic_eigenvalues(T3_num)
    @test length(vals_num) == 3
    
    numeric_eigs = eigvals(float.(T3_num))
    computed_eigs = [float(substitute(v, Dict())) for v in vals_num]
    @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Toeplitz Tridiagonal 5×5" begin
    # Test larger tridiagonal (beyond quartic limit)
    T5_num = [2 1 0 0 0;
              1 2 1 0 0;
              0 1 2 1 0;
              0 0 1 2 1;
              0 0 0 1 2]
    
    params = SymbolicDiagonalization._is_toeplitz_tridiagonal(T5_num)
    @test !isnothing(params)
    @test params == (2, 1, 1)
    
    vals, poly, λ = symbolic_eigenvalues(T5_num)
    @test length(vals) == 5
    
    # Verify against numeric eigenvalues
    numeric_eigs = eigvals(float.(T5_num))
    computed_eigs = [float(substitute(v, Dict())) for v in vals]
    @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
end

# ============================================================================
# Anti-Diagonal Matrices
# ============================================================================

@testset "Anti-Diagonal Matrices" begin
    # Test 3×3 symmetric anti-diagonal matrix
    @variables a b c
    A3 = [0 0 a;
          0 b 0;
          a 0 0]
    
    @test SymbolicDiagonalization._is_antidiagonal(A3)
    @test SymbolicDiagonalization._is_symmetric(A3)
    
    vals, poly, λ = symbolic_eigenvalues(A3)
    @test length(vals) == 3
    
    # Numeric verification
    A3_num = [0 0 3.0; 0 2.0 0; 3.0 0 0]
    vals_num, _, _ = symbolic_eigenvalues(A3_num)
    @test length(vals_num) == 3
    
    numeric_eigs = eigvals(float.(A3_num))
    computed_eigs = [float(substitute(v, Dict())) for v in vals_num]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "Anti-Diagonal 4×4" begin
    # Test 4×4 symmetric anti-diagonal matrix
    @variables a b
    A4 = [0 0 0 a;
          0 0 b 0;
          0 b 0 0;
          a 0 0 0]
    
    @test SymbolicDiagonalization._is_antidiagonal(A4)
    @test SymbolicDiagonalization._is_symmetric(A4)
    
    vals, poly, λ = symbolic_eigenvalues(A4)
    @test length(vals) == 4
    
    # Numeric verification
    A4_num = [0 0 0 3.0; 0 0 2.0 0; 0 2.0 0 0; 3.0 0 0 0]
    vals_num, _, _ = symbolic_eigenvalues(A4_num)
    @test length(vals_num) == 4
    
    numeric_eigs = eigvals(float.(A4_num))
    computed_eigs = [float(substitute(v, Dict())) for v in vals_num]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "Anti-Diagonal 5×5" begin
    # Test 5×5 symmetric anti-diagonal matrix
    @variables a b c
    A5 = [0 0 0 0 a;
          0 0 0 b 0;
          0 0 c 0 0;
          0 b 0 0 0;
          a 0 0 0 0]
    
    @test SymbolicDiagonalization._is_antidiagonal(A5)
    @test SymbolicDiagonalization._is_symmetric(A5)
    
    vals, poly, λ = symbolic_eigenvalues(A5)
    @test length(vals) == 5
    
    # Numeric verification
    A5_num = [0 0 0 0 5.0;
              0 0 0 3.0 0;
              0 0 2.0 0 0;
              0 3.0 0 0 0;
              5.0 0 0 0 0]
    vals_num, _, _ = symbolic_eigenvalues(A5_num)
    @test length(vals_num) == 5
    
    numeric_eigs = eigvals(float.(A5_num))
    computed_eigs = [float(substitute(v, Dict())) for v in vals_num]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

# ============================================================================
# Permutation Matrix Pattern
# ============================================================================

@testset "Permutation Matrix - Identity" begin
    # Identity matrix: all fixed points
    I3 = [1 0 0; 0 1 0; 0 0 1]
    vals, poly, _ = symbolic_eigenvalues(I3)
    @test length(vals) == 3
    @test all(isapprox.(vals, 1.0))
end

@testset "Permutation Matrix - Single 2-Cycle" begin
    # Permutation: 1↔2, 3 fixed
    P = [0 1 0; 1 0 0; 0 0 1]
    vals, poly, _ = symbolic_eigenvalues(P)
    @test length(vals) == 3
    vals_sorted = sort(real.(vals))
    @test isapprox(vals_sorted[1], -1.0, atol=1e-10)
    @test isapprox(vals_sorted[2], 1.0, atol=1e-10)
    @test isapprox(vals_sorted[3], 1.0, atol=1e-10)
end

@testset "Permutation Matrix - Single 3-Cycle" begin
    # Permutation: 1→2→3→1
    P = [0 0 1; 1 0 0; 0 1 0]
    vals, poly, _ = symbolic_eigenvalues(P)
    @test length(vals) == 3
    
    # Check that eigenvalues are 3rd roots of unity
    vals_cubed = [v^3 for v in vals]
    @test all(isapprox.(vals_cubed, 1.0, atol=1e-10))
    
    # Check sum is 0 (property of 3rd roots of unity)
    @test isapprox(sum(vals), 0.0, atol=1e-10)
end

@testset "Permutation Matrix - Single 4-Cycle" begin
    # Permutation: 1→2→3→4→1
    P = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0]
    vals, poly, _ = symbolic_eigenvalues(P)
    @test length(vals) == 4
    
    # Check that eigenvalues are 4th roots of unity
    vals_fourth = [v^4 for v in vals]
    @test all(isapprox.(vals_fourth, 1.0, atol=1e-10))
    
    # Check sum is 0
    @test isapprox(sum(vals), 0.0, atol=1e-10)
end

# ============================================================================
# Block Circulant Matrix Tests (Numeric Only)
# ============================================================================

@testset "Block Circulant - 4×4 with 2×2 blocks" begin
    # [A B]
    # [B A]
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    M = [A B; B A]
    
    # Verify detection
    info = SymbolicDiagonalization._is_block_circulant(M)
    @test !isnothing(info)
    n_blocks, block_size, blocks = info
    @test n_blocks == 2
    @test block_size == 2
    @test blocks[1] == A
    @test blocks[2] == B
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Verify against Julia's eigvals
    true_vals = eigvals(complex(float.(M)))
    computed_vals = [complex(float(v)) for v in vals]
    @test sort(real(true_vals)) ≈ sort(real(computed_vals)) atol=1e-10
    @test sort(imag(true_vals)) ≈ sort(imag(computed_vals)) atol=1e-10
end

@testset "Block Circulant - 6×6 with 2×2 blocks" begin
    # 3 blocks of size 2×2: A, B, C
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    C = [9 10; 11 12]
    M = [A B C; C A B; B C A]
    
    # Verify detection
    info = SymbolicDiagonalization._is_block_circulant(M)
    @test !isnothing(info)
    n_blocks, block_size, blocks = info
    @test n_blocks == 3
    @test block_size == 2
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 6
    
    # Verify against Julia's eigvals
    true_vals = eigvals(complex(float.(M)))
    computed_vals = [complex(float(v)) for v in vals]
    @test sort(real(true_vals)) ≈ sort(real(computed_vals)) atol=1e-10
    @test sort(imag(true_vals)) ≈ sort(imag(computed_vals)) atol=1e-10
end

@testset "Block Circulant - Diagonal Blocks" begin
    # Blocks are diagonal matrices
    A = [2.0 0.0; 0.0 3.0]
    B = [4.0 0.0; 0.0 5.0]
    M = [A B; B A]
    
    info = SymbolicDiagonalization._is_block_circulant(M)
    @test !isnothing(info)
    
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Theory: D₀ = A+B = diag(6, 8), D₁ = A-B = diag(-2, -2)
    # Eigenvalues should be: 6, 8, -2, -2
    expected = sort([6.0, 8.0, -2.0, -2.0])
    computed = sort(real.(vals))
    @test computed ≈ expected atol=1e-10
end

# ============================================================================
# Kronecker Product Tests
# ============================================================================

@testset "Kronecker Product - 4×4 = 2×2 ⊗ 2×2 (numeric)" begin
    # Simple numeric test: A ⊗ B where both are 2×2
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    M = kron(A, B)  # 4×4 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    A_detected, B_detected, m, n = kron_info
    @test m == 2
    @test n == 2
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Verify against Julia's eigvals
    true_vals = eigvals(M)
    computed_vals = [num_to_float64(real(v)) for v in vals]
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - 6×6 = 2×2 ⊗ 3×3 (numeric)" begin
    # Non-square factors
    A = [1 2; 3 4]
    B = [1 0 0; 0 2 0; 0 0 3]  # Diagonal for easier verification
    M = kron(A, B)  # 6×6 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    A_detected, B_detected, m, n = kron_info
    @test m == 2
    @test n == 3
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 6
    
    # Verify against Julia's eigvals
    true_vals = eigvals(M)
    computed_vals = [num_to_float64(real(v)) for v in vals]
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - 6×6 = 3×3 ⊗ 2×2 (numeric)" begin
    # Reverse order from previous test
    A = [1 0 0; 0 2 0; 0 0 3]  # Diagonal
    B = [1 2; 3 4]
    M = kron(A, B)  # 6×6 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 6
    
    # Verify against Julia's eigvals
    true_vals = eigvals(M)
    computed_vals = [num_to_float64(real(v)) for v in vals]
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - 9×9 = 3×3 ⊗ 3×3 (numeric)" begin
    # Larger example
    A = [1 2 0; 0 3 1; 0 0 4]  # Upper triangular
    B = [5 0 0; 6 7 0; 8 9 10]  # Lower triangular
    M = kron(A, B)  # 9×9 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    A_detected, B_detected, m, n = kron_info
    @test m == 3
    @test n == 3
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 9
    
    true_vals = eigvals(M)
    computed_vals = [num_to_float64(real(v)) for v in vals]
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - 8×8 = 4×4 ⊗ 2×2 (numeric)" begin
    # Use diagonal 4×4 factor for easier computation
    A = [1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 4]
    B = [1 2; 3 4]
    M = kron(A, B)  # 8×8 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 8
    
    # Verify against Julia's eigvals
    true_vals = eigvals(M)
    computed_vals = [num_to_float64(real(v)) for v in vals]
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - 8×8 = 2×2 ⊗ 4×4 (numeric)" begin
    # Reverse order - small ⊗ large
    A = [1 2; 3 4]
    B = [1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 4]
    M = kron(A, B)  # 8×8 matrix
    
    # Detect Kronecker structure  
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 8
    
    # Verify against Julia's eigvals
    true_vals = eigvals(M)
    computed_vals = [num_to_float64(real(v)) for v in vals]
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - 4×4 = 2×2 ⊗ 2×2 (symbolic)" begin
    # Symbolic test with variables
    @variables a b c d
    A = [a 0; 0 b]  # Diagonal for simpler expressions
    B = [c 0; 0 d]  # Diagonal for simpler expressions
    M = kron(A, B)  # 4×4 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # For diagonal matrices, eigenvalues should be products: {ac, ad, bc, bd}
    expected = [a*c, a*d, b*c, b*d]
    
    # Check that computed eigenvalues match expected
    vals_simplified = [Symbolics.simplify(v) for v in vals]
    @test Set(vals_simplified) == Set(expected)
end

@testset "Kronecker Product - Identity factors" begin
    # I₂ ⊗ A should have eigenvalues [λ, λ] for each eigenvalue λ of A
    I2 = [1 0; 0 1]
    A = [1 2; 3 4]
    M = kron(I2, A)  # 4×4 matrix
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Each eigenvalue of A should appear twice (multiplicity 2)
    true_vals = eigvals(M)
    computed_vals = Float64.(real.(vals))
    @test sort(true_vals) ≈ sort(computed_vals) atol=1e-10
end

@testset "Kronecker Product - Complex numeric 4×4 = 2×2 ⊗ 2×2" begin
    # Test with complex numeric matrices
    A = ComplexF64[1.0+0.0im 2.0-1.0im; 2.0+1.0im -0.5+0.0im]
    B = ComplexF64[-1.5+0.0im 1.0+1.0im; 1.0-1.0im 1.0+0.0im]
    M = kron(A, B)
    
    # Detect Kronecker structure
    kron_info = SymbolicDiagonalization._is_kronecker_product(M)
    @test !isnothing(kron_info)
    
    # Compute eigenvalues
    vals, poly, _ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Verify against Julia's eigvals
    true_vals = eigvals(M)
    # Sort by absolute value for comparison
    computed_sorted = sort(vals, by=abs)
    true_sorted = sort(true_vals, by=abs)
    @test computed_sorted ≈ true_sorted atol=1e-10
end

# ============================================================================
# Edge Cases
# ============================================================================

@testset "Edge Cases" begin
    # 1×1 matrix (trivially circulant)
    @variables a
    M1 = [a]
    @test SymbolicDiagonalization._is_circulant(M1)
    
    # 2×2 symmetric circulant
    @variables b
    C2 = [a b; b a]
    @test SymbolicDiagonalization._is_circulant(C2)
    
    # 2×2 symmetric tridiagonal
    T2 = [a b; b a]
    params = SymbolicDiagonalization._is_toeplitz_tridiagonal(T2)
    @test !isnothing(params)
    @test isequal(params[1], a)
    @test isequal(params[2], b)
    @test isequal(params[3], b)
end

@testset "Non-Special Matrices" begin
    # Test that non-circulant matrices are correctly rejected
    @variables a b c d
    M1 = [a b c; d a b; b c a]  # Not circulant
    @test !SymbolicDiagonalization._is_circulant(M1)
    
    # Test that non-Toeplitz tridiagonal is rejected
    M2 = [a b 0; b c b; 0 b a]  # Diagonal not constant
    @test isnothing(SymbolicDiagonalization._is_toeplitz_tridiagonal(M2))
    
    M3 = [a b 0; b a c; 0 b a]  # Not symmetric
    @test isnothing(SymbolicDiagonalization._is_toeplitz_tridiagonal(M3))
    
    # Jordan block - tridiagonal but not symmetric
    M4 = [a 1; 0 a]
    @test isnothing(SymbolicDiagonalization._is_toeplitz_tridiagonal(M4))
    
    # Non-anti-diagonal matrices
    M5 = [a 0 0; 0 b 0; 0 0 c]  # Diagonal, not anti-diagonal
    @test !SymbolicDiagonalization._is_antidiagonal(M5)
end
# ============================================================================
# Group Theory Pattern Tests
# Tests for hypercube graphs, strongly regular graphs, etc.
# ============================================================================

@testset "Hypercube Graphs" begin
    # Test Q_1: 2×2 hypercube (single edge)
    Q1 = [0 1; 1 0]
    @test SymbolicDiagonalization._is_hypercube_graph(Q1) == 1
    
    vals, poly, λ = symbolic_eigenvalues(Q1)
    @test length(vals) == 2
    # Q_1 eigenvalues: 1 - 2*0 = 1, 1 - 2*1 = -1
    @test sort([1, -1]) == sort(vals)
    
    # Test Q_2: 4×4 hypercube (square)
    Q2 = [0 1 1 0;
          1 0 0 1;
          1 0 0 1;
          0 1 1 0]
    @test SymbolicDiagonalization._is_hypercube_graph(Q2) == 2
    
    vals, poly, λ = symbolic_eigenvalues(Q2)
    @test length(vals) == 4
    # Q_2 eigenvalues: 2, 0, 0, -2 (with multiplicities from binomial coefficients)
    expected = [2, 0, 0, -2]  # λ_k = 2 - 2k for k=0,1,1,2
    @test sort(vals) == sort(expected)
    
    # Test Q_3: 8×8 hypercube (cube)
    # Structure: [Q_2  I_4; I_4  Q_2]
    I4 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    Q3 = [Q2 I4; I4 Q2]
    @test SymbolicDiagonalization._is_hypercube_graph(Q3) == 3
    
    vals, poly, λ = symbolic_eigenvalues(Q3)
    @test length(vals) == 8
    # Q_3 eigenvalues: 3, 1, 1, 1, -1, -1, -1, -3
    # λ_k = 3 - 2k for k=0,1,2,3 with mult binomial(3,k) = 1,3,3,1
    expected = [3, 1, 1, 1, -1, -1, -1, -3]
    @test sort(vals) == sort(expected)
    
    # Verify eigenvalues numerically
    numeric_eigs = eigvals(float.(Q3))
    @test sort(real(numeric_eigs)) ≈ sort(real(vals)) atol=1e-10
end

@testset "Hypercube Q_4 (16×16)" begin
    # Build Q_4 recursively
    I8 = Matrix{Int}(I, 8, 8)
    Q2 = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0]
    I4 = Matrix{Int}(I, 4, 4)
    Q3 = [Q2 I4; I4 Q2]
    Q4 = [Q3 I8; I8 Q3]
    
    @test SymbolicDiagonalization._is_hypercube_graph(Q4) == 4
    
    vals, poly, λ = symbolic_eigenvalues(Q4)
    @test length(vals) == 16
    
    # Q_4 eigenvalues: λ_k = 4 - 2k for k=0..4 with mult binomial(4,k)
    # Multiplicities: 1, 4, 6, 4, 1
    expected = vcat([4], fill(2, 4), fill(0, 6), fill(-2, 4), [-4])
    @test sort(vals) == sort(expected)
    
    # Verify numerically
    numeric_eigs = eigvals(float.(Q4))
    @test sort(real(numeric_eigs)) ≈ sort(real(vals)) atol=1e-10
end

@testset "Strongly Regular Graphs" begin
    # Test Petersen graph: srg(10, 3, 0, 1)
    # The Petersen graph is a well-known strongly regular graph
    # Adjacency matrix (vertex numbering may vary)
    Petersen = [
        0 1 0 0 1 1 0 0 0 0;
        1 0 1 0 0 0 1 0 0 0;
        0 1 0 1 0 0 0 1 0 0;
        0 0 1 0 1 0 0 0 1 0;
        1 0 0 1 0 0 0 0 0 1;
        1 0 0 0 0 0 0 1 1 0;
        0 1 0 0 0 0 0 0 1 1;
        0 0 1 0 0 1 0 0 0 1;
        0 0 0 1 0 1 1 0 0 0;
        0 0 0 0 1 0 1 1 0 0
    ]
    
    # Check if detected as strongly regular
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(Petersen)
    @test !isnothing(srg_params)
    n, k, λ_param, μ = srg_params
    @test n == 10
    @test k == 3
    @test λ_param == 0
    @test μ == 1
    
    # Compute eigenvalues
    vals, poly, λ = symbolic_eigenvalues(Petersen)
    @test length(vals) == 10
    
    # Petersen graph has eigenvalues: 3 (mult 1), 1 (mult 5), -2 (mult 4)
    # Count each eigenvalue
    vals_int = round.(Int, real.(vals))
    @test count(==(3), vals_int) == 1
    @test count(==(1), vals_int) == 5
    @test count(==(-2), vals_int) == 4
    
    # Verify numerically
    numeric_eigs = eigvals(float.(Petersen))
    @test sort(real(numeric_eigs)) ≈ sort(real(vals)) atol=1e-10
end

@testset "Complete Graph K_5 as SRG" begin
    # K_5 is srg(5, 4, 3, 4) - complete graph on 5 vertices
    K5 = ones(Int, 5, 5) - I(5)
    
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(K5)
    @test !isnothing(srg_params)
    n, k, λ_param, μ = srg_params
    @test n == 5
    @test k == 4
    @test λ_param == 3
    @test μ == 4
    
    vals, poly, λ = symbolic_eigenvalues(K5)
    @test length(vals) == 5
    
    # K_5 eigenvalues: 4 (mult 1), -1 (mult 4)
    vals_int = round.(Int, real.(vals))
    @test count(==(4), vals_int) == 1
    @test count(==(-1), vals_int) == 4
end


@testset "Hypercube Q_5 (32×32)" begin
    # Build Q_5 recursively from Q_4
    I16 = Matrix{Int}(I, 16, 16)
    Q2 = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0]
    I4 = Matrix{Int}(I, 4, 4)
    Q3 = [Q2 I4; I4 Q2]
    I8 = Matrix{Int}(I, 8, 8)
    Q4 = [Q3 I8; I8 Q3]
    Q5 = [Q4 I16; I16 Q4]
    
    @test SymbolicDiagonalization._is_hypercube_graph(Q5) == 5
    
    vals, poly, λ = symbolic_eigenvalues(Q5)
    @test length(vals) == 32
    
    # Q_5 eigenvalues: λ_k = 5 - 2k for k=0..5 with mult binomial(5,k)
    # Multiplicities: 1, 5, 10, 10, 5, 1
    expected = vcat([5], fill(3, 5), fill(1, 10), fill(-1, 10), fill(-3, 5), [-5])
    @test sort(vals) == sort(expected)
    
    # Verify numerically
    numeric_eigs = eigvals(float.(Q5))
    @test sort(real(numeric_eigs)) ≈ sort(real(vals)) atol=1e-10
end

@testset "Hypercube Edge Cases" begin
    # Q_0: 1×1 zero matrix
    Q0 = zeros(Int, 1, 1)
    @test SymbolicDiagonalization._is_hypercube_graph(Q0) == 0
    vals, poly, λ = symbolic_eigenvalues(Q0)
    @test vals == [0]
    
    # Non-hypercube: 3×3 matrix (not a power of 2)
    M3 = [0 1 1; 1 0 1; 1 1 0]
    @test isnothing(SymbolicDiagonalization._is_hypercube_graph(M3))
    
    # Non-hypercube: 4×4 matrix that's not Q_2
    M4 = [0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0]  # K_4 complete graph
    @test isnothing(SymbolicDiagonalization._is_hypercube_graph(M4))
end

@testset "Cycle Graphs as SRG" begin
    # C_5 (pentagon) is srg(5, 2, 0, 1)
    C5 = [
        0 1 0 0 1;
        1 0 1 0 0;
        0 1 0 1 0;
        0 0 1 0 1;
        1 0 0 1 0
    ]
    
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(C5)
    @test !isnothing(srg_params)
    n, k, λ_param, μ = srg_params
    @test n == 5
    @test k == 2
    @test λ_param == 0
    @test μ == 1
    
    vals, poly, λ = symbolic_eigenvalues(C5)
    @test length(vals) == 5
    
    # C_5 has eigenvalues involving golden ratio (can't use round to Int)
    # Just verify they match numerically
    numeric_eigs = eigvals(float.(C5))
    @test sort(real(numeric_eigs)) ≈ sort(real(vals)) atol=1e-10
end

@testset "Complete Bipartite Graphs K_{m,n}" begin
    # K_{3,3} is srg(6, 3, 0, 3)
    K33 = [
        0 0 0 1 1 1;
        0 0 0 1 1 1;
        0 0 0 1 1 1;
        1 1 1 0 0 0;
        1 1 1 0 0 0;
        1 1 1 0 0 0
    ]
    
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(K33)
    @test !isnothing(srg_params)
    n, k, λ_param, μ = srg_params
    @test n == 6
    @test k == 3
    @test λ_param == 0
    @test μ == 3
    
    vals, poly, λ = symbolic_eigenvalues(K33)
    @test length(vals) == 6
    
    # K_{m,n} has eigenvalues: ±√(mn) and 0 with multiplicity m+n-2
    vals_sorted = sort(real.(vals))
    # K_{3,3}: eigenvalues are 3, -3, 0, 0, 0, 0
    @test count(x -> abs(x - 3) < 1e-10, vals_sorted) == 1
    @test count(x -> abs(x + 3) < 1e-10, vals_sorted) == 1
    @test count(x -> abs(x) < 1e-10, vals_sorted) == 4
end

@testset "Additional SRG Examples" begin
    # Test a few simpler SRGs without specific parameter assertions
    
    # Complete graph K_4 as SRG
    K4 = ones(Int, 4, 4) - I(4)
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(K4)
    @test !isnothing(srg_params)
    vals, poly, λ = symbolic_eigenvalues(K4)
    @test length(vals) == 4
    numeric_eigs = eigvals(float.(K4))
    @test sort(real(numeric_eigs)) ≈ sort(real(vals)) atol=1e-10
    
    # 4-cycle (square) is srg(4, 2, 0, 2)
    C4 = [
        0 1 0 1;
        1 0 1 0;
        0 1 0 1;
        1 0 1 0
    ]
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(C4)
    @test !isnothing(srg_params)
    n, k, λ_param, μ = srg_params
    @test n == 4
    @test k == 2
    @test λ_param == 0
    @test μ == 2
    vals, poly, λ = symbolic_eigenvalues(C4)
    @test length(vals) == 4
end

@testset "SRG Non-Examples" begin
    # Path graph P_4 (not strongly regular - different vertex pairs have different common neighbors)
    P4 = [
        0 1 0 0;
        1 0 1 0;
        0 1 0 1;
        0 0 1 0
    ]
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(P4)
    @test isnothing(srg_params)  # P_4 is not strongly regular
    
    # Star graph K_{1,3} (not even regular)
    Star = [
        0 1 1 1;
        1 0 0 0;
        1 0 0 0;
        1 0 0 0
    ]
    srg_params = SymbolicDiagonalization._is_strongly_regular_graph(Star)
    @test isnothing(srg_params)  # Not regular
end

@testset "Performance: SRG vs Generic Method" begin
    # For small SRG, verify that pattern detection is used
    K5 = ones(Int, 5, 5) - I(5)
    
    # Time pattern-based method
    vals1, poly1, λ1 = symbolic_eigenvalues(K5; structure=:auto)
    
    # Verify the pattern method gives correct results
    @test length(vals1) == 5
    vals_int = round.(Int, real.(vals1))
    @test count(==(4), vals_int) == 1
    @test count(==(-1), vals_int) == 4
end
