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
# Circulant Eigenvectors (DFT Basis)
# ============================================================================

@testset "Circulant Eigenvectors" begin
    @testset "Circulant 3×3 Eigenvectors" begin
        C3 = [1 2 3; 3 1 2; 2 3 1]
        pairs, poly, λ = symbolic_eigenpairs(C3)
        
        @test length(pairs) == 3
        
        # Verify each eigenpair: C * v = λ * v
        for (val, vecs) in pairs
            v = complex.(vecs[1])
            Cv = C3 * v
            λv = complex(val) * v
            @test norm(Cv - λv) < 1e-10
        end
    end
    
    @testset "Circulant 4×4 Eigenvectors" begin
        C4 = [1 2 3 4; 4 1 2 3; 3 4 1 2; 2 3 4 1]
        pairs, poly, λ = symbolic_eigenpairs(C4)
        
        @test length(pairs) == 4
        
        # Verify eigenpairs
        for (val, vecs) in pairs
            v = complex.(vecs[1])
            Cv = C4 * v
            λv = complex(val) * v
            @test norm(Cv - λv) < 1e-10
        end
        
        # Eigenvectors should be DFT columns (orthogonal after normalization)
        vecs = [complex.(pairs[k][2][1]) for k in 1:4]
        for i in 1:4
            for j in i+1:4
                # DFT columns are orthogonal
                inner = dot(vecs[i], vecs[j])
                @test abs(inner) < 1e-10
            end
        end
    end
    
    @testset "Circulant 6×6 Eigenvectors" begin
        C6 = [1 2 3 4 5 6; 6 1 2 3 4 5; 5 6 1 2 3 4; 
              4 5 6 1 2 3; 3 4 5 6 1 2; 2 3 4 5 6 1]
        pairs, poly, λ = symbolic_eigenpairs(C6)
        
        @test length(pairs) == 6
        
        # Verify eigenpairs
        for (val, vecs) in pairs
            v = complex.(vecs[1])
            Cv = C6 * v
            λv = complex(val) * v
            @test norm(Cv - λv) < 1e-10
        end
    end
    
    @testset "Circulant Eigenvector Orthogonality" begin
        # For any size n, DFT columns should be orthogonal
        for n in [3, 5, 7, 8]
            # Build n×n circulant with first row [1, 2, ..., n]
            first_row = collect(1:n)
            C = zeros(Int, n, n)
            for i in 1:n
                for j in 1:n
                    C[i, j] = first_row[mod1(j - i + 1, n)]
                end
            end
            
            pairs, _, _ = symbolic_eigenpairs(C)
            @test length(pairs) == n
            
            # Verify orthogonality of eigenvectors
            vecs = [complex.(pairs[k][2][1]) for k in 1:n]
            for i in 1:n
                for j in i+1:n
                    inner = dot(vecs[i], vecs[j])
                    @test abs(inner) < 1e-10
                end
            end
        end
    end
    
    @testset "DFT Basis Helper Functions" begin
        # Test _dft_column directly
        for n in [4, 5, 6]
            for k in 0:n-1
                col = SymbolicDiagonalization._dft_column(n, k)
                @test length(col) == n
                @test col[1] == 1  # First entry is always 1
            end
        end
        
        # Test _circulant_eigenvectors
        for n in [3, 4, 5]
            vecs = SymbolicDiagonalization._circulant_eigenvectors(n)
            @test length(vecs) == n
            @test all(length(v) == n for v in vecs)
        end
    end
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


# ============================================================================
# Dihedral Group Dₙ - Symmetric Circulant Matrices
# ============================================================================

@testset "Symmetric Circulant - Dihedral D₃" begin
    # 3×3 symmetric circulant: first row [a, b, b] (palindromic)
    @variables a b
    D3 = [a b b;
          b a b;
          b b a]
    
    # Check detection
    first_row = SymbolicDiagonalization._is_symmetric_circulant(D3)
    @test !isnothing(first_row)
    @test length(first_row) == 3
    
    vals, poly, λ = symbolic_eigenvalues(D3)
    @test length(vals) == 3
    
    # Numeric verification
    D3_num = [1.0 2.0 2.0; 2.0 1.0 2.0; 2.0 2.0 1.0]
    vals_num, _, _ = symbolic_eigenvalues(D3_num)
    @test length(vals_num) == 3
    
    numeric_eigs = eigvals(D3_num)
    computed_eigs = [real(float(substitute(v, Dict()))) for v in vals_num]
    @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Symmetric Circulant - Dihedral D₄" begin
    # 4×4 symmetric circulant: first row [a, b, c, b]
    @variables a b c
    D4 = [a b c b;
          b a b c;
          c b a b;
          b c b a]
    
    first_row = SymbolicDiagonalization._is_symmetric_circulant(D4)
    @test !isnothing(first_row)
    
    vals, poly, λ = symbolic_eigenvalues(D4)
    @test length(vals) == 4
    
    # Numeric verification
    D4_num = [1.0 2.0 3.0 2.0; 2.0 1.0 2.0 3.0; 3.0 2.0 1.0 2.0; 2.0 3.0 2.0 1.0]
    vals_num, _, _ = symbolic_eigenvalues(D4_num)
    
    numeric_eigs = eigvals(D4_num)
    computed_eigs = [real(float(substitute(v, Dict()))) for v in vals_num]
    @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Symmetric Circulant - Dihedral D₅" begin
    # 5×5 symmetric circulant (larger than quartic limit)
    # First row [a, b, c, c, b]
    @variables a b c
    D5 = [a b c c b;
          b a b c c;
          c b a b c;
          c c b a b;
          b c c b a]
    
    first_row = SymbolicDiagonalization._is_symmetric_circulant(D5)
    @test !isnothing(first_row)
    
    vals, poly, λ = symbolic_eigenvalues(D5)
    @test length(vals) == 5
    
    # Numeric verification
    D5_num = [1.0 2.0 3.0 3.0 2.0; 2.0 1.0 2.0 3.0 3.0; 3.0 2.0 1.0 2.0 3.0; 
              3.0 3.0 2.0 1.0 2.0; 2.0 3.0 3.0 2.0 1.0]
    vals_num, _, _ = symbolic_eigenvalues(D5_num)
    
    numeric_eigs = eigvals(D5_num)
    computed_eigs = [real(float(substitute(v, Dict()))) for v in vals_num]
    @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Symmetric Circulant - Non-Example" begin
    # Not symmetric circulant: first row is not palindromic
    @variables a b c
    M = [a b c; c a b; b c a]  # Regular circulant, not symmetric
    
    first_row = SymbolicDiagonalization._is_symmetric_circulant(M)
    @test isnothing(first_row)  # Should fail palindrome check
    
    # But should still be detected as circulant
    @test SymbolicDiagonalization._is_circulant(M)
end

# ============================================================================
# BCCB - Block Circulant with Circulant Blocks (Zₙ × Zₘ)
# ============================================================================

@testset "BCCB - 4×4 with 2×2 Circulant Blocks" begin
    # BCCB with 2 blocks of size 2: both block circulant and each block circulant
    # Block 0: [a b; b a], Block 1: [c d; d c]
    @variables a b c d
    B0 = [a b; b a]
    B1 = [c d; d c]
    M = [B0 B1; B1 B0]  # Block circulant with circulant blocks
    
    # Check detection
    bccb_info = SymbolicDiagonalization._is_bccb(M)
    @test !isnothing(bccb_info)
    n, m, first_rows = bccb_info
    @test n == 2
    @test m == 2
    
    vals, poly, λ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Numeric verification
    B0_num = [1.0 2.0; 2.0 1.0]
    B1_num = [3.0 4.0; 4.0 3.0]
    M_num = [B0_num B1_num; B1_num B0_num]
    vals_num, _, _ = symbolic_eigenvalues(M_num)
    
    numeric_eigs = eigvals(M_num)
    computed_eigs = [complex(float(substitute(v, Dict()))) for v in vals_num]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "BCCB - 9×9 with 3×3 Circulant Blocks" begin
    # 3 blocks of size 3, each block is circulant
    B0_num = [1.0 2.0 3.0; 3.0 1.0 2.0; 2.0 3.0 1.0]
    B1_num = [4.0 5.0 6.0; 6.0 4.0 5.0; 5.0 6.0 4.0]
    B2_num = [7.0 8.0 9.0; 9.0 7.0 8.0; 8.0 9.0 7.0]
    
    M = [B0_num B1_num B2_num;
         B2_num B0_num B1_num;
         B1_num B2_num B0_num]
    
    # Check detection
    bccb_info = SymbolicDiagonalization._is_bccb(M)
    @test !isnothing(bccb_info)
    n, m, first_rows = bccb_info
    @test n == 3
    @test m == 3
    
    vals, poly, λ = symbolic_eigenvalues(M)
    @test length(vals) == 9
    
    # Verify against Julia's eigvals
    numeric_eigs = eigvals(M)
    computed_eigs = [complex(float(substitute(v, Dict()))) for v in vals]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "BCCB - Non-Example (Block Circulant but Blocks Not Circulant)" begin
    # Block circulant, but individual blocks are not circulant
    A = [1 2; 3 4]  # Not circulant
    B = [5 6; 7 8]  # Not circulant
    M = [A B; B A]
    
    # Should be detected as block circulant
    block_info = SymbolicDiagonalization._is_block_circulant(M)
    @test !isnothing(block_info)
    
    # But NOT as BCCB
    bccb_info = SymbolicDiagonalization._is_bccb(M)
    @test isnothing(bccb_info)
end

# ============================================================================
# Polygon Adjacency - Cycle Graphs Cₙ
# ============================================================================

@testset "Polygon Adjacency - Triangle C₃" begin
    C3 = [0 1 1; 1 0 1; 1 1 0]
    
    polygon_n = SymbolicDiagonalization._is_polygon_adjacency(C3)
    @test polygon_n == 3
    
    vals, poly, λ = symbolic_eigenvalues(C3)
    @test length(vals) == 3
    
    # Triangle eigenvalues: 2, -1, -1
    vals_sorted = sort(real.(vals))
    @test vals_sorted ≈ [-1.0, -1.0, 2.0] atol=1e-10
end

@testset "Polygon Adjacency - Square C₄" begin
    C4 = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]
    
    polygon_n = SymbolicDiagonalization._is_polygon_adjacency(C4)
    @test polygon_n == 4
    
    vals, poly, λ = symbolic_eigenvalues(C4)
    @test length(vals) == 4
    
    # Square eigenvalues: 2, 0, -2, 0
    vals_sorted = sort(real.(vals))
    @test vals_sorted ≈ [-2.0, 0.0, 0.0, 2.0] atol=1e-10
end

@testset "Polygon Adjacency - Pentagon C₅" begin
    C5 = [0 1 0 0 1;
          1 0 1 0 0;
          0 1 0 1 0;
          0 0 1 0 1;
          1 0 0 1 0]
    
    polygon_n = SymbolicDiagonalization._is_polygon_adjacency(C5)
    @test polygon_n == 5
    
    vals, poly, λ = symbolic_eigenvalues(C5)
    @test length(vals) == 5
    
    # Pentagon eigenvalues: 2, (√5-1)/2, (√5-1)/2, -(√5+1)/2, -(√5+1)/2
    # ≈ 2, 0.618, 0.618, -1.618, -1.618
    numeric_eigs = eigvals(float.(C5))
    computed_eigs = [real(float(substitute(v, Dict()))) for v in vals]
    @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
end

@testset "Polygon Adjacency - Hexagon C₆" begin
    C6 = [0 1 0 0 0 1;
          1 0 1 0 0 0;
          0 1 0 1 0 0;
          0 0 1 0 1 0;
          0 0 0 1 0 1;
          1 0 0 0 1 0]
    
    polygon_n = SymbolicDiagonalization._is_polygon_adjacency(C6)
    @test polygon_n == 6
    
    vals, poly, λ = symbolic_eigenvalues(C6)
    @test length(vals) == 6
    
    # Hexagon eigenvalues: 2, 1, -1, -2, -1, 1
    vals_sorted = sort(real.(vals))
    @test vals_sorted ≈ [-2.0, -1.0, -1.0, 1.0, 1.0, 2.0] atol=1e-10
end

@testset "Polygon Adjacency - Larger Cycles" begin
    # C₇ and C₈ to verify formulas work for larger n
    for n in [7, 8, 10, 12]
        # Build Cₙ adjacency matrix
        Cn = zeros(Int, n, n)
        for i in 1:n
            Cn[i, mod1(i+1, n)] = 1
            Cn[i, mod1(i-1, n)] = 1
        end
        
        polygon_n = SymbolicDiagonalization._is_polygon_adjacency(Cn)
        @test polygon_n == n
        
        vals, poly, λ = symbolic_eigenvalues(Cn)
        @test length(vals) == n
        
        # Verify numerically
        numeric_eigs = eigvals(float.(Cn))
        computed_eigs = [real(float(substitute(v, Dict()))) for v in vals]
        @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
    end
end

@testset "Polygon Adjacency - Non-Examples" begin
    # Complete graph K₄ (not a cycle)
    K4 = [0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0]
    @test isnothing(SymbolicDiagonalization._is_polygon_adjacency(K4))
    
    # Path graph P₄ (not closed)
    P4 = [0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0]
    @test isnothing(SymbolicDiagonalization._is_polygon_adjacency(P4))
    
    # Weighted cycle (not 0-1)
    W4 = [0 2 0 2; 2 0 2 0; 0 2 0 2; 2 0 2 0]
    @test isnothing(SymbolicDiagonalization._is_polygon_adjacency(W4))
end

# ============================================================================
# Quaternion Group Q₈ - Quaternion Matrices
# ============================================================================

@testset "Quaternion Units" begin
    # Verify quaternion unit matrices satisfy i² = j² = k² = ijk = -1
    I2, Qi, Qj, Qk = SymbolicDiagonalization._quaternion_unit_matrices()
    
    @test Qi^2 ≈ -I2
    @test Qj^2 ≈ -I2
    @test Qk^2 ≈ -I2
    @test Qi*Qj*Qk ≈ -I2
    
    # Verify ij = k, jk = i, ki = j
    @test Qi*Qj ≈ Qk
    @test Qj*Qk ≈ Qi
    @test Qk*Qi ≈ Qj
end

@testset "Quaternion Matrix Detection" begin
    # 2×2 quaternion matrix: q = a + bi + cj + dk
    # Matrix form: [a+bi  c+di; -c+di  a-bi]
    
    # Numeric test
    Q_num = [1.0+2.0im 3.0+4.0im; -3.0+4.0im 1.0-2.0im]
    quat = SymbolicDiagonalization._is_quaternion_matrix(Q_num)
    @test !isnothing(quat)
    a, b, c, d = quat
    @test a ≈ 1.0
    @test b ≈ 2.0
    @test c ≈ 3.0
    @test d ≈ 4.0
    
    # Pure quaternion (a = 0)
    Q_pure = [2.0im 3.0+4.0im; -3.0+4.0im -2.0im]
    quat_pure = SymbolicDiagonalization._is_quaternion_matrix(Q_pure)
    @test !isnothing(quat_pure)
    @test quat_pure[1] ≈ 0.0  # a = 0
    
    # Real scalar (b = c = d = 0)
    Q_real = [5.0+0.0im 0.0+0.0im; 0.0+0.0im 5.0-0.0im]
    quat_real = SymbolicDiagonalization._is_quaternion_matrix(Q_real)
    @test !isnothing(quat_real)
    @test quat_real[1] ≈ 5.0
    @test all(≈(0.0), quat_real[2:4])
    
    # Non-quaternion matrix
    M_bad = [1.0+2.0im 3.0+4.0im; 5.0+6.0im 7.0-8.0im]  # Doesn't satisfy conjugate structure
    @test isnothing(SymbolicDiagonalization._is_quaternion_matrix(M_bad))
end

@testset "Quaternion Eigenvalues - Numeric" begin
    # Test eigenvalue formula: λ = a ± i√(b² + c² + d²)
    
    # General quaternion
    Q_num = [1.0+2.0im 3.0+4.0im; -3.0+4.0im 1.0-2.0im]
    vals, poly, λ = symbolic_eigenvalues(Q_num)
    @test length(vals) == 2
    
    # Verify against LinearAlgebra
    numeric_eigs = eigvals(Q_num)
    computed_eigs = [complex(float(v)) for v in vals]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
    @test sort(imag(numeric_eigs)) ≈ sort(imag(computed_eigs)) atol=1e-10
    
    # Pure imaginary eigenvalues (a = 0)
    Q_pure = [0.0+1.0im 0.0+0.0im; 0.0+0.0im 0.0-1.0im]  # Just i-component
    vals_pure, _, _ = symbolic_eigenvalues(Q_pure)
    @test all(v -> abs(real(v)) < 1e-10, vals_pure)  # Real part should be 0
    
    # Double eigenvalue (b = c = d = 0, real scalar)
    Q_scalar = [3.0+0.0im 0.0+0.0im; 0.0+0.0im 3.0-0.0im]
    vals_scalar, _, _ = symbolic_eigenvalues(Q_scalar)
    @test length(vals_scalar) == 2
    @test all(v -> abs(v - 3.0) < 1e-10, vals_scalar)  # Both eigenvalues = 3
end

@testset "Block Quaternion Matrix" begin
    # 4×4 block-diagonal with two quaternion blocks
    Q1 = [1.0+0.0im 0.0+1.0im; 0.0+1.0im 1.0-0.0im]  # q = 1 + k
    Q2 = [2.0+1.0im 1.0+0.0im; -1.0+0.0im 2.0-1.0im]  # q = 2 + i + j
    
    M = zeros(ComplexF64, 4, 4)
    M[1:2, 1:2] = Q1
    M[3:4, 3:4] = Q2
    
    # Verify detection
    quaternions = SymbolicDiagonalization._is_block_quaternion(M)
    @test !isnothing(quaternions)
    @test length(quaternions) == 2
    
    # Compute eigenvalues
    vals, poly, λ = symbolic_eigenvalues(M)
    @test length(vals) == 4
    
    # Verify numerically
    numeric_eigs = eigvals(M)
    computed_eigs = [complex(float(v)) for v in vals]
    @test sort(real(numeric_eigs)) ≈ sort(real(computed_eigs)) atol=1e-10
end

@testset "Quaternion Algebra Operations" begin
    # Test quaternion multiplication
    q1 = (1.0, 2.0, 3.0, 4.0)  # 1 + 2i + 3j + 4k
    q2 = (5.0, 6.0, 7.0, 8.0)  # 5 + 6i + 7j + 8k
    
    # Manual calculation: (1+2i+3j+4k)(5+6i+7j+8k)
    # a = 1*5 - 2*6 - 3*7 - 4*8 = 5 - 12 - 21 - 32 = -60
    # b = 1*6 + 2*5 + 3*8 - 4*7 = 6 + 10 + 24 - 28 = 12
    # c = 1*7 - 2*8 + 3*5 + 4*6 = 7 - 16 + 15 + 24 = 30
    # d = 1*8 + 2*7 - 3*6 + 4*5 = 8 + 14 - 18 + 20 = 24
    
    result = SymbolicDiagonalization._symbolic_quaternion_multiply(q1, q2)
    @test result[1] ≈ -60.0  # a
    @test result[2] ≈ 12.0   # b
    @test result[3] ≈ 30.0   # c
    @test result[4] ≈ 24.0   # d
    
    # Test norm squared
    norm_sq = SymbolicDiagonalization._symbolic_quaternion_norm_squared(1.0, 2.0, 3.0, 4.0)
    @test norm_sq ≈ 30.0  # 1 + 4 + 9 + 16
    
    # Test inverse: q * q⁻¹ = 1
    q_inv = SymbolicDiagonalization._symbolic_quaternion_inverse(1.0, 2.0, 3.0, 4.0)
    product = SymbolicDiagonalization._symbolic_quaternion_multiply(q1, q_inv)
    @test product[1] ≈ 1.0 atol=1e-10
    @test abs(product[2]) < 1e-10
    @test abs(product[3]) < 1e-10
    @test abs(product[4]) < 1e-10
end

@testset "SO(2) ⊗ SO(2) Kronecker Products" begin
    @variables θ φ
    
    @testset "Basic SO(2) ⊗ SO(2)" begin
        R1 = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        R2 = [cos(φ) -sin(φ); sin(φ) cos(φ)]
        R = kron(R1, R2)
        
        vals, _, _ = symbolic_eigenvalues(R)
        @test length(vals) == 4
        
        # Verify numerically
        θ_val, φ_val = 0.3, 0.7
        
        # Expected eigenvalues: e^{±i(θ±φ)}
        expected = [
            cos(θ_val + φ_val) + im*sin(θ_val + φ_val),
            cos(θ_val - φ_val) + im*sin(θ_val - φ_val),
            cos(-θ_val + φ_val) + im*sin(-θ_val + φ_val),
            cos(-θ_val - φ_val) + im*sin(-θ_val - φ_val),
        ]
        
        # Evaluate symbolic eigenvalues
        computed = ComplexF64[]
        for v in vals
            re_func = Symbolics.build_function(real(v), [θ, φ], expression=Val{false})
            im_func = Symbolics.build_function(imag(v), [θ, φ], expression=Val{false})
            push!(computed, ComplexF64(re_func([θ_val, φ_val]), im_func([θ_val, φ_val])))
        end
        
        @test isapprox(sort(computed, by=angle), sort(expected, by=angle), atol=1e-10)
    end
    
    @testset "Clean output without sqrt" begin
        R1 = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        R2 = [cos(φ) -sin(φ); sin(φ) cos(φ)]
        R = kron(R1, R2)
        
        vals, _, _ = symbolic_eigenvalues(R)
        
        # Check that no eigenvalue contains sqrt - clean trig expressions
        for v in vals
            str = string(v)
            @test !occursin("sqrt", str) 
        end
    end
    
    # Note: Same-angle case (θ = φ) is a known limitation.
    # When using the same symbolic variable for both rotations, 
    # _is_symbolic_orthogonal may fail because Symbolics doesn't 
    # simplify sin²(θ) + cos²(θ) = 1. The different-angle case works correctly.
end

# ============================================================================
# Coxeter/Weyl Groups - Reflection Groups
# ============================================================================

@testset "Cartan Matrix Type A" begin
    # Type Aₙ Cartan matrix: tridiagonal with 2 on diagonal, -1 off-diagonal
    for n in 1:6
        C = SymbolicDiagonalization._cartan_matrix_A(n)
        @test SymbolicDiagonalization._is_cartan_matrix_A(C) == n
        
        vals, _, _ = symbolic_eigenvalues(C)
        @test length(vals) == n
        
        # Verify eigenvalues numerically
        numeric_eigs = eigvals(Float64.(C))
        computed_eigs = [real(complex(v)) for v in vals]
        @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
    end
end

@testset "Cartan Matrix Type D" begin
    # Type Dₙ Cartan matrix: forked structure (constructor only, no symbolic eigenvalues)
    for n in 4:6
        C = SymbolicDiagonalization._cartan_matrix_D(n)
        @test SymbolicDiagonalization._is_cartan_matrix_D(C) == n
        @test SymbolicDiagonalization._detect_cartan_type(C) == (:D, n)
        # No symbolic eigenvalue support - matrix constructors only
    end
end

@testset "Cartan Matrix Type E" begin
    # Type Eₙ Cartan matrix: exceptional types (constructor only, no symbolic eigenvalues)
    for n in [6, 7, 8]
        C = SymbolicDiagonalization._cartan_matrix_E(n)
        @test SymbolicDiagonalization._is_cartan_matrix_E(C) == n
        @test SymbolicDiagonalization._detect_cartan_type(C) == (:E, n)
        # No symbolic eigenvalue support - matrix constructors only
    end
end

@testset "Cartan Matrix Type B" begin
    # Type Bₙ Cartan matrix: SO(2n+1), non-symmetric (constructor only)
    for n in [2, 3, 4, 5]
        C = SymbolicDiagonalization._cartan_matrix_B(n)
        @test SymbolicDiagonalization._is_cartan_matrix_B(C) == n
        @test SymbolicDiagonalization._detect_cartan_type(C) == (:B, n)
        
        # Verify asymmetric bond: C[1,2] = -2, C[2,1] = -1
        @test C[1, 2] == -2
        @test C[2, 1] == -1
        # No symbolic eigenvalue support
    end
end

@testset "Cartan Matrix Type C" begin
    # Type Cₙ Cartan matrix: Sp(2n), non-symmetric (constructor only)
    for n in [2, 3, 4, 5]
        C = SymbolicDiagonalization._cartan_matrix_C(n)
        @test SymbolicDiagonalization._is_cartan_matrix_C(C) == n
        @test SymbolicDiagonalization._detect_cartan_type(C) == (:C, n)
        
        # Verify asymmetric bond: C[n,n-1] = -2, C[n-1,n] = -1
        @test C[n, n-1] == -2
        @test C[n-1, n] == -1
        # No symbolic eigenvalue support
    end
end

@testset "Cartan Matrix Type F4" begin
    # Type F₄ Cartan matrix: exceptional, non-symmetric (constructor only)
    C = SymbolicDiagonalization._cartan_matrix_F4()
    @test size(C) == (4, 4)
    @test SymbolicDiagonalization._is_cartan_matrix_F4(C) == 4
    @test SymbolicDiagonalization._detect_cartan_type(C) == (:F, 4)
    
    # Verify double bond: C[2,3] = -2, C[3,2] = -1
    @test C[2, 3] == -2
    @test C[3, 2] == -1
    # No symbolic eigenvalue support
end

@testset "Cartan Matrix Type G2" begin
    # Type G₂ Cartan matrix: exceptional with SYMBOLIC eigenvalue formula
    C = SymbolicDiagonalization._cartan_matrix_G2()
    @test size(C) == (2, 2)
    @test SymbolicDiagonalization._is_cartan_matrix_G2(C) == 2
    @test SymbolicDiagonalization._detect_cartan_type(C) == (:G, 2)
    @test SymbolicDiagonalization._detect_cartan_type_symbolic(C) == (:G, 2)
    
    # Verify triple bond: C[1,2] = -3, C[2,1] = -1
    @test C[1, 2] == -3
    @test C[2, 1] == -1
    
    # Symbolic eigenvalues: 2 ± √3
    vals, _, _ = symbolic_eigenvalues(C)
    @test length(vals) == 2
    @test sort([Float64(v) for v in vals]) ≈ [2 - sqrt(3), 2 + sqrt(3)] atol=1e-10
end

@testset "Cartan Symbolic Detection" begin
    # Only types A and G₂ should be detected for symbolic eigenvalues
    A5 = SymbolicDiagonalization._cartan_matrix_A(5)
    @test SymbolicDiagonalization._detect_cartan_type_symbolic(A5) == (:A, 5)
    
    G2 = SymbolicDiagonalization._cartan_matrix_G2()
    @test SymbolicDiagonalization._detect_cartan_type_symbolic(G2) == (:G, 2)
    
    # Types B, C, D, E, F should NOT be detected for symbolic eigenvalues
    B3 = SymbolicDiagonalization._cartan_matrix_B(3)
    @test isnothing(SymbolicDiagonalization._detect_cartan_type_symbolic(B3))
    
    C3 = SymbolicDiagonalization._cartan_matrix_C(3)
    @test isnothing(SymbolicDiagonalization._detect_cartan_type_symbolic(C3))
    
    D4 = SymbolicDiagonalization._cartan_matrix_D(4)
    @test isnothing(SymbolicDiagonalization._detect_cartan_type_symbolic(D4))
    
    E6 = SymbolicDiagonalization._cartan_matrix_E(6)
    @test isnothing(SymbolicDiagonalization._detect_cartan_type_symbolic(E6))
    
    F4 = SymbolicDiagonalization._cartan_matrix_F4()
    @test isnothing(SymbolicDiagonalization._detect_cartan_type_symbolic(F4))
end

@testset "Cartan Eigenvalues Dispatcher" begin
    # Test _cartan_eigenvalues dispatcher for symbolic types only
    @test length(SymbolicDiagonalization._cartan_eigenvalues(:A, 5)) == 5
    @test length(SymbolicDiagonalization._cartan_eigenvalues(:G, 2)) == 2
    
    # Non-symbolic types should error
    @test_throws ErrorException SymbolicDiagonalization._cartan_eigenvalues(:B, 4)
    @test_throws ErrorException SymbolicDiagonalization._cartan_eigenvalues(:C, 4)
    @test_throws ErrorException SymbolicDiagonalization._cartan_eigenvalues(:D, 5)
    @test_throws ErrorException SymbolicDiagonalization._cartan_eigenvalues(:E, 6)
    @test_throws ErrorException SymbolicDiagonalization._cartan_eigenvalues(:F, 4)
end

@testset "Public API Cartan Constructors" begin
    # Test exported functions (constructors work for all types)
    @test cartan_matrix_A(3) == SymbolicDiagonalization._cartan_matrix_A(3)
    @test cartan_matrix_B(3) == SymbolicDiagonalization._cartan_matrix_B(3)
    @test cartan_matrix_C(3) == SymbolicDiagonalization._cartan_matrix_C(3)
    @test cartan_matrix_D(4) == SymbolicDiagonalization._cartan_matrix_D(4)
    @test cartan_matrix_E(6) == SymbolicDiagonalization._cartan_matrix_E(6)
    @test cartan_matrix_F4() == SymbolicDiagonalization._cartan_matrix_F4()
    @test cartan_matrix_G2() == SymbolicDiagonalization._cartan_matrix_G2()
    
    # Test Coxeter functions (these are just lookup tables, no numerics)
    @test coxeter_number(:B, 4) == 8
    @test coxeter_number(:C, 4) == 8
    @test coxeter_number(:F, 4) == 12
    @test coxeter_number(:G, 2) == 6
    
    @test coxeter_exponents(:B, 3) == [1, 3, 5]
    @test coxeter_exponents(:F, 4) == [1, 5, 7, 11]
    @test coxeter_exponents(:G, 2) == [1, 5]
    
    # Test graph Laplacians
    @test size(path_laplacian(5)) == (5, 5)
    @test size(cycle_laplacian(5)) == (5, 5)
end

@testset "Coxeter Numbers and Exponents" begin
    # Coxeter numbers
    @test SymbolicDiagonalization._coxeter_number(:A, 4) == 5
    @test SymbolicDiagonalization._coxeter_number(:D, 5) == 8
    @test SymbolicDiagonalization._coxeter_number(:E, 6) == 12
    @test SymbolicDiagonalization._coxeter_number(:E, 7) == 18
    @test SymbolicDiagonalization._coxeter_number(:E, 8) == 30
    
    # Exponents
    @test SymbolicDiagonalization._coxeter_exponents(:A, 4) == [1, 2, 3, 4]
    @test SymbolicDiagonalization._coxeter_exponents(:E, 6) == [1, 4, 5, 7, 8, 11]
    @test SymbolicDiagonalization._coxeter_exponents(:E, 8) == [1, 7, 11, 13, 17, 19, 23, 29]
end

@testset "Coxeter Element Eigenvalues" begin
    # Coxeter element eigenvalues should be on unit circle
    for (type, n) in [(:A, 3), (:A, 4), (:D, 4), (:E, 6)]
        h = SymbolicDiagonalization._coxeter_number(type, n)
        exps = SymbolicDiagonalization._coxeter_exponents(type, n)
        eigs = SymbolicDiagonalization._coxeter_element_eigenvalues(type, n)
        
        # All eigenvalues on unit circle
        @test all(abs.(eigs) .≈ 1.0)
        
        # Check eigenvalues are e^{2πi·m/h}
        expected = [exp(2π*im * m / h) for m in exps]
        @test sort(eigs, by=angle) ≈ sort(expected, by=angle) atol=1e-10
    end
end

@testset "Path Graph Laplacian" begin
    for n in [3, 5, 8, 10]
        # Build path Laplacian
        L = zeros(Int, n, n)
        for i in 1:n
            L[i, i] = (i == 1 || i == n) ? 1 : 2
            if i > 1
                L[i, i-1] = -1
            end
            if i < n
                L[i, i+1] = -1
            end
        end
        
        @test SymbolicDiagonalization._is_path_laplacian(L) == n
        
        vals, _, _ = symbolic_eigenvalues(L)
        @test length(vals) == n
        
        numeric_eigs = eigvals(Float64.(L))
        computed_eigs = [real(complex(v)) for v in vals]
        @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
    end
end

@testset "Cycle Graph Laplacian" begin
    for n in [4, 5, 6, 8]
        # Build cycle Laplacian
        L = zeros(Int, n, n)
        for i in 1:n
            L[i, i] = 2
            L[i, mod1(i+1, n)] = -1
            L[i, mod1(i-1, n)] = -1
        end
        
        @test SymbolicDiagonalization._is_cycle_laplacian(L) == n
        
        vals, _, _ = symbolic_eigenvalues(L)
        @test length(vals) == n
        
        numeric_eigs = eigvals(Float64.(L))
        computed_eigs = [real(complex(v)) for v in vals]
        @test sort(numeric_eigs) ≈ sort(computed_eigs) atol=1e-10
    end
end

@testset "Cartan Matrix Non-Examples" begin
    # Not a Cartan matrix (wrong diagonal)
    M1 = [1 -1; -1 1]
    @test isnothing(SymbolicDiagonalization._is_cartan_matrix_A(M1))
    
    # Not type A (wrong off-diagonal)
    M2 = [2 -2 0; -2 2 -2; 0 -2 2]
    @test isnothing(SymbolicDiagonalization._is_cartan_matrix_A(M2))
    
    # Generic symmetric matrix
    M3 = [2 1 0; 1 2 1; 0 1 2]
    @test isnothing(SymbolicDiagonalization._is_cartan_matrix_A(M3))
end

@testset "Reflection Matrix Detection" begin
    # Householder reflection
    v = [1.0, 0.0, 0.0]
    H = SymbolicDiagonalization._householder_reflection(v)
    @test SymbolicDiagonalization._is_reflection_matrix(H)
    
    # Eigenvalues: {1, 1, -1}
    eigs = SymbolicDiagonalization._reflection_eigenvalues(3)
    @test count(==(1), eigs) == 2
    @test count(==(-1), eigs) == 1
    
    # Non-reflection (rotation)
    θ = π/4
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    @test !SymbolicDiagonalization._is_reflection_matrix(R)
end

# ============================================================================
# Hadamard Matrices
# ============================================================================

@testset "Hadamard Matrices" begin
    @testset "Hadamard H₁ (2×2)" begin
        H = hadamard_matrix(1)
        @test size(H) == (2, 2)
        @test all(abs.(H) .== 1)
        
        vals, _, _ = symbolic_eigenvalues(H)
        @test length(vals) == 2
        
        # Eigenvalues should be ±√2
        vals_sorted = sort(Float64.(vals))
        @test isapprox(vals_sorted[1], -sqrt(2), atol=1e-10)
        @test isapprox(vals_sorted[2], sqrt(2), atol=1e-10)
    end
    
    @testset "Hadamard H₂ (4×4)" begin
        H = hadamard_matrix(2)
        @test size(H) == (4, 4)
        @test all(abs.(H) .== 1)
        
        vals, _, _ = symbolic_eigenvalues(H)
        @test length(vals) == 4
        
        # Eigenvalues should be ±2, each with multiplicity 2
        vals_sorted = sort(Float64.(vals))
        @test isapprox(vals_sorted[1], -2.0, atol=1e-10)
        @test isapprox(vals_sorted[2], -2.0, atol=1e-10)
        @test isapprox(vals_sorted[3], 2.0, atol=1e-10)
        @test isapprox(vals_sorted[4], 2.0, atol=1e-10)
    end
    
    @testset "Hadamard H₃ (8×8)" begin
        H = hadamard_matrix(3)
        @test size(H) == (8, 8)
        @test all(abs.(H) .== 1)
        
        vals, _, _ = symbolic_eigenvalues(H)
        @test length(vals) == 8
        
        # Eigenvalues should be ±√8 = ±2√2, each with multiplicity 4
        sqrt8 = sqrt(8)
        @test count(v -> isapprox(v, sqrt8, atol=1e-10), vals) == 4
        @test count(v -> isapprox(v, -sqrt8, atol=1e-10), vals) == 4
    end
    
    @testset "Hadamard Orthogonality" begin
        for n in 1:4
            H = hadamard_matrix(n)
            m = 2^n
            # H * H' = m * I
            @test H * H' ≈ m * I(m)
        end
    end
end

# ============================================================================
# DFT (Discrete Fourier Transform) Matrices
# ============================================================================

@testset "DFT Matrices" begin
    @testset "DFT F₄ (4×4)" begin
        F = dft_matrix(4)
        @test size(F) == (4, 4)
        
        vals, _, _ = symbolic_eigenvalues(F)
        @test length(vals) == 4
        
        # Expected multiplicities for n=4: (√n:2, -√n:1, i√n:0, -i√n:1)
        sqrtn = sqrt(4)
        m1 = count(v -> isapprox(v, sqrtn), vals)
        m_neg1 = count(v -> isapprox(v, -sqrtn), vals)
        m_i = count(v -> isapprox(v, sqrtn*im), vals)
        m_negi = count(v -> isapprox(v, -sqrtn*im), vals)
        @test (m1, m_neg1, m_i, m_negi) == (2, 1, 1, 0)  # positive omega convention
    end
    
    @testset "DFT eigenvalue multiplicities" begin
        # Test multiplicities for various sizes
        # Using positive omega convention: ω = e^{+2πi/n}
        expected_mults = Dict(
            1 => (1, 0, 0, 0),
            2 => (1, 1, 0, 0),
            3 => (1, 1, 1, 0),   # swapped i/-i from negative omega
            4 => (2, 1, 1, 0),   # swapped i/-i from negative omega
            5 => (2, 1, 1, 1),
            6 => (2, 2, 1, 1),
            7 => (2, 2, 2, 1),   # swapped i/-i from negative omega
            8 => (3, 2, 2, 1),   # swapped i/-i from negative omega
        )
        
        for n in 1:8
            F = dft_matrix(n)
            vals, _, _ = symbolic_eigenvalues(F)
            sqrtn = sqrt(n)
            
            m1 = count(v -> isapprox(v, sqrtn), vals)
            m_neg1 = count(v -> isapprox(v, -sqrtn), vals)
            m_i = count(v -> isapprox(v, sqrtn*im), vals)
            m_negi = count(v -> isapprox(v, -sqrtn*im), vals)
            
            @test (m1, m_neg1, m_i, m_negi) == expected_mults[n]
        end
    end
    
    @testset "DFT unitarity (normalized)" begin
        for n in [2, 3, 4, 5]
            F = dft_matrix(n, normalized=true)
            # Normalized DFT is unitary: F * F' = I
            @test F * F' ≈ I(n) atol=1e-10
        end
    end
    
    @testset "Normalized DFT eigenvalues" begin
        for n in [2, 3, 4, 5]
            F = dft_matrix(n, normalized=true)
            vals, _, _ = symbolic_eigenvalues(F)
            
            # Normalized DFT has eigenvalues that are 4th roots of unity: {1, -1, i, -i}
            for v in vals
                is_4th_root = any(isapprox(v, r, atol=1e-10) for r in [1, -1, im, -im])
                @test is_4th_root
            end
        end
    end
end

# ============================================================================
# Q₈ Regular Representation (Quaternion Group)
# ============================================================================

# ============================================================================
# Companion Matrices (Frobenius form)
# ============================================================================

@testset "Companion Matrices" begin
    @testset "Companion Matrix - Quadratic" begin
        # p(x) = x² - 5x + 6 = (x-2)(x-3)
        # Coefficients: a₀ = 6, a₁ = -5 (so [6, -5])
        C = companion_matrix([6, -5])
        @test size(C) == (2, 2)
        
        # Verify structure
        coeffs = SymbolicDiagonalization._is_companion_matrix(C)
        @test !isnothing(coeffs)
        @test coeffs == [6, -5]
        
        vals, _, _ = symbolic_eigenvalues(C)
        @test length(vals) == 2
        
        # Eigenvalues should be 2 and 3
        vals_sorted = sort(Float64.(vals))
        @test vals_sorted ≈ [2.0, 3.0] atol=1e-10
    end
    
    @testset "Companion Matrix - Cubic" begin
        # p(x) = x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)
        # Coefficients: [a₀, a₁, a₂] = [6, -11, 6] (negated from characteristic polynomial)
        # Wait, let's be careful: companion matrix for x³ + a₂x² + a₁x + a₀
        # has last column [-a₀, -a₁, -a₂]
        # For (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6:
        # a₂ = -6, a₁ = 11, a₀ = -6
        # So coeffs = [-6, 11, -6]
        C = companion_matrix([-6, 11, -6])
        @test size(C) == (3, 3)
        
        coeffs = SymbolicDiagonalization._is_companion_matrix(C)
        @test !isnothing(coeffs)
        
        vals, _, _ = symbolic_eigenvalues(C)
        @test length(vals) == 3
        
        # Eigenvalues should be 1, 2, 3
        vals_sorted = sort(Float64.(vals))
        @test vals_sorted ≈ [1.0, 2.0, 3.0] atol=1e-10
    end
    
    @testset "Companion Matrix - Quartic" begin
        # p(x) = x⁴ - 10x³ + 35x² - 50x + 24 = (x-1)(x-2)(x-3)(x-4)
        # Monic polynomial: x⁴ + a₃x³ + a₂x² + a₁x + a₀
        # where a₃ = -10, a₂ = 35, a₁ = -50, a₀ = 24
        # So coeffs = [a₀, a₁, a₂, a₃] = [24, -50, 35, -10]
        C = companion_matrix([24, -50, 35, -10])
        @test size(C) == (4, 4)
        
        coeffs = SymbolicDiagonalization._is_companion_matrix(C)
        @test !isnothing(coeffs)
        
        vals, _, _ = symbolic_eigenvalues(C)
        @test length(vals) == 4
        
        # Eigenvalues should be 1, 2, 3, 4
        vals_sorted = sort(Float64.(vals))
        @test vals_sorted ≈ [1.0, 2.0, 3.0, 4.0] atol=1e-10
    end
    
    @testset "Companion Matrix - Symbolic" begin
        @variables a b c
        # Monic polynomial x³ + ax² + bx + c
        C = companion_matrix([c, b, a])
        @test size(C) == (3, 3)
        
        coeffs = SymbolicDiagonalization._is_companion_matrix(C)
        @test !isnothing(coeffs)
        @test isequal(coeffs[1], c)
        @test isequal(coeffs[2], b)
        @test isequal(coeffs[3], a)
        
        vals, _, _ = symbolic_eigenvalues(C)
        @test length(vals) == 3
        
        # Numeric verification
        a_val, b_val, c_val = -6.0, 11.0, -6.0  # (x-1)(x-2)(x-3)
        C_num = companion_matrix([c_val, b_val, a_val])
        vals_num, _, _ = symbolic_eigenvalues(C_num)
        vals_sorted = sort(Float64.(vals_num))
        @test vals_sorted ≈ [1.0, 2.0, 3.0] atol=1e-10
    end
    
    @testset "Companion Matrix - Complex Roots" begin
        # p(x) = x² + 1 has roots ±i
        C = companion_matrix([1, 0])  # x² + 0x + 1
        @test size(C) == (2, 2)
        
        vals, _, _ = symbolic_eigenvalues(C)
        @test length(vals) == 2
        
        # Eigenvalues should be ±i
        vals_complex = complex.(vals)
        @test sort(imag.(vals_complex)) ≈ [-1.0, 1.0] atol=1e-10
        @test all(abs.(real.(vals_complex)) .< 1e-10)
    end
    
    @testset "Companion Matrix - Non-Example" begin
        # Not a companion matrix (wrong structure)
        @variables a b c d
        M1 = [a b; c d]  # General 2×2
        @test isnothing(SymbolicDiagonalization._is_companion_matrix(M1))
        
        # Wrong subdiagonal
        M2 = [0 0 -1; 2 0 0; 0 1 0]  # subdiagonal not all 1s
        @test isnothing(SymbolicDiagonalization._is_companion_matrix(M2))
        
        # Wrong zeros
        M3 = [0 1 -1; 1 0 0; 0 1 0]  # non-zero in wrong place
        @test isnothing(SymbolicDiagonalization._is_companion_matrix(M3))
    end
    
    @testset "Companion Matrix - Edge Cases" begin
        # 1×1 matrices are not companion matrices (need at least 2×2)
        M1x1 = [5]
        @test isnothing(SymbolicDiagonalization._is_companion_matrix(M1x1))
        
        # 2×2 identity-like structure
        M2 = [0 1; 1 0]  # swap matrix, not companion
        @test isnothing(SymbolicDiagonalization._is_companion_matrix(M2))
    end
end

@testset "Q₈ Regular Representation" begin
    @testset "Q₈ constructor and detection" begin
        # Test symbolic Q₈ matrix construction
        @variables c1 cm1 ci cmi cj cmj ck cmk
        
        M = Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)
        @test size(M) == (8, 8)
        
        # Verify detection works
        coeffs = SymbolicDiagonalization._is_Q8_regular_representation(M)
        @test !isnothing(coeffs)
        @test length(coeffs) == 8
    end
    
    @testset "Q₈ eigenvalue formula" begin
        # Test with symbolic coefficients
        @variables c1 cm1 ci cmi cj cmj ck cmk
        
        M = Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)
        vals, _, _ = symbolic_eigenvalues(M)
        @test length(vals) == 8
        
        # Should have 5 distinct eigenvalues (4 from 1D irreps + 1 from 2D with mult 4)
        unique_vals = unique(vals)
        @test length(unique_vals) == 5
        
        # The 2D irrep eigenvalue (c1 - cm1) should have multiplicity 4
        two_d_eigenval = c1 - cm1
        mult_2d = count(v -> isequal(Symbolics.simplify(v - two_d_eigenval), 0), vals)
        @test mult_2d == 4
    end
    
    @testset "Q₈ numeric verification" begin
        # Numeric test with random coefficients
        c1, cm1, ci, cmi, cj, cmj, ck, cmk = 1.0, 0.5, 0.3, -0.3, 0.2, -0.2, 0.1, -0.1
        
        M = Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)
        vals, _, _ = symbolic_eigenvalues(M)
        
        # Compare with direct eigenvalue computation
        direct_vals = sort(real.(eigvals(M)))
        our_vals = sort(Float64.(vals))
        
        @test isapprox(direct_vals, our_vals, atol=1e-10)
    end
    
    @testset "Q₈ character theory eigenvalues" begin
        # Verify the character-theoretic formulas
        c1, cm1, ci, cmi, cj, cmj, ck, cmk = 2.0, 1.0, 0.5, 0.5, 0.3, 0.3, 0.2, 0.2
        
        # Expected eigenvalues from character table:
        λ1 = c1 + cm1 + ci + cmi + cj + cmj + ck + cmk  # trivial rep
        λ2 = c1 + cm1 + ci + cmi - cj - cmj - ck - cmk  # sign on j,k
        λ3 = c1 + cm1 - ci - cmi + cj + cmj - ck - cmk  # sign on i,k
        λ4 = c1 + cm1 - ci - cmi - cj - cmj + ck + cmk  # sign on i,j
        λ5 = c1 - cm1  # 2D irrep (multiplicity 4)
        
        expected = sort([λ1, λ2, λ3, λ4, λ5, λ5, λ5, λ5])
        
        M = Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)
        vals, _, _ = symbolic_eigenvalues(M)
        # Handle potential small imaginary parts from floating point
        computed = sort(real.(complex.(vals)))
        
        @test isapprox(expected, computed, atol=1e-10)
    end
    
    @testset "Q₈ group algebra identities" begin
        # Test that certain special Q₈ matrices have known eigenstructure
        
        # Identity: all coeffs = 1 except cm1 = 0
        # This is the sum of all right multiplication operators
        M_sum = Q8_invariant_matrix(1, 0, 1, 1, 1, 1, 1, 1)
        vals_sum, _, _ = symbolic_eigenvalues(M_sum)
        # The 2D eigenvalue should be 1-0=1 with mult 4
        @test count(v -> isapprox(real(complex(v)), 1.0, atol=1e-10), vals_sum) == 4
        
        # Projection onto trivial: c_g = 1 for all g
        M_trivial = Q8_invariant_matrix(1, 1, 1, 1, 1, 1, 1, 1)
        vals_trivial, _, _ = symbolic_eigenvalues(M_trivial)
        # Trivial eigenvalue = 8, all others = 0
        @test any(v -> isapprox(real(complex(v)), 8.0, atol=1e-10), vals_trivial)
        @test count(v -> isapprox(abs(complex(v)), 0.0, atol=1e-10), vals_trivial) >= 4  # 2D irrep gives 0
    end
end
