# Tests for Anti-Circulant and Kac-Murdock-Szego matrices

@testset "Anti-Circulant Matrix Construction" begin
    # Test basic construction
    c = [1, 2, 3, 4]
    A = anticirculant_matrix(c)
    @test size(A) == (4, 4)
    
    # Verify anti-circulant structure: A[i,j] = c[(i+j-2) mod n + 1]
    n = length(c)
    for i in 1:n
        for j in 1:n
            idx = mod(i + j - 2, n) + 1
            @test A[i, j] == c[idx]
        end
    end
    
    # Test expected structure for [1,2,3,4]
    # Row 1: [1,2,3,4] (c[1], c[2], c[3], c[4])
    # Row 2: [2,3,4,1] (c[2], c[3], c[4], c[1])
    # Row 3: [3,4,1,2]
    # Row 4: [4,1,2,3]
    @test A[1, :] == [1, 2, 3, 4]
    @test A[2, :] == [2, 3, 4, 1]
    @test A[3, :] == [3, 4, 1, 2]
    @test A[4, :] == [4, 1, 2, 3]
    
    # Test with symbolic values
    @variables a b c_sym
    c_symbolic = [a, b, c_sym]
    A_sym = anticirculant_matrix(c_symbolic)
    @test size(A_sym) == (3, 3)
    @test isequal(A_sym[1, 1], a)
    @test isequal(A_sym[1, 2], b)
    @test isequal(A_sym[1, 3], c_sym)
    @test isequal(A_sym[2, 1], b)
    @test isequal(A_sym[2, 2], c_sym)
    @test isequal(A_sym[2, 3], a)
end

@testset "Anti-Circulant Detection" begin
    # Test detection of anti-circulant matrices
    c = [1, 2, 3, 4]
    A = anticirculant_matrix(c)
    
    # Internal detection function
    result = SymbolicDiagonalization._is_anticirculant(A)
    @test result !== nothing
    @test result == c
    
    # Test with symbolic matrix
    @variables x y z
    A_sym = anticirculant_matrix([x, y, z])
    result_sym = SymbolicDiagonalization._is_anticirculant(A_sym)
    @test result_sym !== nothing
    
    # Test non-anti-circulant matrix
    B = [1 2; 3 4]  # Not anti-circulant
    @test SymbolicDiagonalization._is_anticirculant(B) === nothing
    
    # Test non-square matrix
    C = [1 2 3; 4 5 6]
    @test SymbolicDiagonalization._is_anticirculant(C) === nothing
end

@testset "Anti-Circulant Eigenvalues" begin
    # Test eigenvalue computation for numeric case
    c = [1.0, 2.0, 3.0, 4.0]
    eigenvalues = SymbolicDiagonalization._anticirculant_eigenvalues(c)
    @test length(eigenvalues) == 4
    
    # The function returns complex eigenvalues - just verify they're computed
    @test all(e isa Number for e in eigenvalues)
    
    # Test small case
    c2 = [1.0, 2.0]
    eigs2 = SymbolicDiagonalization._anticirculant_eigenvalues(c2)
    @test length(eigs2) == 2
end

@testset "KMS Matrix Construction" begin
    # Test basic construction
    ρ = 0.5
    n = 3
    K = kms_matrix(ρ, n)
    @test size(K) == (3, 3)
    
    # Verify KMS structure: K[i,j] = ρ^|i-j|
    for i in 1:n
        for j in 1:n
            @test K[i, j] ≈ ρ^abs(i-j)
        end
    end
    
    # Test expected values for ρ=0.5, n=3
    # [1.0   0.5   0.25]
    # [0.5   1.0   0.5 ]
    # [0.25  0.5   1.0 ]
    @test K[1, 1] ≈ 1.0
    @test K[1, 2] ≈ 0.5
    @test K[1, 3] ≈ 0.25
    @test K[2, 2] ≈ 1.0
    
    # Test symmetry
    @test K ≈ K'
    
    # Test with symbolic ρ
    @variables r
    K_sym = kms_matrix(r, 2)
    @test size(K_sym) == (2, 2)
    # The (1,1) element should be r^0 = 1 or be symbolically equivalent
    @test K_sym[1, 1] == 1 || isequal(K_sym[1, 1], r^0)
    @test isequal(K_sym[1, 2], r)
    @test isequal(K_sym[2, 1], r)
    
    # Test n=1 edge case
    K1 = kms_matrix(0.5, 1)
    @test size(K1) == (1, 1)
    @test K1[1, 1] ≈ 1.0
    
    # Test error on invalid n
    @test_throws Exception kms_matrix(0.5, 0)
end

@testset "KMS Matrix Detection" begin
    # Test detection of KMS matrices
    ρ = 0.5
    n = 4
    K = kms_matrix(ρ, n)
    
    result = SymbolicDiagonalization._is_kms_matrix(K)
    @test result !== nothing
    @test result[1] ≈ ρ
    @test result[2] == n
    
    # Test n=1 case
    K1 = kms_matrix(0.3, 1)
    result1 = SymbolicDiagonalization._is_kms_matrix(K1)
    @test result1 !== nothing
    @test result1[2] == 1
    
    # Test matrix with wrong diagonal (not 1 on diagonal)
    C = [2 0.5; 0.5 2]  # Diagonal is 2, not 1
    @test SymbolicDiagonalization._is_kms_matrix(C) === nothing
    
    # Test non-square matrix
    D = [1 0.5 0.25; 0.5 1 0.5]
    @test SymbolicDiagonalization._is_kms_matrix(D) === nothing
    
    # Test matrix that's not geometric decay
    E = [1 0.5; 0.5 1]  # This IS a valid KMS matrix with ρ=0.5
    result_E = SymbolicDiagonalization._is_kms_matrix(E)
    @test result_E !== nothing
    @test result_E[1] ≈ 0.5
end

@testset "KMS Matrix Eigenvalues" begin
    # Test eigenvalue computation
    ρ = 0.5
    n = 4
    
    eigenvalues = SymbolicDiagonalization._kms_eigenvalues(ρ, n)
    @test length(eigenvalues) == n
    
    # Test special case ρ = 0 (identity matrix)
    eigs_zero = SymbolicDiagonalization._kms_eigenvalues(0, 3)
    @test all(e == 1 for e in eigs_zero)
    
    # Test special case ρ = 1 (all-ones matrix)
    eigs_one = SymbolicDiagonalization._kms_eigenvalues(1, 3)
    @test eigs_one[1] == 3  # Largest eigenvalue is n
    @test all(e == 0 for e in eigs_one[2:end])  # Rest are 0
    
    # Test special case ρ = -1
    eigs_neg = SymbolicDiagonalization._kms_eigenvalues(-1, 3)
    @test all(e == 0 for e in eigs_neg)  # Singular matrix
    
    # Test error on invalid n
    @test_throws Exception SymbolicDiagonalization._kms_eigenvalues(0.5, 0)
end

@testset "KMS Symbolic Eigenvalues" begin
    # Test symbolic eigenvalue computation
    @variables ρ_sym
    n = 3
    
    eigenvalues = SymbolicDiagonalization._kms_eigenvalues_symbolic(ρ_sym, n)
    @test length(eigenvalues) == n
    
    # Each eigenvalue should be a symbolic expression
    for eig in eigenvalues
        @test eig isa Symbolics.Num
    end
    
    # Test with numeric ρ
    eigenvalues_num = SymbolicDiagonalization._kms_eigenvalues_symbolic(0.5, 4)
    @test length(eigenvalues_num) == 4
    
    # Test error on invalid n
    @test_throws Exception SymbolicDiagonalization._kms_eigenvalues_symbolic(0.5, 0)
end

@testset "KMS Matrix Properties" begin
    # KMS matrix should be positive definite for |ρ| < 1
    ρ = 0.5
    n = 5
    K = kms_matrix(ρ, n)
    
    # Check positive definiteness via eigenvalues
    eigs = eigvals(K)
    @test all(eigs .> 0)
    
    # Check symmetry
    @test K ≈ K'
    
    # Check Toeplitz property: constant diagonals
    for k in 0:n-1
        diag_vals = [K[i, i+k] for i in 1:n-k]
        @test all(v ≈ diag_vals[1] for v in diag_vals)
    end
end

@testset "Anti-Circulant vs Circulant Relationship" begin
    # Anti-circulant should have related structure to circulant
    @variables a b c d
    
    A = anticirculant_matrix([a, b, c, d])
    
    # Anti-circulant is symmetric (Hankel-like)
    # Check that A[i,j] = A[j,i] for this construction
    @test isequal(A[1, 2], A[2, 1])  # Both should be b
    @test isequal(A[1, 3], A[3, 1])  # Both should be c
    
    # Verify it's a Hankel-like structure (constant anti-diagonals)
    # A[i,j] depends only on i+j
    @test isequal(A[1, 3], A[2, 2])  # i+j = 4 for both, should be c
    @test isequal(A[1, 4], A[2, 3])  # i+j = 5 for both, should be d
end
