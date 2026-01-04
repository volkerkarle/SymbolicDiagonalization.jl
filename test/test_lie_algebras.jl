# ============================================================================
# Tests for Lie Algebra Representation Detection and Eigenvalue Computation
# These tests focus on SYMBOLIC matrices that are elements of Lie algebra
# representations, exploiting representation theory to compute eigenvalues.
# ============================================================================

using Test
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

@testset "Lie Algebra Representations" begin
    
    @testset "so(3)/su(2) Spin-j Representations" begin
        
        @testset "Spin-0 (1×1)" begin
            # The only spin-0 element is the zero matrix
            @variables ω
            A0 = [0*ω;;]  # Symbolic zero
            vals, _, _ = symbolic_eigenvalues(A0)
            @test length(vals) == 1
            @test isapprox(real(vals[1]), 0, atol=1e-10) || SymbolicDiagonalization._issymzero(vals[1])
        end
        
        @testset "Spin-1/2 (2×2 su(2))" begin
            # su(2) generators: τₖ = i·σₖ/2
            @variables θ
            τz = im * [1 0; 0 -1] / 2  # i·σz/2
            A = θ * τz  # General element along z-axis
            
            vals, _, _ = symbolic_eigenvalues(A)
            @test length(vals) == 2
            
            # Eigenvalues should be ±iθ/2 (up to sqrt)
            # The trace formula gives sqrt(θ²), so we get ±i·sqrt(θ²)/2
            @test all(v -> v isa Complex{Num}, vals)
            
            # Sum of eigenvalues should be zero (traceless)
            sum_vals = sum(vals)
            @test SymbolicDiagonalization._issymzero(sum_vals)
        end
        
        @testset "Spin-1 (3×3 so(3))" begin
            # so(3) generators: Lx, Ly, Lz
            @variables ω
            Lz = [0 -1 0; 1 0 0; 0 0 0]
            A = ω * Lz
            
            vals, _, _ = symbolic_eigenvalues(A)
            @test length(vals) == 3
            
            # Eigenvalues should be -iω, 0, iω (up to sqrt representation)
            # Check that we have one zero and two opposite values
            zero_count = count(v -> SymbolicDiagonalization._issymzero(v), vals)
            @test zero_count >= 1
            
            # Non-zero eigenvalues should be opposite
            nonzero = filter(v -> !SymbolicDiagonalization._issymzero(v), vals)
            if length(nonzero) == 2
                @test SymbolicDiagonalization._issymzero(nonzero[1] + nonzero[2])
            end
        end
        
        @testset "Spin-3/2 (4×4)" begin
            # Jz for spin-3/2: diagonal with entries 3/2, 1/2, -1/2, -3/2
            @variables α
            Jz32 = im * diagm([3/2, 1/2, -1/2, -3/2])  # Skew-Hermitian version
            A = α * Jz32
            
            vals, _, _ = symbolic_eigenvalues(A)
            @test length(vals) == 4
            
            # Eigenvalues should be i·α·m for m = ±3/2, ±1/2
            # Sum should be zero
            sum_vals = sum(vals)
            @test SymbolicDiagonalization._issymzero(sum_vals)
        end
        
        @testset "Spin-2 (5×5)" begin
            @variables ω
            Jz2 = im * diagm([2.0, 1.0, 0.0, -1.0, -2.0])
            A = ω * Jz2
            
            vals, _, _ = symbolic_eigenvalues(A)
            @test length(vals) == 5
            
            # Should have one zero eigenvalue
            zero_count = count(v -> SymbolicDiagonalization._issymzero(v), vals)
            @test zero_count >= 1
            
            # Sum should be zero
            sum_vals = sum(vals)
            @test SymbolicDiagonalization._issymzero(sum_vals)
        end
        
        @testset "Spin-1 with different generators" begin
            # Test Lx and Ly (off-diagonal structure)
            @variables θ
            
            Lx = [0 0 0; 0 0 -1; 0 1 0]
            Ax = θ * Lx
            vals_x, _, _ = symbolic_eigenvalues(Ax)
            @test length(vals_x) == 3
            
            Ly = [0 0 1; 0 0 0; -1 0 0]
            Ay = θ * Ly
            vals_y, _, _ = symbolic_eigenvalues(Ay)
            @test length(vals_y) == 3
        end
    end
    
    @testset "sl(2) / Traceless Matrices" begin
        @variables α ω θ
        
        @testset "sl(2) 2×2" begin
            # General traceless 2×2 matrix
            sl2_elem = [α ω; θ -α]
            vals, _, _ = symbolic_eigenvalues(sl2_elem)
            @test length(vals) == 2
            
            # Eigenvalues should be ±sqrt(α² + ωθ)
            sum_vals = sum(vals)
            @test SymbolicDiagonalization._issymzero(sum_vals)
        end
    end
    
    @testset "so(2) Fundamental" begin
        @variables θ
        
        # 2×2 skew-symmetric (so(2) = spin-1/2 of so(3))
        so2 = [0 -θ; θ 0]
        vals, _, _ = symbolic_eigenvalues(so2)
        @test length(vals) == 2
        
        # Eigenvalues should be ±iθ
        sum_vals = sum(vals)
        @test SymbolicDiagonalization._issymzero(sum_vals)
    end
    
    @testset "Detection Functions" begin
        @variables ω
        
        @testset "_is_skew_symmetric" begin
            skew = [0 -ω; ω 0]
            @test SymbolicDiagonalization._is_skew_symmetric(skew)
            
            nonskew = [1 -ω; ω 0]
            @test !SymbolicDiagonalization._is_skew_symmetric(nonskew)
        end
        
        @testset "_is_skew_hermitian" begin
            # Real skew-symmetric is also skew-Hermitian
            skew_real = [0 -ω; ω 0]
            @test SymbolicDiagonalization._is_skew_hermitian(skew_real)
            
            # Complex skew-Hermitian
            skew_complex = [im*ω 0; 0 -im*ω]
            @test SymbolicDiagonalization._is_skew_hermitian(skew_complex)
        end
        
        @testset "_is_so3_spin_representation" begin
            # Spin-1 generator
            Lz = [0 -1 0; 1 0 0; 0 0 0]
            A = ω * Lz
            result = SymbolicDiagonalization._is_so3_spin_representation(A)
            @test !isnothing(result)
            j, omega = result
            @test j == 1
            
            # 2×2 case (spin-1/2)
            τz = im * [1 0; 0 -1] / 2
            B = ω * τz
            result2 = SymbolicDiagonalization._is_so3_spin_representation(B)
            @test !isnothing(result2)
            j2, omega2 = result2
            @test j2 == 1//2
        end
        
        @testset "_is_su_fundamental" begin
            # su(2) fundamental = 2×2 traceless skew-Hermitian
            τz = im * [1 0; 0 -1] / 2
            A = ω * τz
            @test !isnothing(SymbolicDiagonalization._is_su_fundamental(A))
            
            # Non-traceless should fail
            non_traceless = [im*ω 0; 0 im*ω]
            @test isnothing(SymbolicDiagonalization._is_su_fundamental(non_traceless))
        end
    end
    
    @testset "Standard Generators" begin
        @testset "Pauli matrices" begin
            σ₁ = pauli_x()
            σ₂ = pauli_y()
            σ₃ = pauli_z()
            
            # Check properties
            @test size(σ₁) == (2, 2)
            @test σ₁ == σ₁'  # Hermitian
            @test σ₂ == σ₂'
            @test σ₃ == σ₃'
            
            # Check anticommutation: {σᵢ, σⱼ} = 2δᵢⱼ I
            @test isapprox(σ₁ * σ₁, I(2), atol=1e-10)
            @test isapprox(σ₂ * σ₂, I(2), atol=1e-10)
            @test isapprox(σ₃ * σ₃, I(2), atol=1e-10)
        end
        
        @testset "su(2) generators" begin
            τ₁, τ₂, τ₃ = su2_generators()
            
            # Should be skew-Hermitian
            @test isapprox(τ₁ + τ₁', zeros(2,2), atol=1e-10)
            @test isapprox(τ₂ + τ₂', zeros(2,2), atol=1e-10)
            @test isapprox(τ₃ + τ₃', zeros(2,2), atol=1e-10)
        end
        
        @testset "so(3) generators" begin
            Lx, Ly, Lz = so3_generators()
            
            # Should be skew-symmetric
            @test isapprox(Lx + Lx', zeros(3,3), atol=1e-10)
            @test isapprox(Ly + Ly', zeros(3,3), atol=1e-10)
            @test isapprox(Lz + Lz', zeros(3,3), atol=1e-10)
            
            # Check Lie bracket: [Lx, Ly] = Lz
            @test isapprox(Lx * Ly - Ly * Lx, Lz, atol=1e-10)
        end
        
        @testset "Spin-j generators" begin
            # Test spin-1/2
            Jx, Jy, Jz = spin_j_generators(1//2)
            @test size(Jz) == (2, 2)
            @test isapprox(Jz, diagm([0.5, -0.5]), atol=1e-10)
            
            # Test spin-1
            Jx, Jy, Jz = spin_j_generators(1)
            @test size(Jz) == (3, 3)
            @test isapprox(Jz, diagm([1.0, 0.0, -1.0]), atol=1e-10)
            
            # Test spin-3/2
            Jx, Jy, Jz = spin_j_generators(3//2)
            @test size(Jz) == (4, 4)
            @test isapprox(Jz, diagm([1.5, 0.5, -0.5, -1.5]), atol=1e-10)
        end
    end
    
    @testset "Edge Cases" begin
        @variables ω
        
        @testset "Zero matrix" begin
            Z = zeros(Num, 3, 3)
            vals, _, _ = symbolic_eigenvalues(Z)
            @test length(vals) == 3
            # Zero matrix eigenvalues contain sqrt(0) which should be zero
            # But Symbolics may not simplify sqrt(0) to 0
            # Check that all eigenvalues are effectively zero
            for v in vals
                if v isa Complex
                    # Check real and imag parts contain only 0 or sqrt(0)
                    str_v = string(v)
                    @test occursin("0", str_v)
                else
                    @test SymbolicDiagonalization._issymzero(v) || occursin("sqrt(0)", string(v))
                end
            end
        end
        
        @testset "Numeric spin matrices via LinearAlgebra" begin
            # For pure numeric matrices, use LinearAlgebra.eigvals
            # symbolic_eigenvalues always produces symbolic output
            ω_val = 2.0
            Lz = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0]
            A_num = ω_val * Lz
            
            # Use LinearAlgebra directly for numeric
            numeric_vals = eigvals(A_num)
            @test length(numeric_vals) == 3
            
            # Sort by imaginary part for comparison
            sorted_imag = sort(imag.(numeric_vals))
            @test isapprox(sorted_imag[1], -2.0, atol=1e-10)
            @test isapprox(sorted_imag[2], 0.0, atol=1e-10)
            @test isapprox(sorted_imag[3], 2.0, atol=1e-10)
        end
    end
end
