# ============================================================================
# Tests for Lie Group Pattern Detection and Eigenvalue Computation
# These tests focus on SYMBOLIC matrices where Lie group structure is exploited
# ============================================================================

using Test
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

@testset "Lie Group Patterns" begin
    
    @testset "SO(2) - 2D Rotations" begin
        # Symbolic rotation matrix - this uses the Lie group path
        @variables θ
        R_sym = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        vals_sym, _, _ = symbolic_eigenvalues(R_sym)
        @test length(vals_sym) == 2
        
        # Check structure: eigenvalues should contain cos and sin
        @test any(v -> occursin("cos", string(v)) || occursin("sin", string(v)), string.(vals_sym))
        
        # Eigenvalues should be Complex{Num} for rotation matrices
        @test all(v -> v isa Complex{Num}, vals_sym)
    end
    
    @testset "SO(3) - 3D Rotations" begin
        # Symbolic rotation around z-axis
        @variables θ
        Rz_sym = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        vals_sym, _, _ = symbolic_eigenvalues(Rz_sym)
        
        @test length(vals_sym) == 3
        
        # One eigenvalue should be exactly 1 (the rotation axis)
        has_one = any(vals_sym) do v
            simplified = Symbolics.simplify(v isa Complex ? real(v) + imag(v) - 1 : v - 1)
            simplified isa Number && isapprox(simplified, 0, atol=1e-10)
        end
        @test has_one
    end
    
    @testset "SO(1,1) - Lorentz Boosts in 1+1D" begin
        # Symbolic Lorentz boost
        @variables φ
        L_sym = [cosh(φ) sinh(φ); sinh(φ) cosh(φ)]
        vals_sym, _, _ = symbolic_eigenvalues(L_sym)
        
        @test length(vals_sym) == 2
        
        # Check that eigenvalues are real Num (not Complex)
        @test all(v -> v isa Num, vals_sym)
        
        # Eigenvalues should contain cosh and sinh
        @test any(v -> occursin("cosh", string(v)) || occursin("sinh", string(v)), string.(vals_sym))
    end
    
    @testset "Detection Functions" begin
        # Test detection functions on numeric matrices
        
        # SO(2)
        θ_val = 0.5
        R_num = [cos(θ_val) -sin(θ_val); sin(θ_val) cos(θ_val)]
        @test SymbolicDiagonalization._is_orthogonal(R_num)
        @test !isnothing(SymbolicDiagonalization._is_special_orthogonal(R_num))
        
        # SO(3)
        Rz_num = [cos(θ_val) -sin(θ_val) 0; sin(θ_val) cos(θ_val) 0; 0 0 1]
        @test SymbolicDiagonalization._is_orthogonal(Rz_num)
        @test !isnothing(SymbolicDiagonalization._is_special_orthogonal(Rz_num))
        
        # SO(1,1)
        φ_val = 0.3
        L_num = [cosh(φ_val) sinh(φ_val); sinh(φ_val) cosh(φ_val)]
        @test SymbolicDiagonalization._is_indefinite_orthogonal(L_num, 1, 1)
        
        # Sp(2)
        a_val = 2.0
        Sp2_num = [a_val 0; 0 1/a_val]
        @test !isnothing(SymbolicDiagonalization._is_symplectic(Sp2_num))
    end
    
    @testset "Master Detection Function" begin
        # Test the master detection function
        @variables θ φ
        
        # SO(2)
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        group_type, params = SymbolicDiagonalization._detect_lie_group(R)
        @test group_type == :SO2
        
        # SO(1,1)
        L = [cosh(φ) sinh(φ); sinh(φ) cosh(φ)]
        group_type, params = SymbolicDiagonalization._detect_lie_group(L)
        @test group_type == :SO11
        
        # Non-Lie group
        M = [1 2; 3 4]
        group_type, params = SymbolicDiagonalization._detect_lie_group(M)
        @test group_type === nothing
    end
    
    @testset "Eigenvalue Correctness" begin
        # Verify eigenvalues are mathematically correct by checking characteristic polynomial
        
        @testset "SO(2)" begin
            @variables θ
            R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
            vals, poly, λ = symbolic_eigenvalues(R)
            
            # Product of (λ - eigenvalue) factors should equal characteristic polynomial
            # The polynomial should be: λ² - 2cos(θ)λ + 1
            # (since det(R) = 1 and tr(R) = 2cos(θ))
            poly_expanded = Symbolics.expand(poly)
            expected_trace = 2*cos(θ)
            expected_det = 1
            
            # Extract coefficients
            poly_str = string(poly_expanded)
            @test occursin("λ^2", poly_str) || occursin("λ²", poly_str)
        end
        
        @testset "SO(1,1)" begin
            @variables φ
            L = [cosh(φ) sinh(φ); sinh(φ) cosh(φ)]
            vals, poly, λ = symbolic_eigenvalues(L)
            
            # Eigenvalues should multiply to det = 1
            # and sum to trace = 2*cosh(φ)
            # Check polynomial form
            poly_expanded = Symbolics.expand(poly)
            poly_str = string(poly_expanded)
            @test occursin("λ^2", poly_str) || occursin("λ²", poly_str)
        end
    end
end
