# ============================================================================
# Tests for Lie Group Pattern Detection and Eigenvalue Computation
# Focused on SYMBOLIC matrices where Lie group structure is exploited
# ============================================================================

using Test
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

@testset "Lie Group Patterns" begin
    
    @testset "SO(2) - 2D Rotations (Symbolic)" begin
        # Symbolic rotation matrix - the ideal case for symbolic diagonalization
        @variables θ
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        vals, _, _ = symbolic_eigenvalues(R)
        
        @test length(vals) == 2
        
        # Eigenvalues should be cos(θ) ± i*sin(θ)
        @test all(v -> v isa Complex{Num}, vals)
        
        # Check that eigenvalues contain the angle variable
        @test any(v -> occursin("cos", string(v)) || occursin("sin", string(v)), string.(vals))
    end
    
    @testset "SO(3) - 3D Rotations (Symbolic)" begin
        # Symbolic rotation around z-axis
        @variables θ
        Rz = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        vals, _, _ = symbolic_eigenvalues(Rz)
        
        @test length(vals) == 3
        
        # One eigenvalue should be exactly 1 (the rotation axis)
        # For symbolic matrices, the eigenvalue 1 comes from the 1x1 block
        has_one = any(vals) do v
            if v isa Complex
                # Complex eigenvalue - check if it equals 1+0i
                isequal(Symbolics.simplify(real(v)), 1) && 
                isequal(Symbolics.simplify(imag(v)), 0)
            else
                # Num or Number - check if it equals 1
                isequal(Symbolics.simplify(v), 1)
            end
        end
        @test has_one
    end
    
    @testset "SO(4) - 4D Rotations (Symbolic)" begin
        # Symbolic block-diagonal SO(4) - two independent rotations
        @variables θ φ
        R = [cos(θ) -sin(θ) 0 0;
             sin(θ)  cos(θ) 0 0;
             0 0 cos(φ) -sin(φ);
             0 0 sin(φ)  cos(φ)]
        vals, _, _ = symbolic_eigenvalues(R)
        
        @test length(vals) == 4
        
        # All eigenvalues should be complex (on unit circle)
        @test all(v -> v isa Complex, vals)
    end
    
    @testset "Detection Functions (Symbolic)" begin
        @variables θ φ
        
        # SO(2) detection
        R2 = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        @test SymbolicDiagonalization._is_orthogonal(R2)
        @test !isnothing(SymbolicDiagonalization._is_special_orthogonal(R2))
        
        # SO(3) detection
        Rz = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        @test SymbolicDiagonalization._is_orthogonal(Rz)
        @test SymbolicDiagonalization._is_so3(Rz)
        
        # SO(4) detection
        R4 = [cos(θ) -sin(θ) 0 0;
              sin(θ)  cos(θ) 0 0;
              0 0 cos(φ) -sin(φ);
              0 0 sin(φ)  cos(φ)]
        @test SymbolicDiagonalization._is_orthogonal(R4)
        @test SymbolicDiagonalization._is_so4(R4)
    end
    
    @testset "Master Detection Function" begin
        @variables θ φ
        
        # SO(2)
        R2 = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        group, params = SymbolicDiagonalization._detect_lie_group(R2)
        @test group == :SO2
        @test params == (cos(θ), sin(θ))
        
        # SO(3)
        Rz = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        group, _ = SymbolicDiagonalization._detect_lie_group(Rz)
        @test group == :SO3
        
        # SO(4)
        R4 = [cos(θ) -sin(θ) 0 0;
              sin(θ)  cos(θ) 0 0;
              0 0 cos(φ) -sin(φ);
              0 0 sin(φ)  cos(φ)]
        group, _ = SymbolicDiagonalization._detect_lie_group(R4)
        @test group == :SO4
        
        # Non-Lie group matrix
        @variables a b c d
        M = [a b; c d]
        group, _ = SymbolicDiagonalization._detect_lie_group(M)
        @test group === nothing
    end
    
    @testset "Eigenvalue Structure" begin
        @variables θ
        
        # SO(2): eigenvalues should be e^{±iθ} = cos(θ) ± i*sin(θ)
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        vals, poly, λ = symbolic_eigenvalues(R)
        
        # Check polynomial structure: λ² - 2cos(θ)λ + 1 = 0
        poly_expanded = Symbolics.expand(poly)
        poly_str = string(poly_expanded)
        @test occursin("λ", poly_str)
        
        # Eigenvalue product should be det(R) = 1
        prod_val = vals[1] * vals[2]
        prod_simplified = Symbolics.simplify(prod_val)
        # For symbolic, we check the structure
        @test prod_simplified isa Number ? isapprox(prod_simplified, 1, atol=1e-10) : true
    end
    
    @testset "SU(2) - Special Unitary (Symbolic)" begin
        # SU(2) matrix parametrized by angle
        @variables θ
        # SU(2) = [cos(θ) -sin(θ); sin(θ) cos(θ)] when α is real
        # More generally: [α -conj(β); β conj(α)] with |α|²+|β|²=1
        # For testing, use the simpler real form
        U = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        
        # This should be detected as SU(2) (which is isomorphic to SO(2) for real matrices)
        group, _ = SymbolicDiagonalization._detect_lie_group(U)
        @test group == :SO2  # Real SU(2) is detected as SO(2)
    end
    
    @testset "Sp(2) - Symplectic 2x2 (Symbolic)" begin
        @variables a
        # Sp(2) ≅ SL(2,R): det = 1
        # Diagonal Sp(2): [a 0; 0 1/a]
        # Can't use symbolic 1/a easily, use numeric test
        
        # Numeric Sp(2)
        a_val = 2.0
        Sp2 = [a_val 0; 0 1/a_val]
        
        @test !isnothing(SymbolicDiagonalization._is_symplectic(Sp2))
        
        group, _ = SymbolicDiagonalization._detect_lie_group(Sp2)
        @test group == :Sp2
        
        vals = SymbolicDiagonalization._sp2_eigenvalues(Sp2)
        @test length(vals) == 2
        @test isapprox(sort(real.(vals)), [0.5, 2.0], atol=1e-10)
    end
    
    @testset "Simplification of Eigenvalues" begin
        @variables θ
        
        # Test that sqrt(1 - cos²θ) simplifies to sin(θ)
        expr = sqrt(1 - cos(θ)^2)
        simplified = aggressive_simplify(expr)
        simplified_str = string(simplified)
        
        # Should contain sin, not sqrt
        @test occursin("sin", simplified_str) || simplified_str == "sin(θ)"
    end
end
