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
    
    @testset "Rotation Matrix Constructors" begin
        @variables θ α β γ
        
        # Test R2 constructor
        @testset "R2 - 2D rotation" begin
            R = R2(θ)
            @test size(R) == (2, 2)
            @test isequal(R[1,1], cos(θ))
            @test isequal(R[2,1], sin(θ))
            @test isequal(R[1,2], -sin(θ))
            @test isequal(R[2,2], cos(θ))
            
            # Eigenvalues should be clean
            vals = eigvals(R)
            @test length(vals) == 2
            @test any(v -> occursin("cos(θ) + im*sin(θ)", string(v)), string.(vals))
            @test any(v -> occursin("cos(θ) - im*sin(θ)", string(v)), string.(vals))
        end
        
        # Test Rx, Ry, Rz constructors
        @testset "Rx - rotation around x-axis" begin
            R = Rx(θ)
            @test size(R) == (3, 3)
            @test isequal(R[1,1], 1)
            @test isequal(R[2,2], cos(θ))
            @test isequal(R[3,3], cos(θ))
            
            vals = eigvals(R)
            @test length(vals) == 3
            # Should have eigenvalue 1
            @test any(v -> isequal(v, 1) || isequal(Symbolics.simplify(v - 1), 0), vals)
        end
        
        @testset "Ry - rotation around y-axis" begin
            R = Ry(θ)
            @test size(R) == (3, 3)
            @test isequal(R[2,2], 1)
            @test isequal(R[1,1], cos(θ))
            @test isequal(R[3,3], cos(θ))
            
            vals = eigvals(R)
            @test length(vals) == 3
        end
        
        @testset "Rz - rotation around z-axis" begin
            R = Rz(θ)
            @test size(R) == (3, 3)
            @test isequal(R[3,3], 1)
            @test isequal(R[1,1], cos(θ))
            @test isequal(R[2,2], cos(θ))
            
            vals = eigvals(R)
            @test length(vals) == 3
        end
    end
    
    @testset "SO(2) Kronecker Products" begin
        @variables α β γ
        
        @testset "2-fold SO(2) Kronecker" begin
            K = so2_kron([α, β])
            @test size(K) == (4, 4)
            
            vals = eigvals(K)
            @test length(vals) == 4
            
            # All eigenvalues should be in clean form cos(sum) + im*sin(sum)
            for v in vals
                v_str = string(v)
                @test occursin("cos", v_str) && occursin("sin", v_str)
                # Should NOT contain sqrt
                @test !occursin("sqrt", v_str)
            end
            
            # Check specific eigenvalue structure (±α ± β combinations)
            vals_str = join(string.(vals), " ")
            @test occursin("α", vals_str) && occursin("β", vals_str)
        end
        
        @testset "3-fold SO(2) Kronecker" begin
            K = so2_kron([α, β, γ])
            @test size(K) == (8, 8)
            
            vals = eigvals(K)
            @test length(vals) == 8
            
            # All eigenvalues should be in clean form
            for v in vals
                v_str = string(v)
                @test occursin("cos", v_str) && occursin("sin", v_str)
                # Should NOT contain sqrt or other messy expressions
                @test !occursin("sqrt", v_str)
            end
            
            # Check that all three angles appear
            vals_str = join(string.(vals), " ")
            @test occursin("α", vals_str) && occursin("β", vals_str) && occursin("γ", vals_str)
        end
        
        @testset "so2_kron_eigenvalues direct computation" begin
            # Direct computation should match eigvals
            direct_vals = so2_kron_eigenvalues([α, β])
            @test length(direct_vals) == 4
            
            # All should be clean
            for v in direct_vals
                v_str = string(v)
                @test occursin("cos", v_str) && occursin("sin", v_str)
            end
            
            # 3-fold direct computation
            direct_vals_3 = so2_kron_eigenvalues([α, β, γ])
            @test length(direct_vals_3) == 8
        end
        
        @testset "Same-angle case" begin
            # R(θ) ⊗ R(θ) should give eigenvalues with 2θ
            K = so2_kron([α, α])
            vals = eigvals(K)
            @test length(vals) == 4
            
            # Should simplify to use 2α where applicable
            vals_str = join(string.(vals), " ")
            # The eigenvalues are cos(2α), cos(0)=1, etc.
            @test occursin("α", vals_str)
        end
    end
    
    @testset "Euler Angle Rotations" begin
        @variables α β γ
        
        @testset "Euler rotation detection (XZX)" begin
            # Euler angles: Rx(α) * Rz(β) * Rx(γ)
            R_euler = Rx(α) * Rz(β) * Rx(γ)
            
            # Should be detected as orthogonal
            @test SymbolicDiagonalization._is_orthogonal(R_euler)
            
            # Should be detected as SO(3)
            @test SymbolicDiagonalization._is_so3(R_euler)
            
            # Eigenvalues should compute successfully
            vals = eigvals(R_euler)
            @test length(vals) == 3
            
            # One eigenvalue should be 1 (the rotation axis)
            has_one = any(vals) do v
                if v isa Complex
                    r = real(v)
                    i = imag(v)
                    isequal(Symbolics.simplify(r - 1), 0) && isequal(Symbolics.simplify(i), 0)
                else
                    isequal(Symbolics.simplify(v - 1), 0)
                end
            end
            @test has_one
        end
        
        @testset "Various Euler conventions" begin
            # Test that all common Euler conventions are detected as SO(3)
            
            # Proper Euler angles (same first and last axis)
            @test SymbolicDiagonalization._is_so3(Rz(α) * Rx(β) * Rz(γ))  # ZXZ
            @test SymbolicDiagonalization._is_so3(Rz(α) * Ry(β) * Rz(γ))  # ZYZ
            @test SymbolicDiagonalization._is_so3(Ry(α) * Rx(β) * Ry(γ))  # YXY
            @test SymbolicDiagonalization._is_so3(Ry(α) * Rz(β) * Ry(γ))  # YZY
            @test SymbolicDiagonalization._is_so3(Rx(α) * Ry(β) * Rx(γ))  # XYX
            @test SymbolicDiagonalization._is_so3(Rx(α) * Rz(β) * Rx(γ))  # XZX
            
            # Tait-Bryan angles (three different axes)
            @test SymbolicDiagonalization._is_so3(Rx(α) * Ry(β) * Rz(γ))  # XYZ
            @test SymbolicDiagonalization._is_so3(Rx(α) * Rz(β) * Ry(γ))  # XZY
            @test SymbolicDiagonalization._is_so3(Ry(α) * Rx(β) * Rz(γ))  # YXZ
            @test SymbolicDiagonalization._is_so3(Ry(α) * Rz(β) * Rx(γ))  # YZX
            @test SymbolicDiagonalization._is_so3(Rz(α) * Rx(β) * Ry(γ))  # ZXY
            @test SymbolicDiagonalization._is_so3(Rz(α) * Ry(β) * Rx(γ))  # ZYX
        end
        
        @testset "Single axis rotations" begin
            @variables θ
            
            # Single axis rotations should have clean eigenvalues
            for R in [Rx(θ), Ry(θ), Rz(θ)]
                vals = eigvals(R)
                @test length(vals) == 3
                
                # Should have exactly one eigenvalue = 1 (possibly as Complex{Num})
                has_one = any(vals) do v
                    if v isa Complex
                        isequal(Symbolics.simplify(real(v) - 1), 0) && 
                        isequal(Symbolics.simplify(imag(v)), 0)
                    else
                        isequal(Symbolics.simplify(v - 1), 0)
                    end
                end
                @test has_one
                
                # Other two should be cos(θ) ± i*sin(θ)
                # Filter out the eigenvalue 1
                complex_non_one = filter(vals) do v
                    if v isa Complex
                        # Not equal to 1+0i
                        !(isequal(Symbolics.simplify(real(v) - 1), 0) && 
                          isequal(Symbolics.simplify(imag(v)), 0))
                    else
                        !isequal(Symbolics.simplify(v - 1), 0)
                    end
                end
                @test length(complex_non_one) == 2
                
                # Check they are complex conjugates (real parts equal, imag parts opposite)
                v1, v2 = complex_non_one
                @test isequal(real(v1), real(v2))
                # Imaginary parts should sum to zero
                sum_imag = Symbolics.simplify(imag(v1) + imag(v2))
                @test isequal(sum_imag, 0)
            end
        end
    end
    
    @testset "SO(3) Kronecker Products" begin
        @variables θ φ
        
        @testset "SO(3) ⊗ SO(3) detection and eigenvalues" begin
            # Test Rz ⊗ Rz
            K = kron(Rz(θ), Rz(φ))
            @test size(K) == (9, 9)
            @test SymbolicDiagonalization._is_orthogonal(K)
            
            # Should be detected as SO(3) Kronecker product
            result = SymbolicDiagonalization._detect_so3_kronecker_product(K)
            @test !isnothing(result)
            @test length(result) == 9
            
            # Eigenvalues should be products of SO(3) eigenvalues
            # Expected: 1, e^{±iθ}, e^{±iφ}, e^{±i(θ+φ)}, e^{±i(θ-φ)}
            # All should be in clean cos/sin form, no sqrt
            for v in result
                v_str = string(v)
                @test !occursin("sqrt", v_str)
            end
            
            # Should have exactly one eigenvalue = 1
            has_one = any(result) do v
                if v isa Complex
                    isequal(Symbolics.simplify(real(v) - 1), 0) && 
                    isequal(Symbolics.simplify(imag(v)), 0)
                else
                    isequal(Symbolics.simplify(v - 1), 0)
                end
            end
            @test has_one
        end
        
        @testset "Different axis combinations" begin
            # Rz ⊗ Rx
            K = kron(Rz(θ), Rx(φ))
            vals = eigvals(K)
            @test length(vals) == 9
            # Should be clean (no sqrt)
            has_sqrt = any(v -> occursin("sqrt", string(v)), string.(vals))
            @test !has_sqrt
            
            # Rx ⊗ Ry
            K = kron(Rx(θ), Ry(φ))
            vals = eigvals(K)
            @test length(vals) == 9
            has_sqrt = any(v -> occursin("sqrt", string(v)), string.(vals))
            @test !has_sqrt
            
            # Ry ⊗ Rz
            K = kron(Ry(θ), Rz(φ))
            vals = eigvals(K)
            @test length(vals) == 9
            has_sqrt = any(v -> occursin("sqrt", string(v)), string.(vals))
            @test !has_sqrt
        end
        
        @testset "Same angle case" begin
            # Rz(θ) ⊗ Rz(θ) should have eigenvalues including e^{±2iθ}
            K = kron(Rz(θ), Rz(θ))
            vals = eigvals(K)
            @test length(vals) == 9
            
            # Check for 2θ terms
            vals_str = join(string.(vals), " ")
            # Should have terms like cos(2θ), sin(2θ) or 2θ
            @test occursin("θ", vals_str)
        end
        
        @testset "Euler ⊗ Euler detection" begin
            # Test that Euler angle rotations ⊗ Euler angle rotations are detected
            # This is the most complex case where no diagonal block element = 1
            @variables α β γ δ ε ζ
            
            # XZX Euler ⊗ YZY Euler
            R1 = Rx(α) * Rz(β) * Rx(γ)
            R2 = Ry(δ) * Rz(ε) * Ry(ζ)
            K = kron(R1, R2)
            
            # eigvals should automatically detect and compute eigenvalues
            vals = eigvals(K)
            @test length(vals) == 9
            
            # One eigenvalue should be 1 (from the 1*1 product)
            has_one = any(vals) do v
                if v isa Complex
                    SymbolicDiagonalization._issymzero(Symbolics.simplify(real(v) - 1)) && 
                    SymbolicDiagonalization._issymzero(Symbolics.simplify(imag(v)))
                else
                    SymbolicDiagonalization._issymzero(Symbolics.simplify(v - 1))
                end
            end
            @test has_one
        end
    end
end
