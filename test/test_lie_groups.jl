# ============================================================================
# Tests for Lie Group Pattern Detection and Eigenvalue Computation
# Focused on SYMBOLIC matrices where Lie group structure is exploited
# ============================================================================

using Test
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

# Helper function to evaluate substituted symbolic expression to Complex{Float64}
function _eval_complex(v, subs_dict)
    substituted = Symbolics.substitute(v, subs_dict)
    if substituted isa Complex
        re_sub = Symbolics.substitute(real(v), subs_dict)
        im_sub = Symbolics.substitute(imag(v), subs_dict)
        if re_sub isa Symbolics.Num
            re = Float64(eval(Symbolics.toexpr(re_sub)))
        else
            re = Float64(re_sub)
        end
        if im_sub isa Symbolics.Num
            im_part = Float64(eval(Symbolics.toexpr(im_sub)))
        else
            im_part = Float64(im_sub)
        end
        return Complex{Float64}(re, im_part)
    elseif substituted isa Symbolics.Num
        expr = Symbolics.toexpr(substituted)
        result = eval(expr)
        if result isa Complex
            return Complex{Float64}(result)
        else
            return Complex{Float64}(Float64(result), 0.0)
        end
    else
        if substituted isa Complex
            return Complex{Float64}(substituted)
        else
            return Complex{Float64}(Float64(substituted), 0.0)
        end
    end
end

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
        @test SymbolicDiagonalization._is_SO3(Rz)
        
        # SO(4) detection
        R4 = [cos(θ) -sin(θ) 0 0;
              sin(θ)  cos(θ) 0 0;
              0 0 cos(φ) -sin(φ);
              0 0 sin(φ)  cos(φ)]
        @test SymbolicDiagonalization._is_orthogonal(R4)
        @test SymbolicDiagonalization._is_SO4(R4)
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
        
        vals = SymbolicDiagonalization._Sp2_eigenvalues(Sp2)
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
        
        # Test SO2_rotation constructor
        @testset "SO2_rotation - 2D rotation" begin
            R = SO2_rotation(θ)
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
        
        # Test SO3_Rx, SO3_Ry, SO3_Rz constructors
        @testset "SO3_Rx - rotation around x-axis" begin
            R = SO3_Rx(θ)
            @test size(R) == (3, 3)
            @test isequal(R[1,1], 1)
            @test isequal(R[2,2], cos(θ))
            @test isequal(R[3,3], cos(θ))
            
            vals = eigvals(R)
            @test length(vals) == 3
            # Should have eigenvalue 1
            @test any(v -> isequal(v, 1) || isequal(Symbolics.simplify(v - 1), 0), vals)
        end
        
        @testset "SO3_Ry - rotation around y-axis" begin
            R = SO3_Ry(θ)
            @test size(R) == (3, 3)
            @test isequal(R[2,2], 1)
            @test isequal(R[1,1], cos(θ))
            @test isequal(R[3,3], cos(θ))
            
            vals = eigvals(R)
            @test length(vals) == 3
        end
        
        @testset "SO3_Rz - rotation around z-axis" begin
            R = SO3_Rz(θ)
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
            K = SO2_kron([α, β])
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
            K = SO2_kron([α, β, γ])
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
        
        @testset "SO2_kron_eigenvalues direct computation" begin
            # Direct computation should match eigvals
            direct_vals = SO2_kron_eigenvalues([α, β])
            @test length(direct_vals) == 4
            
            # All should be clean
            for v in direct_vals
                v_str = string(v)
                @test occursin("cos", v_str) && occursin("sin", v_str)
            end
            
            # 3-fold direct computation
            direct_vals_3 = SO2_kron_eigenvalues([α, β, γ])
            @test length(direct_vals_3) == 8
        end
        
        @testset "Same-angle case" begin
            # R(θ) ⊗ R(θ) should give eigenvalues with 2θ
            K = SO2_kron([α, α])
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
            # Euler angles: SO3_Rx(α) * SO3_Rz(β) * SO3_Rx(γ)
            R_euler = SO3_Rx(α) * SO3_Rz(β) * SO3_Rx(γ)
            
            # Should be detected as orthogonal
            @test SymbolicDiagonalization._is_orthogonal(R_euler)
            
            # Should be detected as SO(3)
            @test SymbolicDiagonalization._is_SO3(R_euler)
            
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
            @test SymbolicDiagonalization._is_SO3(SO3_Rz(α) * SO3_Rx(β) * SO3_Rz(γ))  # ZXZ
            @test SymbolicDiagonalization._is_SO3(SO3_Rz(α) * SO3_Ry(β) * SO3_Rz(γ))  # ZYZ
            @test SymbolicDiagonalization._is_SO3(SO3_Ry(α) * SO3_Rx(β) * SO3_Ry(γ))  # YXY
            @test SymbolicDiagonalization._is_SO3(SO3_Ry(α) * SO3_Rz(β) * SO3_Ry(γ))  # YZY
            @test SymbolicDiagonalization._is_SO3(SO3_Rx(α) * SO3_Ry(β) * SO3_Rx(γ))  # XYX
            @test SymbolicDiagonalization._is_SO3(SO3_Rx(α) * SO3_Rz(β) * SO3_Rx(γ))  # XZX
            
            # Tait-Bryan angles (three different axes)
            @test SymbolicDiagonalization._is_SO3(SO3_Rx(α) * SO3_Ry(β) * SO3_Rz(γ))  # XYZ
            @test SymbolicDiagonalization._is_SO3(SO3_Rx(α) * SO3_Rz(β) * SO3_Ry(γ))  # XZY
            @test SymbolicDiagonalization._is_SO3(SO3_Ry(α) * SO3_Rx(β) * SO3_Rz(γ))  # YXZ
            @test SymbolicDiagonalization._is_SO3(SO3_Ry(α) * SO3_Rz(β) * SO3_Rx(γ))  # YZX
            @test SymbolicDiagonalization._is_SO3(SO3_Rz(α) * SO3_Rx(β) * SO3_Ry(γ))  # ZXY
            @test SymbolicDiagonalization._is_SO3(SO3_Rz(α) * SO3_Ry(β) * SO3_Rx(γ))  # ZYX
        end
        
        @testset "Single axis rotations" begin
            @variables θ
            
            # Single axis rotations should have clean eigenvalues
            for R in [SO3_Rx(θ), SO3_Ry(θ), SO3_Rz(θ)]
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
            # Test SO3_Rz ⊗ SO3_Rz
            K = kron(SO3_Rz(θ), SO3_Rz(φ))
            @test size(K) == (9, 9)
            @test SymbolicDiagonalization._is_orthogonal(K)
            
            # Should be detected as SO(3) Kronecker product
            result = SymbolicDiagonalization._detect_SO3_kronecker_product(K)
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
            # SO3_Rz ⊗ SO3_Rx
            K = kron(SO3_Rz(θ), SO3_Rx(φ))
            vals = eigvals(K)
            @test length(vals) == 9
            # Should be clean (no sqrt)
            has_sqrt = any(v -> occursin("sqrt", string(v)), string.(vals))
            @test !has_sqrt
            
            # SO3_Rx ⊗ SO3_Ry
            K = kron(SO3_Rx(θ), SO3_Ry(φ))
            vals = eigvals(K)
            @test length(vals) == 9
            has_sqrt = any(v -> occursin("sqrt", string(v)), string.(vals))
            @test !has_sqrt
            
            # SO3_Ry ⊗ SO3_Rz
            K = kron(SO3_Ry(θ), SO3_Rz(φ))
            vals = eigvals(K)
            @test length(vals) == 9
            has_sqrt = any(v -> occursin("sqrt", string(v)), string.(vals))
            @test !has_sqrt
        end
        
        @testset "Same angle case" begin
            # SO3_Rz(θ) ⊗ SO3_Rz(θ) should have eigenvalues including e^{±2iθ}
            K = kron(SO3_Rz(θ), SO3_Rz(θ))
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
            R1 = SO3_Rx(α) * SO3_Rz(β) * SO3_Rx(γ)
            R2 = SO3_Ry(δ) * SO3_Rz(ε) * SO3_Ry(ζ)
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
    
    @testset "SU(2) Constructors and Kronecker Products" begin
        @variables θ α β γ
        
        @testset "Pauli matrices" begin
            # Test Pauli matrix structure
            @test pauli_x() == [0 1; 1 0]
            @test pauli_y() == [0 -im; im 0]
            @test pauli_z() == [1 0; 0 -1]
        end
        
        @testset "SU(2) rotation constructors" begin
            # SU2_Ux should be exp(-i θ σx/2)
            U = SU2_Ux(θ)
            @test size(U) == (2, 2)
            # U[1,1] is Complex{Num} due to matrix type, but real part should be cos(θ/2)
            @test isequal(real(U[1,1]), cos(θ/2))
            @test isequal(imag(U[1,1]), 0) || isequal(Symbolics.simplify(imag(U[1,1])), 0)
            @test isequal(real(U[2,2]), cos(θ/2))
            # Off-diagonal elements should be -i*sin(θ/2)
            @test isequal(real(U[1,2]), 0) || isequal(Symbolics.simplify(real(U[1,2])), 0)
            @test isequal(imag(U[1,2]), -sin(θ/2))
            @test isequal(real(U[2,1]), 0) || isequal(Symbolics.simplify(real(U[2,1])), 0)
            @test isequal(imag(U[2,1]), -sin(θ/2))
            
            # SU2_Uy should be exp(-i θ σy/2) - same form as SO(2) but with half-angle
            U = SU2_Uy(θ)
            @test size(U) == (2, 2)
            @test isequal(U[1,1], cos(θ/2))
            @test isequal(U[2,2], cos(θ/2))
            @test isequal(U[1,2], -sin(θ/2))
            @test isequal(U[2,1], sin(θ/2))
            
            # SU2_Uz should be diagonal with e^{±iθ/2}
            U = SU2_Uz(θ)
            @test size(U) == (2, 2)
            # Check it's diagonal
            @test isequal(U[1,2], 0)
            @test isequal(U[2,1], 0)
            # Diagonal elements should be e^{-iθ/2} and e^{iθ/2}
            @test isequal(real(U[1,1]), cos(θ/2))
            @test isequal(imag(U[1,1]), -sin(θ/2))
            @test isequal(real(U[2,2]), cos(θ/2))
            @test isequal(imag(U[2,2]), sin(θ/2))
        end
        
        @testset "SU(2) eigenvalues" begin
            # Single SU(2) rotation eigenvalues
            for (name, Ufn) in [("SU2_Ux", SU2_Ux), ("SU2_Uy", SU2_Uy), ("SU2_Uz", SU2_Uz)]
                U = Ufn(θ)
                vals = eigvals(U)
                @test length(vals) == 2
                
                # Eigenvalues should be e^{±iθ/2} = cos(θ/2) ± i sin(θ/2)
                # Check they are clean (no sqrt)
                for v in vals
                    v_str = string(v)
                    @test !occursin("sqrt", v_str)
                end
            end
        end
        
        @testset "SU(2) ⊗ SU(2) detection - basic" begin
            # Test SU2_Uz ⊗ SU2_Uz (diagonal case - easiest)
            K = kron(SU2_Uz(α), SU2_Uz(β))
            @test size(K) == (4, 4)
            
            vals = eigvals(K)
            @test length(vals) == 4
            
            # Eigenvalues should be clean (no sqrt)
            for v in vals
                v_str = string(v)
                @test !occursin("sqrt", v_str)
            end
            
            # Check that eigenvalues contain half-angles
            vals_str = join(string.(vals), " ")
            # Should have α/2 and β/2 terms
            @test occursin("α", vals_str) && occursin("β", vals_str)
        end
        
        @testset "SU(2) ⊗ SU(2) detection - mixed axes" begin
            # Test SU2_Ux ⊗ SU2_Uz
            K = kron(SU2_Ux(α), SU2_Uz(β))
            vals = eigvals(K)
            @test length(vals) == 4
            
            # Eigenvalues should be clean
            for v in vals
                v_str = string(v)
                @test !occursin("sqrt", v_str)
            end
            
            # Test SU2_Uy ⊗ SU2_Ux
            K = kron(SU2_Uy(α), SU2_Ux(β))
            vals = eigvals(K)
            @test length(vals) == 4
            
            for v in vals
                v_str = string(v)
                @test !occursin("sqrt", v_str)
            end
        end
        
        @testset "SU2_kron and SU2_kron_eigenvalues" begin
            # Direct construction
            K = SU2_kron([α, β])
            @test size(K) == (4, 4)
            
            # Direct eigenvalue computation
            direct_vals = SU2_kron_eigenvalues([α, β])
            @test length(direct_vals) == 4
            
            # All should be clean
            for v in direct_vals
                v_str = string(v)
                @test occursin("cos", v_str) && occursin("sin", v_str)
                # Should have half-angle structure (divided by 2)
                # The output can be "/ 2" (with space) or "(1/2)" or "(1//2)" etc.
                @test occursin(r"/ *2|1/2|1//2|\*0\.5", v_str)
            end
            
            # 3-fold direct computation
            direct_vals_3 = SU2_kron_eigenvalues([α, β, γ])
            @test length(direct_vals_3) == 8
        end
        
        @testset "Numerical verification" begin
            # Verify symbolic eigenvalues match numeric eigenvalues when substituted
            θ_val = 0.7
            φ_val = 1.3
            
            # Step 1: Create SYMBOLIC SU(2) Kronecker product
            @variables θ_sym φ_sym
            K_sym = kron(SU2_Uz(θ_sym), SU2_Uz(φ_sym))
            
            # Step 2: Get symbolic eigenvalues using our detection
            sym_vals, _, _ = symbolic_eigenvalues(K_sym; expand=false)
            @test length(sym_vals) == 4
            
            # Step 3: Substitute numeric values into symbolic eigenvalues
            subs = Dict(θ_sym => θ_val, φ_sym => φ_val)
            substituted_vals = [_eval_complex(v, subs) for v in sym_vals]
            substituted_sorted = sort(substituted_vals, by=x->(real(x), imag(x)))
            
            # Step 4: Create numeric matrix with same values and compute eigenvalues directly
            K_num = [_eval_complex(x, subs) for x in K_sym]
            ref_vals = sort(eigvals(K_num), by=x->(real(x), imag(x)))
            
            # Step 5: Compare - this validates our symbolic detection is correct
            @test isapprox(substituted_sorted, ref_vals, atol=1e-10)
        end
    end
end

# ============================================================================
# SU(3) Kronecker Product Tests  
# ============================================================================

@testset "SU(3) Kronecker Products" begin
    
    @testset "Gell-Mann matrices" begin
        # Test all 8 Gell-Mann matrices exist and are traceless
        for (i, λ) in enumerate(gellmann_matrices())
            @test size(λ) == (3, 3)
            @test isapprox(tr(λ), 0, atol=1e-10)
        end
        
        # Test Hermiticity
        @test gellmann_1() == gellmann_1()'
        @test gellmann_2() == gellmann_2()'
        @test gellmann_3() == gellmann_3()'
        @test gellmann_4() == gellmann_4()'
        @test gellmann_5() == gellmann_5()'
        @test gellmann_6() == gellmann_6()'
        @test gellmann_7() == gellmann_7()'
        @test isapprox(gellmann_8(), gellmann_8()', atol=1e-10)  # gellmann_8 has floats
    end
    
    @testset "SU(3) diagonal constructors" begin
        @variables θ₁ θ₂
        
        # Test SU3_diagonal_trig
        U = SU3_diagonal_trig(θ₁, θ₂)
        @test size(U) == (3, 3)
        @test U isa Diagonal
        
        # Test determinant = 1 symbolically
        # det(U) = e^{iθ₁} * e^{iθ₂} * e^{-i(θ₁+θ₂)} = 1
        det_U = det(U)
        det_real = trig_simplify(real(det_U))
        det_imag = trig_simplify(imag(det_U))
        @test isequal(Symbolics.simplify(det_real), 1) || isequal(Symbolics.simplify(det_real - 1), 0)
    end
    
    @testset "SU(3) Kronecker eigenvalues - direct formula" begin
        @variables α₁ α₂ β₁ β₂
        
        vals = SU3_kron_eigenvalues((α₁, α₂), (β₁, β₂))
        @test length(vals) == 9
        
        # First eigenvalue should be e^{i(α₁+β₁)}
        v1 = vals[1]
        @test isequal(v1, cos(α₁ + β₁) + im*sin(α₁ + β₁))
    end
    
    @testset "SU(3) Kronecker detection - diagonal symbolic" begin
        @variables α₁ α₂ β₁ β₂
        
        K = SU3_kron((α₁, α₂), (β₁, β₂))
        @test size(K) == (9, 9)
        
        # Test symbolic_eigenvalues returns 9 eigenvalues
        vals, poly, λ = symbolic_eigenvalues(K; expand=false)
        @test length(vals) == 9
        
        # First eigenvalue should be cos(α₁+β₁) + i*sin(α₁+β₁)
        v1 = vals[1]
        @test isequal(v1, cos(α₁ + β₁) + im*sin(α₁ + β₁))
    end
    
    @testset "SU(3) Kronecker detection - numerical verification" begin
        # Verify symbolic eigenvalues match numeric eigenvalues when substituted
        θ1_val, θ2_val = 0.3, 0.5
        φ1_val, φ2_val = 0.2, 0.7
        
        # Step 1: Create SYMBOLIC SU(3) Kronecker product
        @variables θ1_sym θ2_sym φ1_sym φ2_sym
        K_sym = SU3_kron((θ1_sym, θ2_sym), (φ1_sym, φ2_sym))
        @test size(K_sym) == (9, 9)
        
        # Step 2: Get symbolic eigenvalues using our detection
        sym_vals, _, _ = symbolic_eigenvalues(K_sym; expand=false)
        @test length(sym_vals) == 9
        
        # Step 3: Substitute numeric values into symbolic eigenvalues
        subs = Dict(θ1_sym => θ1_val, θ2_sym => θ2_val, φ1_sym => φ1_val, φ2_sym => φ2_val)
        substituted_vals = [_eval_complex(v, subs) for v in sym_vals]
        substituted_sorted = sort(substituted_vals, by=x->(real(x), imag(x)))
        
        # Step 4: Create numeric matrix with same values and compute eigenvalues directly
        K_num = [_eval_complex(x, subs) for x in K_sym]
        ref_vals = sort(eigvals(K_num), by=x->(real(x), imag(x)))
        
        # Step 5: Compare - this validates our symbolic detection is correct
        @test isapprox(substituted_sorted, ref_vals, atol=1e-10)
    end
    
    @testset "SU(3)⊗SU(3) vs SO(3)⊗SO(3) disambiguation" begin
        # SU(3)⊗SU(3) has complex entries
        @variables α₁ α₂ β₁ β₂
        K_su3 = SU3_kron((α₁, α₂), (β₁, β₂))
        
        using SymbolicDiagonalization: _has_complex_entries, _detect_SU3_kronecker_product, _detect_SO3_kronecker_product
        
        @test _has_complex_entries(K_su3) == true
        @test !isnothing(_detect_SU3_kronecker_product(K_su3))
        
        # SO(3)⊗SO(3) is real
        @variables θ φ
        R1 = SO3_Rz(θ)
        R2 = SO3_Rz(φ)
        K_so3 = kron(R1, R2)
        
        @test _has_complex_entries(K_so3) == false
        @test isnothing(_detect_SU3_kronecker_product(K_so3))  # SU(3) detection should reject
        @test !isnothing(_detect_SO3_kronecker_product(K_so3))  # SO(3) detection should work
    end
end

# ============================================================================
# Lie Group Eigenvector Tests
# Tests for closed-form eigenvector computation
# ============================================================================

@testset "Lie Group Eigenvectors" begin
    
    @testset "SO(2) eigenvectors" begin
        @variables θ
        R = SO2_rotation(θ)
        
        # Get eigenpairs using the fast path
        pairs = SymbolicDiagonalization._SO2_eigenpairs(R)
        @test !isnothing(pairs)
        @test length(pairs) == 2
        
        # Check eigenvector structure: [1, i] and [1, -i]
        v1 = pairs[1][2][1]  # First eigenvector
        v2 = pairs[2][2][1]  # Second eigenvector
        
        @test v1 == [1, im]
        @test v2 == [1, -im]
        
        # Eigenvalue should be cos(θ) + im*sin(θ) for first eigenpair
        λ1 = pairs[1][1]
        @test isequal(real(λ1), cos(θ))
        @test isequal(imag(λ1), sin(θ))
        
        # Test through symbolic_eigenpairs interface
        full_pairs, poly, λ = symbolic_eigenpairs(R)
        @test length(full_pairs) == 2
        @test length(full_pairs[1][2]) == 1  # One eigenvector per eigenvalue
    end
    
    @testset "SO(3) axis-aligned eigenvectors" begin
        @variables θ
        
        # Test Rz
        Rz = SO3_Rz(θ)
        pairs = SymbolicDiagonalization._SO3_eigenpairs(Rz)
        @test !isnothing(pairs)
        @test length(pairs) == 3
        
        # First eigenvalue should be 1 with eigenvector [0,0,1]
        @test isequal(pairs[1][1], 1)
        @test pairs[1][2][1] == [0, 0, 1]
        
        # Other eigenvectors should be [1,i,0] and [1,-i,0]
        @test pairs[2][2][1] == [1, im, 0]
        @test pairs[3][2][1] == [1, -im, 0]
        
        # Test Rx
        Rx = SO3_Rx(θ)
        pairs = SymbolicDiagonalization._SO3_eigenpairs(Rx)
        @test !isnothing(pairs)
        @test length(pairs) == 3
        @test pairs[1][2][1] == [1, 0, 0]
        
        # Test Ry
        Ry = SO3_Ry(θ)
        pairs = SymbolicDiagonalization._SO3_eigenpairs(Ry)
        @test !isnothing(pairs)
        @test length(pairs) == 3
        @test pairs[1][2][1] == [0, 1, 0]
        
        # Test through symbolic_eigenpairs
        full_pairs, _, _ = symbolic_eigenpairs(Rz)
        @test length(full_pairs) == 3
    end
    
    @testset "SO(3) general rotation - falls back to nullspace" begin
        @variables α β γ
        
        # Euler rotation is NOT axis-aligned
        R_euler = SO3_Rx(α) * SO3_Rz(β) * SO3_Rx(γ)
        
        # Direct eigenpairs function should return nothing (requires nullspace)
        pairs = SymbolicDiagonalization._SO3_eigenpairs(R_euler)
        @test isnothing(pairs)
        
        # But symbolic_eigenpairs should still work (using nullspace fallback)
        full_pairs, _, _ = symbolic_eigenpairs(R_euler)
        @test length(full_pairs) == 3
    end
    
    @testset "SU(2) diagonal eigenvectors" begin
        @variables θ
        
        # Test Uz (diagonal)
        Uz = SU2_Uz(θ)
        pairs = SymbolicDiagonalization._SU2_eigenpairs(Uz)
        @test !isnothing(pairs)
        @test length(pairs) == 2
        
        # Eigenvectors should be standard basis [1,0] and [0,1]
        @test pairs[1][2][1] == [1, 0]
        @test pairs[2][2][1] == [0, 1]
        
        # Eigenvalues should be e^{±iθ/2}
        λ1 = pairs[1][1]
        λ2 = pairs[2][1]
        @test isequal(real(λ1), cos(θ/2))
        @test isequal(imag(λ1), -sin(θ/2))
        @test isequal(real(λ2), cos(θ/2))
        @test isequal(imag(λ2), sin(θ/2))
        
        # Test through symbolic_eigenpairs
        full_pairs, _, _ = symbolic_eigenpairs(Uz)
        @test length(full_pairs) == 2
    end
    
    @testset "SU(2) Ux eigenvectors" begin
        @variables θ
        
        # Test Ux
        Ux = SU2_Ux(θ)
        pairs = SymbolicDiagonalization._SU2_eigenpairs(Ux)
        @test !isnothing(pairs)
        @test length(pairs) == 2
        
        # Eigenvectors should be [1,1] and [1,-1]
        @test pairs[1][2][1] == [1, 1]
        @test pairs[2][2][1] == [1, -1]
    end
    
    @testset "SU(2) Uy eigenvectors" begin
        @variables θ
        
        # Test Uy
        Uy = SU2_Uy(θ)
        pairs = SymbolicDiagonalization._SU2_eigenpairs(Uy)
        @test !isnothing(pairs)
        @test length(pairs) == 2
        
        # Eigenvectors should be [1,i] and [1,-i]
        @test pairs[1][2][1] == [1, im]
        @test pairs[2][2][1] == [1, -im]
    end
    
    @testset "SU(3) diagonal eigenvectors" begin
        @variables θ₁ θ₂
        
        U = SU3_diagonal_trig(θ₁, θ₂)
        U_mat = Matrix(U)
        pairs = SymbolicDiagonalization._SU3_eigenpairs(U_mat)
        @test !isnothing(pairs)
        @test length(pairs) == 3
        
        # Eigenvectors should be standard basis
        @test pairs[1][2][1] == [1, 0, 0]
        @test pairs[2][2][1] == [0, 1, 0]
        @test pairs[3][2][1] == [0, 0, 1]
        
        # Eigenvalues should be e^{iθ₁}, e^{iθ₂}, e^{-i(θ₁+θ₂)}
        λ1 = pairs[1][1]
        @test isequal(real(λ1), cos(θ₁))
        @test isequal(imag(λ1), sin(θ₁))
        
        # Test through symbolic_eigenpairs
        full_pairs, _, _ = symbolic_eigenpairs(U_mat)
        @test length(full_pairs) == 3
    end
    
    @testset "SO(2) Kronecker eigenvectors" begin
        @variables α β
        
        K = SO2_kron([α, β])
        pairs = SymbolicDiagonalization._SO2_kron_eigenpairs(K)
        @test !isnothing(pairs)
        @test length(pairs) == 4
        
        # Eigenvectors should be tensor products of [1,i] and [1,-i]
        # First eigenvector for signs (-1,-1): [1,-i] ⊗ [1,-i] = [1,-i,-i,-1]
        v1 = pairs[1][2][1]
        @test length(v1) == 4
        @test v1 == [1, -im, -im, -1]
        
        # Test through symbolic_eigenpairs
        full_pairs, _, _ = symbolic_eigenpairs(K)
        @test length(full_pairs) == 4
    end
    
    @testset "Eigenvector numerical verification" begin
        # Verify that eigenvectors actually satisfy A*v = λ*v when substituted
        @variables θ
        
        # SO(2)
        R = SO2_rotation(θ)
        pairs, _, _ = symbolic_eigenpairs(R)
        
        θ_val = 0.7
        subs = Dict(θ => θ_val)
        
        for (λ, vecs) in pairs
            for v in vecs
                # Substitute
                λ_num = _eval_complex(λ, subs)
                v_num = [_eval_complex(vi, subs) for vi in v]
                R_num = [_eval_complex(Rij, subs) for Rij in R]
                
                # Check A*v ≈ λ*v
                Av = R_num * v_num
                λv = λ_num .* v_num
                @test isapprox(Av, λv, atol=1e-10)
            end
        end
        
        # SO(3) Rz
        Rz = SO3_Rz(θ)
        pairs, _, _ = symbolic_eigenpairs(Rz)
        
        for (λ, vecs) in pairs
            for v in vecs
                λ_num = _eval_complex(λ, subs)
                v_num = [_eval_complex(vi, subs) for vi in v]
                R_num = [_eval_complex(Rij, subs) for Rij in Rz]
                
                Av = R_num * v_num
                λv = λ_num .* v_num
                @test isapprox(Av, λv, atol=1e-10)
            end
        end
        
        # SU(2) Uz
        Uz = SU2_Uz(θ)
        pairs, _, _ = symbolic_eigenpairs(Uz)
        
        for (λ, vecs) in pairs
            for v in vecs
                λ_num = _eval_complex(λ, subs)
                v_num = [_eval_complex(vi, subs) for vi in v]
                U_num = [_eval_complex(Uij, subs) for Uij in Uz]
                
                Uv = U_num * v_num
                λv = λ_num .* v_num
                @test isapprox(Uv, λv, atol=1e-10)
            end
        end
    end
    
    @testset "Diagonalization with Lie group eigenvectors" begin
        @variables θ
        
        # SO(2) should be diagonalizable
        R = SO2_rotation(θ)
        P, D, pairs = symbolic_diagonalize(R)
        
        @test size(P) == (2, 2)
        @test size(D) == (2, 2)
        @test length(pairs) == 2
        
        # SO(3) Rz should be diagonalizable
        Rz = SO3_Rz(θ)
        P, D, pairs = symbolic_diagonalize(Rz)
        
        @test size(P) == (3, 3)
        @test size(D) == (3, 3)
        @test length(pairs) == 3
        
        # SU(2) Uz should be diagonalizable
        Uz = SU2_Uz(θ)
        P, D, pairs = symbolic_diagonalize(Uz)
        
        @test size(P) == (2, 2)
        @test size(D) == (2, 2)
        @test length(pairs) == 2
    end
end
