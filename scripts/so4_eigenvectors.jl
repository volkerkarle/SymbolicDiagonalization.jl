# SO(4) Eigenvector Algorithm - Symbolic Version
# 
# The key insight: we don't need orthonormal basis, just any basis for the 
# invariant 2-plane. We can use columns of the projection matrix directly!

using LinearAlgebra
using Symbolics

"""
Symbolic SO(4) eigenvector algorithm.

Given an SO(4) matrix A, computes all four eigenvectors symbolically.
"""
function so4_eigenpairs_symbolic(A::Matrix)
    n = size(A, 1)
    @assert n == 4 "Matrix must be 4×4"
    
    # Step 1: Compute the two rotation angles from trace invariants
    t1 = tr(A)
    t2 = tr(A * A)
    
    sum_cos = t1 / 2
    prod_cos = (t1^2 - t2 - 4) / 8
    
    discriminant = sum_cos^2 - 4 * prod_cos
    sqrt_disc = sqrt(discriminant)
    
    cos_θ1 = (sum_cos + sqrt_disc) / 2
    cos_θ2 = (sum_cos - sqrt_disc) / 2
    
    sin_θ1 = sqrt(1 - cos_θ1^2)
    sin_θ2 = sqrt(1 - cos_θ2^2)
    
    # Step 2: Compute projection operators onto invariant 2-planes
    I4 = Matrix{eltype(A)}(I, 4, 4)
    P1 = A^2 - 2*cos_θ2*A + I4  # Range = θ₁ eigenspace
    P2 = A^2 - 2*cos_θ1*A + I4  # Range = θ₂ eigenspace
    
    # Step 3: Get basis vectors from columns of P
    # P has rank 2, so we need to find 2 linearly independent columns
    # For a generic SO(4), columns 1 and 2 should work
    u1_1 = P1[:, 1]
    u1_2 = P1[:, 2]
    u2_1 = P2[:, 1]
    u2_2 = P2[:, 2]
    
    # Step 4: Compute how A acts on this basis
    # A*u1_1 = α*u1_1 + β*u1_2  (for some α, β)
    # We need to find the 2×2 matrix representation of A in this basis
    # 
    # Actually, easier approach: 
    # The eigenvalue equation A*v = λ*v with v in span{u1_1, u1_2}
    # means v = c1*u1_1 + c2*u1_2 and we need (A - λI)(c1*u1_1 + c2*u1_2) = 0
    #
    # Even simpler: just solve the linear system!
    # (A - λI) has a 2D null space. We can find it by:
    # v = adj(A - λI) * any_column, or use the formula for 2D nullspace
    
    # For eigenvalue λ, a null vector of (A - λI) is any column of adj(A - λI)
    # that is nonzero. For a rank-2 matrix (A - λI), we can also use:
    # v = (A - μI) * w for any w not in the μ-eigenspace, where μ ≠ λ
    
    λ1_plus = cos_θ1 + im*sin_θ1
    λ1_minus = cos_θ1 - im*sin_θ1
    λ2_plus = cos_θ2 + im*sin_θ2
    λ2_minus = cos_θ2 - im*sin_θ2
    
    # Use the projection trick:
    # P1 = (A - λ2_plus*I)(A - λ2_minus*I) projects onto λ1 eigenspace
    # So any column of P1 is in the λ1 eigenspace!
    # Then (A - λ1_minus*I) * P1[:,1] is in the λ1_plus eigenspace (1D)
    
    # For λ1_plus eigenspace:
    M1_plus = A - λ1_minus*I4
    v1_plus = M1_plus * u1_1  # Project P1 column onto λ1_plus eigenspace
    
    # For λ1_minus eigenspace:
    M1_minus = A - λ1_plus*I4
    v1_minus = M1_minus * u1_1
    
    # For λ2_plus eigenspace:
    M2_plus = A - λ2_minus*I4
    v2_plus = M2_plus * u2_1
    
    # For λ2_minus eigenspace:
    M2_minus = A - λ2_plus*I4
    v2_minus = M2_minus * u2_1
    
    return [
        (λ1_plus, v1_plus),
        (λ1_minus, v1_minus),
        (λ2_plus, v2_plus),
        (λ2_minus, v2_minus)
    ]
end

"""
Verify eigenpairs numerically.
"""
function verify_eigenpairs(A, pairs; atol=1e-10)
    all_passed = true
    for (λ, v) in pairs
        Av = A * v
        λv = λ * v
        passed = isapprox(Av, λv, atol=atol)
        if !passed
            println("FAILED: λ = $λ")
            println("  ||Av - λv|| = ", norm(Av - λv))
            all_passed = false
        end
    end
    return all_passed
end

# ============================================================================
# Test: Numerical verification
# ============================================================================

function test_numeric()
    println("=" ^ 70)
    println("Testing SO(4) symbolic algorithm (numerical evaluation)")
    println("=" ^ 70)
    println()
    
    n_passed = 0
    n_trials = 10
    
    for trial in 1:n_trials
        θ1 = rand() * π
        θ2 = rand() * π
        
        A_block = [cos(θ1) -sin(θ1) 0 0;
                   sin(θ1)  cos(θ1) 0 0;
                   0 0 cos(θ2) -sin(θ2);
                   0 0 sin(θ2)  cos(θ2)]
        
        Q = qr(randn(4,4)).Q |> Matrix
        A = Q * A_block * Q'
        
        pairs = so4_eigenpairs_symbolic(A)
        passed = verify_eigenpairs(A, pairs)
        
        if passed
            n_passed += 1
            println("Trial $trial: ✓ PASSED")
        else
            println("Trial $trial: ✗ FAILED")
        end
    end
    
    println()
    println("Results: $n_passed / $n_trials passed")
    return n_passed == n_trials
end

# ============================================================================
# Test: Fully symbolic
# ============================================================================

function test_symbolic()
    println()
    println("=" ^ 70)
    println("Testing with fully symbolic SO(4) matrix")
    println("=" ^ 70)
    println()
    
    @variables θ φ
    
    # Block diagonal SO(4) with symbolic angles
    A = [cos(θ) -sin(θ) 0 0;
         sin(θ)  cos(θ) 0 0;
         0 0 cos(φ) -sin(φ);
         0 0 sin(φ)  cos(φ)]
    
    println("Symbolic SO(4) matrix:")
    display(A)
    println()
    
    pairs = so4_eigenpairs_symbolic(A)
    
    println("Eigenpairs:")
    for (i, (λ, v)) in enumerate(pairs)
        println("  λ$i = ", Symbolics.simplify(λ))
        println("  v$i = ", Symbolics.simplify.(v))
        println()
    end
    
    # Verify symbolically: A*v - λ*v should simplify to zero
    println("Verification (should be zero vectors):")
    for (i, (λ, v)) in enumerate(pairs)
        residual = A * v - λ * v
        simplified = Symbolics.simplify.(residual)
        println("  A*v$i - λ$i*v$i = ", simplified)
    end
end

# Run tests
test_numeric()
test_symbolic()
