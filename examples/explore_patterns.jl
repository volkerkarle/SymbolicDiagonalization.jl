"""
Exploration script to discover closed-form eigenvalue patterns in tridiagonal matrices.

Strategy:
1. Test simple symmetric tridiagonal patterns with constant diagonals
2. Try perturbations of constant patterns
3. Look for eigenvalues that can be expressed in closed form
4. Analyze what makes certain patterns have simple eigenvalues
"""

using LinearAlgebra

"""
Create a symmetric tridiagonal matrix with given diagonal and off-diagonal values.
"""
function make_symmetric_tridiagonal(diag_val, off_diag_vals)
    n = length(off_diag_vals) + 1
    mat = zeros(Float64, n, n)
    
    # Set diagonal
    for i in 1:n
        mat[i, i] = diag_val
    end
    
    # Set off-diagonals (symmetric)
    for i in 1:length(off_diag_vals)
        mat[i, i+1] = off_diag_vals[i]
        mat[i+1, i] = off_diag_vals[i]
    end
    
    return mat
end

"""
Try to identify simple patterns in eigenvalues.
"""
function analyze_eigenvalues(eigs, params)
    println("  Eigenvalues: ", [round(e, digits=6) for e in eigs])
    
    # Check for simple patterns
    a_val = params[:a]
    b_val = get(params, :b, nothing)
    d_val = get(params, :d, nothing)
    c_val = get(params, :c, nothing)
    
    tol = 1e-10
    
    # Check if any eigenvalues are exactly a
    for (i, λ) in enumerate(eigs)
        if abs(λ - a_val) < tol
            println("    λ[$i] = a")
        end
    end
    
    # Check for a ± b patterns
    if !isnothing(b_val)
        for (i, λ) in enumerate(eigs)
            if abs(λ - (a_val + b_val)) < tol
                println("    λ[$i] = a + b")
            elseif abs(λ - (a_val - b_val)) < tol
                println("    λ[$i] = a - b")
            elseif abs(λ - (a_val + 2*b_val)) < tol
                println("    λ[$i] = a + 2b")
            elseif abs(λ - (a_val - 2*b_val)) < tol
                println("    λ[$i] = a - 2b")
            end
        end
    end
    
    # Check for square root patterns
    if !isnothing(b_val) && !isnothing(d_val)
        # Try various combinations
        combinations = [
            ("2b² + d²", 2*b_val^2 + d_val^2),
            ("b² + d²", b_val^2 + d_val^2),
            ("b² + 2d²", b_val^2 + 2*d_val^2),
            ("3b² + d²", 3*b_val^2 + d_val^2),
            ("b² + d² + 2bd", b_val^2 + d_val^2 + 2*b_val*d_val),
            ("2(b² + d²)", 2*(b_val^2 + d_val^2)),
            ("b² + 3d²", b_val^2 + 3*d_val^2),
            ("4b² + d²", 4*b_val^2 + d_val^2),
        ]
        
        for (i, λ) in enumerate(eigs)
            diff_from_a = λ - a_val
            for (name, val) in combinations
                if abs(diff_from_a - sqrt(val)) < tol
                    println("    λ[$i] = a + √($name)")
                elseif abs(diff_from_a + sqrt(val)) < tol
                    println("    λ[$i] = a - √($name)")
                end
            end
        end
    end
    
    #  Check for three-parameter patterns
    if !isnothing(b_val) && !isnothing(d_val) && !isnothing(c_val)
        combinations = [
            ("b² + d² + c²", b_val^2 + d_val^2 + c_val^2),
            ("2b² + d² + c²", 2*b_val^2 + d_val^2 + c_val^2),
        ]
        
        for (i, λ) in enumerate(eigs)
            diff_from_a = λ - a_val
            for (name, val) in combinations
                if abs(diff_from_a - sqrt(val)) < tol
                    println("    λ[$i] = a + √($name)")
                elseif abs(diff_from_a + sqrt(val)) < tol
                    println("    λ[$i] = a - √($name)")
                end
            end
        end
    end
end

println("="^70)
println("EXPLORING TRIDIAGONAL MATRIX PATTERNS")
println("="^70)

# Test parameters
test_params = Dict(:a => 1.0, :b => 2.0, :d => 3.0, :c => 4.0)

println("\n" * "="^70)
println("1. CONSTANT TRIDIAGONAL (Toeplitz)")
println("="^70)
println("Pattern: All diagonals = a, all off-diagonals = b")
println("Known: These have closed-form eigenvalues involving trigonometric functions")

for n in 3:7
    println("\n[Size $n×$n]")
    mat = make_symmetric_tridiagonal(test_params[:a], fill(test_params[:b], n-1))
    eigs = sort(real.(eigvals(mat)))
    analyze_eigenvalues(eigs, test_params)
end

println("\n" * "="^70)
println("2. SINGLE PERTURBATION PATTERNS")
println("="^70)
println("Pattern: Mostly constant off-diagonals [b,b,...,b] with one different value")

println("\n[5×5: off-diagonal = [d, b, b, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:d], test_params[:b], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[5×5: off-diagonal = [b, d, b, b]] ← Our known pattern!")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:d], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[5×5: off-diagonal = [b, b, d, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:b], test_params[:d], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[5×5: off-diagonal = [b, b, b, d]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:b], test_params[:b], test_params[:d]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[6×6: off-diagonal = [b, d, b, b, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:d], test_params[:b], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[6×6: off-diagonal = [d, b, b, b, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:d], test_params[:b], test_params[:b], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[7×7: off-diagonal = [b, d, b, b, b, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:d], test_params[:b], test_params[:b], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n" * "="^70)
println("3. SYMMETRIC PERTURBATION PATTERNS")
println("="^70)
println("Pattern: Symmetric placement of different values")

println("\n[6×6: off-diagonal = [b, d, b, b, d, b]] - symmetric")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:d], test_params[:b], test_params[:b], test_params[:d], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[7×7: off-diagonal = [b, d, b, b, b, d, b]] - symmetric")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:d], test_params[:b], test_params[:b], test_params[:b], test_params[:d], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[7×7: off-diagonal = [d, b, b, b, b, b, d]] - symmetric")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:d], test_params[:b], test_params[:b], test_params[:b], test_params[:b], test_params[:b], test_params[:d]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n" * "="^70)
println("4. TWO-PERTURBATION PATTERNS (Adjacent)")
println("="^70)

println("\n[6×6: off-diagonal = [b, d, d, b, b, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:d], test_params[:d], test_params[:b], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n[6×6: off-diagonal = [b, b, d, d, b, b]]")
mat = make_symmetric_tridiagonal(test_params[:a], [test_params[:b], test_params[:b], test_params[:d], test_params[:d], test_params[:b], test_params[:b]])
eigs = sort(real.(eigvals(mat)))
analyze_eigenvalues(eigs, test_params)

println("\n" * "="^70)
println("ANALYSIS COMPLETE")
println("="^70)
println("\nLook for patterns where multiple eigenvalues match simple formulas!")
println("These are candidates for closed-form solutions.")
println("\nKey observations:")
println("  - Constant tridiagonal has known closed-form (trigonometric)")
println("  - [b,d,b,b] in 5×5 has closed-form: a ± √(2b²+d²), a±b, a")
println("  - Check if other patterns show similar structure!")
