using Test
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra
using Aqua

@testset "SymbolicDiagonalization.jl" begin
    @testset "Aqua.jl Quality Assurance" begin
        # Test for common package quality issues
        Aqua.test_all(SymbolicDiagonalization;
            ambiguities = false,  # Disabled: Symbolics has many ambiguities we can't control
            stale_deps = false,   # Disabled: Test deps handled separately
            deps_compat = false,  # Disabled: LinearAlgebra is stdlib
            piracies = false,     # Disabled: Intentionally extending eigen/eigvals for Num matrices
        )
    end
    
    @testset "Basic Functionality" begin
        include("test_basic.jl")
    end
    
    @testset "Edge Cases and Error Handling" begin
        include("test_edge_cases.jl")
    end
    
    @testset "Structure Detection" begin
        include("test_structure.jl")
    end
    
    @testset "Special Patterns" begin
        include("test_patterns.jl")
    end
    
    @testset "LinearAlgebra Interface" begin
        include("test_interface.jl")
    end
    
    @testset "Lie Group Patterns" begin
        include("test_lie_groups.jl")
    end
    
    @testset "Lie Algebra Representations" begin
        include("test_lie_algebras.jl")
    end
end
