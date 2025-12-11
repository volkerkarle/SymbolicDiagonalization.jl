using Test
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

@testset "SymbolicDiagonalization.jl" begin
    @testset "Basic Functionality" begin
        include("test_basic.jl")
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
end
