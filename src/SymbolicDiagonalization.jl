module SymbolicDiagonalization

using LinearAlgebra
using Symbolics

include("rref.jl")
include("charpoly.jl")
include("roots.jl")
include("diagonalize.jl")

# Export original API for advanced use
export characteristic_polynomial, symbolic_roots, symbolic_eigenvalues, symbolic_eigenpairs, symbolic_diagonalize

# Export exception types for error handling
export ExpressionComplexityError, ComputationTimeoutError

# LinearAlgebra.eigen and LinearAlgebra.eigvals are automatically available when LinearAlgebra is imported

end
