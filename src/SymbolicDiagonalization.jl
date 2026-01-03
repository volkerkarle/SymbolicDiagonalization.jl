module SymbolicDiagonalization

using LinearAlgebra
using Symbolics


# ============================================================================
# Configuration Constants
# ============================================================================

"""Default threshold for warning about symbolic variable count."""
const DEFAULT_COMPLEXITY_THRESHOLD = 5

"""Default computation timeout in seconds."""
const DEFAULT_TIMEOUT_SECONDS = 300

"""Default maximum expression size (number of terms/operations)."""
const DEFAULT_MAX_TERMS = 10000

"""Maximum matrix size for adjugate-based eigenvector computation (O(n!) complexity)."""
const MAX_ADJUGATE_SIZE = 3

"""Valid structure hints for eigenvalue computation."""
const VALID_STRUCTURES = (:auto, :hermitian, :symmetric, :unitary, :general, :none, :diagonal, :triangular)

# ============================================================================
# Core Utilities (must be loaded first as other modules depend on them)
# ============================================================================

include("rref.jl")        # RREF and _issymzero
include("charpoly.jl")    # Characteristic polynomial computation
include("roots.jl")       # Polynomial root solvers

# ============================================================================
# Structure Detection and Matrix Analysis
# ============================================================================

include("structure.jl")   # Structure detection, matrix properties, adjugate

# ============================================================================
# Pattern Detection Modules
# Each module provides closed-form eigenvalue formulas for specific patterns
# ============================================================================

include("patterns/graphs.jl")       # Hypercube and strongly regular graphs
include("patterns/circulant.jl")    # Circulant and block circulant matrices
include("patterns/kronecker.jl")    # Kronecker product detection (matrix structure)

include("patterns/tridiagonal.jl")  # Toeplitz tridiagonal, special 5×5, anti-diagonal
include("patterns/permutation.jl")  # Permutation matrices
include("patterns/lie_groups.jl")   # SO(n), SU(n), Sp(2n), Lorentz groups

# ============================================================================
# Public API
# ============================================================================

include("diagonalize.jl") # Main API: symbolic_eigenvalues, symbolic_eigenpairs, etc.

# Export original API for advanced use
export characteristic_polynomial, symbolic_roots, symbolic_eigenvalues, symbolic_eigenpairs, symbolic_diagonalize

# Export exception types for error handling
export ExpressionComplexityError, ComputationTimeoutError



# LinearAlgebra.eigen and LinearAlgebra.eigvals are automatically available when LinearAlgebra is imported

end
