module SymbolicDiagonalization

using LinearAlgebra
using Symbolics
using SymbolicUtils
using SymbolicUtils.Rewriters

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
include("simplify.jl")    # Trigonometric and expression simplification rules
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
include("patterns/kronecker.jl")    # Kronecker product detection (generic)

include("patterns/tridiagonal.jl")  # Toeplitz tridiagonal, special 5×5, anti-diagonal
include("patterns/permutation.jl")  # Permutation matrices

# ============================================================================
# Lie Groups and Algebras (modular organization)
# ============================================================================

include("patterns/lie/common.jl")   # Shared detection helpers
include("patterns/lie/SO2.jl")      # SO(2) rotations and Kronecker products
include("patterns/lie/SO3.jl")      # SO(3) rotations and Kronecker products
include("patterns/lie/SO4.jl")      # SO(4) detection and eigenvalues
include("patterns/lie/SU2.jl")      # SU(2), Pauli matrices, Kronecker products
include("patterns/lie/SU3.jl")      # SU(3), Gell-Mann matrices, Kronecker products
include("patterns/lie/Sp.jl")       # Sp(2), Sp(4) symplectic groups
include("patterns/lie/algebras.jl") # Lie algebra representations (spin-j, etc.)

# ============================================================================
# Public API
# ============================================================================

include("diagonalize.jl") # Main API: symbolic_eigenvalues, symbolic_eigenpairs, etc.

# Export original API for advanced use
export characteristic_polynomial, symbolic_roots, symbolic_eigenvalues, symbolic_eigenpairs, symbolic_diagonalize

# Export simplification utilities
export trig_simplify, aggressive_simplify, simplify_eigenvalue, simplify_eigenvalues

# Export exception types for error handling
export ExpressionComplexityError, ComputationTimeoutError

# ============================================================================
# Lie Group Exports
# ============================================================================

# SO(2) - 2D Rotations
export SO2_rotation, SO2_kron, SO2_kron_eigenvalues

# SO(3) - 3D Rotations
export SO3_Rx, SO3_Ry, SO3_Rz

# SU(2) - Special Unitary 2D
export pauli_x, pauli_y, pauli_z
export SU2_Ux, SU2_Uy, SU2_Uz, SU2_kron, SU2_kron_eigenvalues

# SU(3) - Special Unitary 3D
export gellmann_matrices, gellmann_1, gellmann_2, gellmann_3, gellmann_4
export gellmann_5, gellmann_6, gellmann_7, gellmann_8
export SU3_diagonal, SU3_diagonal_trig, SU3_kron, SU3_kron_eigenvalues

# Lie Algebra Generators
export spin_j_generators, so3_generators, su2_generators

# LinearAlgebra.eigen and LinearAlgebra.eigvals are automatically available when LinearAlgebra is imported

end
