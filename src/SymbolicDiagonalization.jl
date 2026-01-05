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
include("roots.jl")       # Polynomial root solvers (degrees 1-4)

# ============================================================================
# Structure Detection and Matrix Analysis
# ============================================================================

include("structure.jl")   # Structure detection, matrix properties, adjugate

# ============================================================================
# Finite Group Patterns
# Patterns arising from matrices invariant under finite group actions.
# These have closed-form eigenvalues due to solvable Galois groups.
# ============================================================================

# Cyclic group Zₙ: Circulant matrices (DFT diagonalizes)
include("patterns/finite/circulant.jl")

# Dihedral group Dₙ: Symmetric circulant matrices (real eigenvalues via cosine formula)
include("patterns/finite/dihedral.jl")

# Quaternion group Q₈: Quaternion-structured matrices (eigenvalues from norm)
include("patterns/finite/quaternion.jl")

# Symmetric group Sₙ: Permutation matrices (roots of unity from cycles)
include("patterns/finite/permutation.jl")

# (Z₂)ⁿ and strongly regular: Graph patterns with algebraic eigenvalues
include("patterns/finite/graphs.jl")

# Coxeter/Weyl groups: Reflection groups with root system eigenvalues
include("patterns/finite/coxeter.jl")

# Hadamard and DFT matrices: Orthogonal/unitary with closed-form eigenvalues
include("patterns/finite/hadamard.jl")

# Special Toeplitz: Anti-circulant, Kac-Murdock-Szegő
include("patterns/finite/toeplitz_special.jl")

# ============================================================================
# Compact Lie Group Patterns
# Patterns arising from continuous symmetry groups (rotations, unitaries).
# Eigenvalues on unit circle, computed via trace invariants.
# ============================================================================

# Common detection utilities (orthogonal, unitary, symplectic checks)
include("patterns/lie/common.jl")

# SO(2): 2D rotations, e^{±iθ}
include("patterns/lie/SO2.jl")

# SO(3): 3D rotations, {1, e^{±iθ}}
include("patterns/lie/SO3.jl")

# SO(4): 4D double rotations, {e^{±iθ₁}, e^{±iθ₂}}
include("patterns/lie/SO4.jl")

# SU(2): Special unitary 2D, e^{±iθ/2} (spin-1/2)
include("patterns/lie/SU2.jl")

# SU(3): Special unitary 3D, cubic formula
include("patterns/lie/SU3.jl")

# Sp(2n): Symplectic groups, reciprocal eigenvalue pairs
include("patterns/lie/Sp.jl")

# Lie algebra representations (spin-j, etc.)
include("patterns/lie/algebras.jl")

# ============================================================================
# Tensor Product Patterns
# Kronecker products A ⊗ B where eigenvalues are products of factor eigenvalues.
# For Lie group Kronecker products, see SO2.jl, SU2.jl (clean trig forms).
# ============================================================================

include("patterns/kronecker.jl")

# ============================================================================
# Structural Patterns
# Non-group patterns with solvable structure (tridiagonal, anti-diagonal).
# ============================================================================

include("patterns/tridiagonal.jl")

# ============================================================================
# Public API
# ============================================================================

include("diagonalize.jl") # Main API: symbolic_eigenvalues, symbolic_eigenpairs, etc.
include("latex.jl")       # LaTeX display wrapper for pretty output

# Export original API for advanced use
export characteristic_polynomial, symbolic_roots, symbolic_eigenvalues, symbolic_eigenpairs, symbolic_diagonalize

# Export LaTeX display wrapper
export LaTeX

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

# SU(3) - Special Unitary 3D (Gell-Mann matrices only, no trivial diagonal constructors)
export gellmann_matrices, gellmann_1, gellmann_2, gellmann_3, gellmann_4
export gellmann_5, gellmann_6, gellmann_7, gellmann_8

# Lie Algebra Generators
export spin_j_generators, so3_generators, su2_generators

# ============================================================================
# Coxeter/Weyl Group Exports
# ============================================================================

# Cartan matrix constructors (all classical and exceptional types)
export cartan_matrix_A, cartan_matrix_B, cartan_matrix_C, cartan_matrix_D, cartan_matrix_E
export cartan_matrix_F4, cartan_matrix_G2

# Coxeter group theory
export coxeter_number, coxeter_exponents, coxeter_element_eigenvalues

# Graph Laplacians
export path_laplacian, cycle_laplacian

# Reflection matrices
export householder_reflection

# ============================================================================
# Structured Matrix Exports
# ============================================================================

# Hadamard and DFT matrices
export hadamard_matrix, dft_matrix

# Anti-circulant and Kac-Murdock-Szegő matrices
export anticirculant_matrix, kms_matrix

# Q₈ group algebra matrices
export Q8_invariant_matrix

# Companion matrices (Frobenius form)
export companion_matrix

# LinearAlgebra.eigen and LinearAlgebra.eigvals are automatically available when LinearAlgebra is imported

end
