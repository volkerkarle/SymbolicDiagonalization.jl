# ============================================================================
# SU(3) - Special Unitary Group in 3D
# ============================================================================
#
# This file contains SU(3) related functionality:
# - Gell-Mann matrices: gellmann_1 through gellmann_8 (the 8 generators)
# - Detection: _is_SU3(A)
# - Eigenvalues: _SU3_eigenvalues(A) - uses trace-based cubic formula
#
# SU(3) is the group of 3×3 special unitary matrices (unitary with det=1).
# Eigenvalues lie on the unit circle with product = 1.
#
# NOTE: We do NOT provide:
# - Diagonal SU(3) constructors (trivial - eigenvalues are on diagonal)
# - Eigenvector computation (for diagonal: trivial; for general: need nullspace)
# - Kronecker products of diagonal SU(3) (trivial)
#
# The non-trivial eigenvalue computation uses the characteristic polynomial:
#   λ³ - tr(A)·λ² + e₂·λ - 1 = 0
# where e₂ = (tr²(A) - tr(A²))/2, which is solved via the cubic formula.
# ============================================================================

using LinearAlgebra
using Symbolics

# ============================================================================
# Gell-Mann Matrices (Generators of SU(3))
# ============================================================================

"""
    gellmann_1() -> Matrix{Int}

Gell-Mann matrix λ₁:
    [0 1 0]
    [1 0 0]
    [0 0 0]
"""
gellmann_1() = [0 1 0; 1 0 0; 0 0 0]

"""
    gellmann_2() -> Matrix{Complex{Int}}

Gell-Mann matrix λ₂:
    [0 -i 0]
    [i  0 0]
    [0  0 0]
"""
gellmann_2() = [0 -im 0; im 0 0; 0 0 0]

"""
    gellmann_3() -> Matrix{Int}

Gell-Mann matrix λ₃:
    [1  0 0]
    [0 -1 0]
    [0  0 0]
"""
gellmann_3() = [1 0 0; 0 -1 0; 0 0 0]

"""
    gellmann_4() -> Matrix{Int}

Gell-Mann matrix λ₄:
    [0 0 1]
    [0 0 0]
    [1 0 0]
"""
gellmann_4() = [0 0 1; 0 0 0; 1 0 0]

"""
    gellmann_5() -> Matrix{Complex{Int}}

Gell-Mann matrix λ₅:
    [0 0 -i]
    [0 0  0]
    [i 0  0]
"""
gellmann_5() = [0 0 -im; 0 0 0; im 0 0]

"""
    gellmann_6() -> Matrix{Int}

Gell-Mann matrix λ₆:
    [0 0 0]
    [0 0 1]
    [0 1 0]
"""
gellmann_6() = [0 0 0; 0 0 1; 0 1 0]

"""
    gellmann_7() -> Matrix{Complex{Int}}

Gell-Mann matrix λ₇:
    [0 0  0]
    [0 0 -i]
    [0 i  0]
"""
gellmann_7() = [0 0 0; 0 0 -im; 0 im 0]

"""
    gellmann_8() -> Matrix{Float64}

Gell-Mann matrix λ₈:
    [1/√3   0     0   ]
    [  0  1/√3    0   ]
    [  0    0  -2/√3  ]
"""
gellmann_8() = [1/sqrt(3) 0 0; 0 1/sqrt(3) 0; 0 0 -2/sqrt(3)]

"""
    gellmann_matrices() -> Vector{Matrix}

Return all 8 Gell-Mann matrices [λ₁, λ₂, ..., λ₈].

These form a basis for the Lie algebra su(3) (traceless skew-Hermitian 3×3 matrices).
SU(3) elements can be written as exp(i·Σⱼ θⱼλⱼ/2).
"""
gellmann_matrices() = [gellmann_1(), gellmann_2(), gellmann_3(), gellmann_4(), 
                        gellmann_5(), gellmann_6(), gellmann_7(), gellmann_8()]

# ============================================================================
# Detection
# ============================================================================

"""
    _is_SU3(A)

Check if A is a 3×3 special unitary matrix (SU(3)).
Returns true if A ∈ SU(3), false otherwise.
"""
function _is_SU3(A)
    size(A) == (3, 3) || return false
    return !isnothing(_is_special_unitary(A))
end

"""
    _is_SU3_trig(A)

Check if A is an SU(3) matrix using trig-aware simplification.
Useful when matrix entries contain trigonometric expressions.
"""
function _is_SU3_trig(A)
    size(A) == (3, 3) || return false
    
    # Check unitarity: A * A^H = I with trig simplification
    AAH = A * adjoint(A)
    for i in 1:3, j in 1:3
        target = i == j ? 1 : 0
        diff = AAH[i, j] - target
        if !_issymzero_trig(real(diff)) || !_issymzero_trig(imag(diff))
            return false
        end
    end
    
    # Check det = 1 with trig simplification
    d = det(A)
    if !_issymzero_trig(real(d) - 1) || !_issymzero_trig(imag(d))
        return false
    end
    
    return true
end

# ============================================================================
# Eigenvalues
# ============================================================================

"""
    _SU3_eigenvalues(A)

Compute eigenvalues of an SU(3) matrix using the cubic formula.

For SU(3), eigenvalues are on the unit circle with product = 1.
The characteristic polynomial is: λ³ - tr(A)·λ² + e₂·λ - 1 = 0

This is non-trivial even for diagonal matrices if expressed symbolically,
as it extracts the structure from the trace rather than reading diagonal entries.
"""
function _SU3_eigenvalues(A)
    _is_SU3(A) || return nothing
    
    t1 = tr(A)
    t2 = tr(A * A)
    e2 = (t1^2 - t2) / 2
    
    # Characteristic polynomial: λ³ - t1·λ² + e2·λ - 1 = 0
    coeffs = [-1, e2, -t1, 1]
    return symbolic_roots(coeffs)
end
