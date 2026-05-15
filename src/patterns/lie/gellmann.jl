# ============================================================================
# Gell-Mann Matrices (Generators of SU(3))
# ============================================================================
#
# The Gell-Mann matrices λ₁,...,λ₈ form a basis for the Lie algebra su(3).
# SU(3) elements can be written as exp(i·Σⱼ θⱼλⱼ/2).
#
# These were previously in SU3.jl alongside SU(3) eigenvalue computation.
# The SU(3) eigenvalue computation was removed as it added minimal value over
# the general cubic solver — the only unique contribution was a 3-line
# trace-to-coefficient formula.
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
"""
gellmann_matrices() = [gellmann_1(), gellmann_2(), gellmann_3(), gellmann_4(),
                        gellmann_5(), gellmann_6(), gellmann_7(), gellmann_8()]
