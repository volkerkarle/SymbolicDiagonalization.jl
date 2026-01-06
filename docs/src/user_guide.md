# User Guide

This guide covers practical workflows for using SymbolicDiagonalization.jl.

## Basic Usage

### Getting Eigenvalues

```@example guide
using SymbolicDiagonalization, Symbolics, LinearAlgebra
using Main: LaTeX

# Define symbolic variables
@variables a b c

# Create a matrix
M = [a b; b c]

# Get eigenvalues
vals = eigvals(M)
nothing # hide
```

**Matrix:**

```@example guide
LaTeX(M)
```

**Eigenvalues:**

```@example guide
LaTeX(vals)
```

### Getting Eigenvectors

```@example guide
@variables x y

# Symmetric matrix
S = [x y; y x]
E = eigen(S)
nothing # hide
```

**Eigenvectors (columns of E.vectors):**

```@example guide
LaTeX(E.vectors)
```

### Full Diagonalization

For `M = P * D * P⁻¹` decomposition:

```@example guide
P, D, pairs = symbolic_diagonalize(S)
nothing # hide
```

**Diagonal matrix D:**

```@example guide
LaTeX(D)
```

## Working with Larger Matrices

### Block-Diagonal Matrices

The package automatically detects block structure:

```@example guide
@variables a b c d e f

# 4×4 block-diagonal
B = [a  b  0  0;
     b  a  0  0;
     0  0  c  d;
     0  0  d  c]

vals = eigvals(B)
nothing # hide
```

**Matrix:**

```@example guide
LaTeX(B)
```

**Eigenvalues (solved block-by-block):**

```@example guide
LaTeX(vals)
```

### Circulant Matrices

Any n×n circulant matrix has closed-form eigenvalues via the DFT:

```@example guide
@variables a b c d

# 4×4 circulant
C = [a b c d;
     d a b c;
     c d a b;
     b c d a]

vals = eigvals(C)
nothing # hide
```

**Matrix:**

```@example guide
LaTeX(C)
```

**Eigenvalues:**

```@example guide
LaTeX(vals)
```

### Kronecker Products

Eigenvalues of A ⊗ B are products of individual eigenvalues:

```@example guide
@variables α β γ δ

A = [α 0; 0 β]
B = [γ 0; 0 δ]
K = kron(A, B)

vals = eigvals(K)
nothing # hide
```

**The 4×4 Kronecker product:**

```@example guide
LaTeX(K)
```

**Eigenvalues {αγ, αδ, βγ, βδ}:**

```@example guide
LaTeX(vals)
```

## Rotation Matrices

### SO(2): 2D Rotations

```@example guide
@variables θ

R = SO2_rotation(θ)
vals = eigvals(R)
nothing # hide
```

**Rotation matrix:**

```@example guide
LaTeX(R)
```

**Eigenvalues e^{±iθ}:**

```@example guide
LaTeX(vals)
```

### SO(3): 3D Rotations

```@example guide
Rz = SO3_Rz(θ)
vals = eigvals(Rz)
nothing # hide
```

**Rotation about z-axis:**

```@example guide
LaTeX(Rz)
```

**Eigenvalues {1, e^{±iθ}}:**

```@example guide
LaTeX(vals)
```

### Kronecker Products of Rotations

```@example guide
@variables θ φ

K = kron(SO2_rotation(θ), SO2_rotation(φ))
vals = eigvals(K)
nothing # hide
```

**Eigenvalues e^{i(±θ±φ)}:**

```@example guide
LaTeX(vals)
```

## Tips and Best Practices

### 1. Check for Structure First

Before calling `eigvals()`, consider if your matrix has exploitable structure:

```julia
# Good: Circulant structure detected automatically
C = [1 2 3; 3 1 2; 2 3 1]
eigvals(C)  # Instant via DFT

# Slower: Generic 3×3 symmetric
M = [a b c; b d e; c e f]
eigvals(M)  # Uses Cardano formula (still exact, but larger expressions)
```

### 2. Use Constructors for Lie Groups

The built-in constructors ensure proper structure detection:

```julia
# Good: Uses SO2_rotation constructor
R = SO2_rotation(θ)

# Works but less clear: Manual construction
R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
```

### 3. Handle Large Expressions

For complex matrices, expressions can grow large:

```julia
# Set timeout for expensive computations
eigvals(M; timeout=60)  # 60 seconds max

# Limit expression complexity
eigvals(M; max_terms=1000)
```

### 4. Numeric Verification

Verify symbolic results numerically:

```julia
@variables a b
M = [a b; b a]
vals = eigvals(M)

# Substitute numeric values
using Symbolics: substitute
val_numeric = substitute(vals[1], Dict(a => 1.0, b => 2.0))
# Should match: eigvals([1.0 2.0; 2.0 1.0])
```

## Common Patterns

### Tridiagonal Matrices

Symmetric Toeplitz tridiagonal matrices have closed-form eigenvalues:

```@example guide
# Path Laplacian is symmetric Toeplitz tridiagonal
L = path_laplacian(5)
vals = eigvals(L)
nothing # hide
```

**5×5 Path Laplacian eigenvalues:**

```@example guide
LaTeX(vals)
```

### Hadamard Matrices

Sylvester-Hadamard matrices have eigenvalues ±√(2ⁿ):

```@example guide
H = hadamard_matrix(3)  # 8×8
vals = eigvals(H)
nothing # hide
```

**Hadamard eigenvalues (exact):**

```@example guide
LaTeX(vals)
```

### DFT Matrices

Fourier matrices have eigenvalues in {±√n, ±i√n}:

```@example guide
F = dft_matrix(4)
vals = eigvals(F)
nothing # hide
```

**DFT eigenvalues:**

```@example guide
LaTeX(vals)
```

## Next Steps

- [Pattern Library](pattern_library.md) - Complete list of supported patterns
- [Mathematical Background](mathematical_background.md) - Theory behind the algorithms
