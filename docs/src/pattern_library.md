# Pattern Library

Complete catalog of all special matrix patterns with closed-form eigenvalue solutions implemented in `SymbolicDiagonalization.jl`.

## Table of Contents

1. [Basic Patterns](#basic-patterns)
   - [Diagonal Matrices](#diagonal-matrices)
   - [Triangular Matrices](#triangular-matrices)
   - [Jordan Blocks](#jordan-blocks)
2. [Structure-Based Patterns](#structure-based-patterns)
   - [Block-Diagonal](#block-diagonal-matrices)
   - [Persymmetric](#persymmetric-matrices)
3. [Symmetry-Based Patterns](#symmetry-based-patterns)
   - [Circulant Matrices](#circulant-matrices)
   - [Block Circulant](#block-circulant-matrices)
   - [Anti-Diagonal](#anti-diagonal-matrices)
4. [Tensor Product Patterns](#tensor-product-patterns)
   - [Kronecker Products](#kronecker-products)
5. [Tridiagonal Patterns](#tridiagonal-patterns)
   - [Symmetric Toeplitz Tridiagonal](#symmetric-toeplitz-tridiagonal)
   - [Special 5×5 Patterns](#special-5×5-tridiagonal-patterns)
6. [Permutation Patterns](#permutation-patterns)
   - [Permutation Matrices](#permutation-matrices)

---

## Basic Patterns

### Diagonal Matrices

**Structure**: Non-zero only on main diagonal

**Matrix Form**:
```math
\mathbf{D} = \begin{bmatrix}
d_1 & 0 & 0 & \cdots & 0 \\
0 & d_2 & 0 & \cdots & 0 \\
0 & 0 & d_3 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & d_n
\end{bmatrix}
```

**Eigenvalues**: $\{d_1, d_2, d_3, \ldots, d_n\}$ (diagonal entries)

**Eigenvectors**: Standard basis vectors $\mathbf{e}_i$

**Complexity**: $O(n)$ - direct reading

**Example**:
```julia
@variables a b c
D = [a 0 0; 0 b 0; 0 0 c]

eigvals(D)  # [a, b, c]
```

**Why this works**: Diagonal matrices are already in eigenvalue form. The eigenvalue equation $\mathbf{D}\mathbf{v} = \lambda\mathbf{v}$ is satisfied by standard basis vectors.

**Mathematical properties**:
- Always diagonalizable
- Eigenvectors are orthogonal
- Simple eigenvalue structure

---

### Triangular Matrices

**Structure**: Upper or lower triangular (zeros below or above diagonal)

**Upper Triangular Form**:
```math
\mathbf{U} = \begin{bmatrix}
u_{11} & u_{12} & u_{13} & \cdots & u_{1n} \\
0 & u_{22} & u_{23} & \cdots & u_{2n} \\
0 & 0 & u_{33} & \cdots & u_{3n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & u_{nn}
\end{bmatrix}
```

**Eigenvalues**: $\{u_{11}, u_{22}, u_{33}, \ldots, u_{nn}\}$ (diagonal entries)

**Eigenvectors**: Computed via back-substitution

**Complexity**: $O(n)$ for eigenvalues, $O(n^2)$ for eigenvectors

**Example**:
```julia
@variables a b c
U = [a 1 2; 0 b 3; 0 0 c]

eigvals(U)  # [a, b, c]
```

**Why this works**: For triangular matrices, $\det(\lambda\mathbf{I} - \mathbf{U})$ is a product of diagonal terms $(\lambda - u_{ii})$, so eigenvalues are exactly the diagonal entries.

**Mathematical properties**:
- May not be diagonalizable (if repeated eigenvalues)
- Eigenvalues independent of off-diagonal entries
- Characteristic polynomial factorizes immediately

---

### Jordan Blocks

**Structure**: Single repeated eigenvalue on diagonal, ones on superdiagonal

**Matrix Form**:
```math
\mathbf{J}(\lambda, n) = \begin{bmatrix}
\lambda & 1 & 0 & \cdots & 0 \\
0 & \lambda & 1 & \cdots & 0 \\
0 & 0 & \lambda & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \lambda
\end{bmatrix}
```

**Eigenvalues**: $\{\lambda\}$ with algebraic multiplicity $n$

**Eigenvectors**: Only ONE eigenvector (not diagonalizable for $n > 1$)

**Complexity**: $O(1)$ for eigenvalues

**Example**:
```julia
@variables λ
J = [λ 1 0; 0 λ 1; 0 0 λ]

eigvals(J)  # [λ, λ, λ]
```

**Why this works**: $\det(t\mathbf{I} - \mathbf{J}) = (t - \lambda)^n$, giving repeated root $\lambda$. Geometric multiplicity is 1 (rank deficiency).

**Mathematical properties**:
- **Not diagonalizable** (unless $n = 1$)
- Prototype of non-diagonalizable matrices
- Appears in Jordan normal form decomposition

**Caution**: `symbolic_diagonalize()` will throw error since Jordan blocks ($n > 1$) are not diagonalizable.

---

## Structure-Based Patterns

### Block-Diagonal Matrices

**Structure**: Non-zero only in diagonal blocks

**Matrix Form**:
```math
\mathbf{B} = \begin{bmatrix}
\mathbf{A}_1 & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \\
\mathbf{0} & \mathbf{A}_2 & \mathbf{0} & \cdots & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{A}_3 & \cdots & \mathbf{0} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{A}_k
\end{bmatrix}
```
where each $\mathbf{A}_i$ is a square matrix (possibly different sizes).

**Eigenvalues**: $\lambda(\mathbf{B}) = \lambda(\mathbf{A}_1) \cup \lambda(\mathbf{A}_2) \cup \cdots \cup \lambda(\mathbf{A}_k)$

**Complexity**: $O(\sum_i f(n_i))$ where $f(n)$ is cost to solve size-$n$ block

**Example**:
```julia
@variables a b c d
B = [a  b  0  0;
     b  a  0  0;
     0  0  c  d;
     0  0  d  c]

eigvals(B)  # [a+b, a-b, c+d, c-d]
```

**Why this works**: The characteristic polynomial factorizes:
```math
\det(\lambda\mathbf{I} - \mathbf{B}) = \det(\lambda\mathbf{I} - \mathbf{A}_1) \cdot \det(\lambda\mathbf{I} - \mathbf{A}_2) \cdot \ldots \cdot \det(\lambda\mathbf{I} - \mathbf{A}_k)
```

So eigenvalues are union of eigenvalues of each block.

**Detection algorithm**:
1. Find connected components of non-zero structure
2. Recursively solve each block independently
3. Concatenate eigenvalues

**Practical importance**: An $8 \times 8$ matrix with four $2 \times 2$ blocks requires solving 4 quadratics, not one degree-8 polynomial!

**Mathematical properties**:
- Diagonalizable iff all blocks are diagonalizable
- Eigenvectors are block-structured
- Extremely common in practice (physical subsystems, symmetry blocks)

---

### Persymmetric Matrices

**Structure**: Symmetric about anti-diagonal: $Q[i,j] = Q[n+1-j, n+1-i]$

**Matrix Form** ($4 \times 4$ example):
```math
\mathbf{P} = \begin{bmatrix}
a & b & c & d \\
b & e & f & c \\
c & f & e & b \\
d & c & b & a
\end{bmatrix}
```

**Eigenvalues**: Computed by splitting into half-sized problems

**Complexity**: $O(f(n/2))$ where $f$ is cost for size $n/2$

**Example**:
```julia
@variables a b c d e f
P = [a b c d; b e f c; c f e b; d c b a]

eigvals(P)  # Splits into two 2×2 problems
```

**Why this works**: Persymmetric structure allows transformation:
```math
\mathbf{Q} = \mathbf{J} \mathbf{P} \mathbf{J} \quad \text{(where } \mathbf{J} \text{ is anti-diagonal identity)}
```

This creates block structure in eigenvector space, splitting problem in half.

**Transformation details**:
- Let $\mathbf{J}$ be the flip matrix (anti-diagonal identity)
- Form symmetric combinations: $(\mathbf{P} \pm \mathbf{JPJ})/2$
- These commute and can be simultaneously diagonalized
- Reduces to two half-sized eigenvalue problems

**Detection**: Check if $Q[i,j] = Q[n+1-j, n+1-i]$ for all $i,j$

**Mathematical properties**:
- Always diagonalizable (if symmetric)
- Eigenspaces split into symmetric/antisymmetric subspaces
- Often appears in signal processing, centrosymmetric matrices

**Limitation**: Current implementation handles some but not all persymmetric cases.

---

## Symmetry-Based Patterns

### Circulant Matrices

**Structure**: Each row is cyclic shift of previous row

**Matrix Form**:
```math
\mathbf{C} = \begin{bmatrix}
c_0 & c_1 & c_2 & \cdots & c_{n-1} \\
c_{n-1} & c_0 & c_1 & \cdots & c_{n-2} \\
c_{n-2} & c_{n-1} & c_0 & \cdots & c_{n-3} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
c_1 & c_2 & c_3 & \cdots & c_0
\end{bmatrix}
```

**Eigenvalues** (closed-form for any $n$):
```math
\lambda_k = \sum_{j=0}^{n-1} c_j \cdot \omega^{kj} \quad \text{for } k = 0, 1, \ldots, n-1
```
where $\omega = \exp(2\pi i/n)$ is the $n$th root of unity.

**Eigenvectors**: Columns of DFT matrix (independent of entries!)

**Complexity**: $O(n)$ to compute all $n$ eigenvalues

**Example** ($3 \times 3$):
```julia
@variables a b c
C = [a b c; c a b; b c a]

# ω = exp(2πi/3)
eigvals(C)  # [a+b+c, a+b*ω+c*ω², a+b*ω²+c*ω]
```

**Why this works**: Circulant matrices commute with cyclic shifts. The DFT diagonalizes ALL circulant matrices simultaneously. The eigenvalues are exactly the Discrete Fourier Transform of the first row.

**Mathematical foundation**:
- Circulant matrices form a commutative algebra
- All share the same eigenvectors (DFT basis)
- Multiplication in circulant space = convolution in sequence space
- Diagonalization in DFT space = pointwise multiplication

**Detection algorithm**:
1. Check if $C[i,j] = C[(i-1) \bmod n, (j-1) \bmod n]$ for all $i,j$
2. Extract first row $[c_0, c_1, \ldots, c_{n-1}]$
3. Apply DFT formula

**Applications**:
- Signal processing (circular convolution)
- Time series analysis (periodic processes)
- Graph Laplacians (cycle graphs)
- Coding theory

**Practical notes**:
- Works for ANY size $n$ (even $n = 100+$!)
- Generic $n \times n$ requires degree-$n$ polynomial (impossible for $n \geq 5$)
- Circulant structure bypasses Abel-Ruffini limitation

**Variants**:
- **Skew-circulant**: Last element negated on each shift
- **g-circulant**: Generalized with group operation
- **Block circulant**: See next section

---

### Block Circulant Matrices

**Structure**: Block version of circulant (blocks shift instead of elements)

**Matrix Form** (4 blocks, size $k \times k$ each):
```math
\mathbf{BC} = \begin{bmatrix}
\mathbf{A}_0 & \mathbf{A}_1 & \mathbf{A}_2 & \mathbf{A}_3 \\
\mathbf{A}_3 & \mathbf{A}_0 & \mathbf{A}_1 & \mathbf{A}_2 \\
\mathbf{A}_2 & \mathbf{A}_3 & \mathbf{A}_0 & \mathbf{A}_1 \\
\mathbf{A}_1 & \mathbf{A}_2 & \mathbf{A}_3 & \mathbf{A}_0
\end{bmatrix}
```

**Eigenvalues**: Compute by solving $n$ eigenvalue problems of size $k \times k$:
```math
\lambda(\mathbf{BC}) = \bigcup_{m=0}^{n-1} \lambda\left(\sum_{j=0}^{n-1} \mathbf{A}_j \cdot \omega^{mj}\right)
```
where $\omega = \exp(2\pi i/n)$.

**Complexity**: $O(n \cdot f(k))$ where $f(k)$ is cost for $k \times k$ matrix

**Example** ($4 \times 4$ with $2 \times 2$ blocks):
```julia
@variables a b c d
A = [a b; c d]
B = [1 0; 0 1]

# 4×4 block circulant [A B; B A]
BC = [a b 1 0;
      c d 0 1;
      1 0 a b;
      0 1 c d]

eigvals(BC)  # Solves 2 problems of size 2×2
```

**Why this works**: Block circulant matrices are diagonalized by block DFT. The reduction formula creates $n$ different $k \times k$ matrices (linear combinations of blocks with DFT coefficients), which can be solved independently.

**Mathematical foundation**:
- Generalization of circulant to matrix entries
- Block DFT plays same role as regular DFT
- Eigenvalue problem decouples into $n$ independent $k \times k$ problems

**Detection algorithm**:
1. Check block structure (equal-sized blocks)
2. Verify block circulant property
3. Extract blocks $[\mathbf{A}_0, \mathbf{A}_1, \ldots, \mathbf{A}_{n-1}]$
4. Form $n$ combinations with DFT weights
5. Solve each $k \times k$ system

**Key examples**:
- **$12 \times 12$ with $3 \times 3$ blocks**: 4 blocks → solve 4 cubic equations (feasible)
- **$8 \times 8$ with $4 \times 4$ blocks**: 2 blocks → solve 2 quartic equations (feasible but large)
- **$16 \times 16$ with $2 \times 2$ blocks**: 8 blocks → solve 8 quadratic equations (easy)

Without block structure, these would require degree-12, degree-8, and degree-16 polynomials (impossible)!

**Applications**:
- Multi-channel signal processing
- Block Toeplitz systems
- Vectorized circulant operations
- Kronecker-structured problems

**Limitation**: Requires $k \leq 4$ for closed form (unless further structure in blocks).

---

### Anti-Diagonal Matrices

**Structure**: Non-zero only on anti-diagonal, with symmetry

**Matrix Form** ($n \times n$, symmetric about anti-diagonal):
```math
\mathbf{A} = \begin{bmatrix}
0 & 0 & 0 & \cdots & 0 & a_n \\
0 & 0 & 0 & \cdots & a_{n-1} & 0 \\
0 & 0 & 0 & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
a_{n-1} & 0 & 0 & \cdots & 0 & 0 \\
a_n & 0 & 0 & \cdots & 0 & 0
\end{bmatrix}
```
with symmetry: $A[i,j] = A[n+1-j, n+1-i]$

**Eigenvalues** (closed-form for any $n$):
- **Odd $n$**: $\{a_k, -a_k\}$ for $k = 1, \ldots, (n-1)/2$, plus $a_{(n+1)/2}$ (center element)
- **Even $n$**: $\{a_k, -a_k\}$ for $k = 1, \ldots, n/2$

**Eigenvector structure**: $\pm$ pairs with specific symmetry

**Complexity**: $O(n)$ - direct computation

**Example** ($5 \times 5$):
```julia
@variables a b c
A = [0 0 0 0 a;
     0 0 0 b 0;
     0 0 c 0 0;
     0 b 0 0 0;
     a 0 0 0 0]

eigvals(A)  # {c, ±a, ±b}
```

**Why this works**: The persymmetric symmetry combined with anti-diagonal structure forces eigenvalues into $\pm$ pairs. The transformation $\mathbf{J}\mathbf{A}$ (flip matrix) equals $\pm\mathbf{A}$, creating symmetric/antisymmetric subspaces.

**Mathematical foundation**:
- Anti-diagonal flip symmetry: $\mathbf{JAJ}^\dagger = \pm\mathbf{A}$
- Eigenvectors come in symmetric/antisymmetric pairs
- Eigenvalue equation: If $\mathbf{A}\mathbf{v} = \lambda\mathbf{v}$, then $\mathbf{A}(\mathbf{J}\mathbf{v}) = -\lambda(\mathbf{J}\mathbf{v})$

**Detection algorithm**:
1. Check if zeros everywhere except anti-diagonal
2. Verify symmetry: $A[i, n+1-i] = A[n+1-i, i]$
3. Extract anti-diagonal elements
4. Form $\pm$ pairs (center unpaired if odd $n$)

**Applications**:
- Exchange matrices
- Reversal operators
- Symmetric perturbations of identity
- Quantum mechanics (parity operators)

**Mathematical properties**:
- Always diagonalizable
- Real eigenvalues if entries are real
- Orthogonal eigenvectors (pairs)

---

## Tensor Product Patterns

### Kronecker Products

**Structure**: Tensor product $\mathbf{A} \otimes \mathbf{B} = [a_{ij}\mathbf{B}]$

**Matrix Form** ($2 \times 2 \otimes 2 \times 2$ example):
```math
\mathbf{A} \otimes \mathbf{B} = \begin{bmatrix}
a_{11}\mathbf{B} & a_{12}\mathbf{B} \\
a_{21}\mathbf{B} & a_{22}\mathbf{B}
\end{bmatrix}
```

where each entry $a_{ij}\mathbf{B}$ is a block.

**Eigenvalues** (closed-form for any sizes $m$, $n$):
```math
\lambda(\mathbf{A} \otimes \mathbf{B}) = \{\lambda_i(\mathbf{A}) \cdot \lambda_j(\mathbf{B}) : i = 1\ldots m, j = 1\ldots n\}
```

All $m \cdot n$ products of eigenvalues of $\mathbf{A}$ and $\mathbf{B}$.

**Eigenvectors**: $\mathbf{v}_{ij} = \mathbf{v}_i(\mathbf{A}) \otimes \mathbf{v}_j(\mathbf{B})$ (tensor products)

**Complexity**: $O(f(m) + f(n))$ where $f$ is cost for each factor

**Example** ($4 \times 4 = 2 \times 2 \otimes 2 \times 2$):
```julia
@variables a b c d
A = [a 0; 0 b]  # eigenvalues {a, b}
B = [c 0; 0 d]  # eigenvalues {c, d}

K = kron(A, B)  # 4×4 matrix

eigvals(K)  # {ac, ad, bc, bd}
```

**Why this works**: The eigenvalue equation for Kronecker products:
```math
(\mathbf{A} \otimes \mathbf{B})(\mathbf{u} \otimes \mathbf{v}) = (\mathbf{A}\mathbf{u}) \otimes (\mathbf{B}\mathbf{v}) = (\lambda_u\mathbf{u}) \otimes (\lambda_v\mathbf{v}) = (\lambda_u\lambda_v)(\mathbf{u} \otimes \mathbf{v})
```

So if $\mathbf{u}$ is eigenvector of $\mathbf{A}$ with eigenvalue $\lambda_u$, and $\mathbf{v}$ is eigenvector of $\mathbf{B}$ with eigenvalue $\lambda_v$, then $\mathbf{u}\otimes\mathbf{v}$ is eigenvector of $\mathbf{A}\otimes\mathbf{B}$ with eigenvalue $\lambda_u\lambda_v$.

**Mathematical foundation**:
- Tensor product structure in linear algebra
- Factorization of eigenspaces: $E(\mathbf{A}\otimes\mathbf{B}) = E(\mathbf{A}) \otimes E(\mathbf{B})$
- Multiplicative property of eigenvalues

**Detection algorithm**:
1. Check block pattern: equal-sized blocks
2. Verify scaling relationship between blocks
3. Extract factor matrices $\mathbf{A}$ and $\mathbf{B}$
4. Solve eigenvalues of $\mathbf{A}$ and $\mathbf{B}$ independently
5. Form all products $\lambda_i(\mathbf{A}) \cdot \lambda_j(\mathbf{B})$

**Key examples**:
- **$6 \times 6 = 2 \times 2 \otimes 3 \times 3$**: Quadratic $\times$ cubic = feasible
- **$12 \times 12 = 3 \times 3 \otimes 4 \times 4$**: Cubic $\times$ quartic = feasible (but large)
- **$12 \times 12 = 4 \times 4 \otimes 3 \times 3$**: Same as above

Without Kronecker structure, degree-6 and degree-12 would be impossible!

**Applications**:
- Multi-dimensional systems (separable PDEs)
- Quantum mechanics (composite systems)
- Graph products (Cartesian product graphs)
- Multi-way tensor analysis

**Detection challenges**:
- Requires exact scaling pattern recognition
- Floating-point errors can break detection
- Current implementation: basic pattern matching

---

## Tridiagonal Patterns

### Symmetric Toeplitz Tridiagonal

**Structure**: Tridiagonal with constant diagonals (Toeplitz + tridiagonal)

**Matrix Form**:
```math
\mathbf{T} = \begin{bmatrix}
a & b & 0 & 0 & \cdots & 0 \\
b & a & b & 0 & \cdots & 0 \\
0 & b & a & b & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots
\end{bmatrix}
```

**Eigenvalues** (closed-form for any $n$):
```math
\lambda_k = a + 2b\cos\left(\frac{k\pi}{n+1}\right) \quad \text{for } k = 1, 2, \ldots, n
```

**Eigenvectors** (known analytically):
```math
v_k(j) = \sin\left(\frac{jk\pi}{n+1}\right) \quad \text{for } j = 1, \ldots, n
```

**Complexity**: $O(n)$ to compute all eigenvalues

**Example** ($4 \times 4$):
```julia
@variables a b
T = [a b 0 0;
     b a b 0;
     0 b a b;
     0 0 b a]

# λₖ = a + 2b·cos(kπ/5) for k=1,2,3,4
eigvals(T)
```

**Why this works**: The eigenvectors are discrete sine functions, which are known eigenfunctions of the discrete Laplacian operator. The eigenvalues come from the cosine terms in the eigenvalue equation.

**Mathematical foundation**:
- Related to discrete Laplacian operator
- Eigenvectors are Discrete Sine Transform (DST) basis
- Connection to Chebyshev polynomials
- Comes from finite-difference discretization of $d^2/dx^2$

**Derivation**:
The eigenvalue equation $\mathbf{T}\mathbf{v} = \lambda\mathbf{v}$ with $\mathbf{v} = [\sin(k\pi/(n+1)), \sin(2k\pi/(n+1)), \ldots, \sin(nk\pi/(n+1))]$ gives:
```math
a\sin\left(\frac{jk\pi}{n+1}\right) + b\sin\left(\frac{(j-1)k\pi}{n+1}\right) + b\sin\left(\frac{(j+1)k\pi}{n+1}\right) = \lambda_k\sin\left(\frac{jk\pi}{n+1}\right)
```

Using trigonometric identity $\sin(\alpha)+\sin(\beta) = 2\sin((\alpha+\beta)/2)\cos((\alpha-\beta)/2)$, this simplifies to $\lambda_k = a + 2b\cos(k\pi/(n+1))$.

**Detection algorithm**:
1. Check tridiagonal structure (zeros beyond sub/super-diagonals)
2. Verify constant diagonals: all diagonal = $a$, all off-diagonal = $b$
3. Verify symmetry: sub-diagonal = super-diagonal
4. Apply formula for $k = 1, \ldots, n$

**Applications**:
- Finite difference methods (1D Laplacian)
- Vibrating string (discrete model)
- Heat equation discretization
- Graph Laplacian (path graph)
- Quantum mechanics (tight-binding model)

**Generalizations**:
- **Non-symmetric**: Different super/sub-diagonals (no closed form in general)
- **Periodic**: Adds corner elements (circulant tridiagonal)
- **Non-Toeplitz**: Variable diagonals (generally no closed form)

**Why it scales to any $n$**: The eigenvector basis is KNOWN and independent of matrix size. We don't need to solve polynomial—we directly compute eigenvalues from formula.

---

### Special 5×5 Tridiagonal Patterns

**Structure**: Tridiagonal $5 \times 5$ with one "perturbation" from constant pattern

Two specific patterns have been discovered with closed-form eigenvalues.

#### Pattern 1: [b, d, b, b]

**Matrix Form**:
```math
\mathbf{M} = \begin{bmatrix}
a & b & 0 & 0 & 0 \\
b & a & d & 0 & 0 \\
0 & d & a & b & 0 \\
0 & 0 & b & a & b \\
0 & 0 & 0 & b & a
\end{bmatrix}
```

Off-diagonal sequence: $[b, d, b, b]$

**Eigenvalues** (closed-form):
```math
\left\{a + \sqrt{2b^2 + d^2}, \; a - \sqrt{2b^2 + d^2}, \; a + b, \; a - b, \; a\right\}
```

**Example**:
```julia
@variables a b d
M = [a b 0 0 0; b a d 0 0; 0 d a b 0; 0 0 b a b; 0 0 0 b a]

eigvals(M)  # {a ± √(2b² + d²), a ± b, a}
```

#### Pattern 2: [b, b, d, b]

**Matrix Form**:
```math
\mathbf{M} = \begin{bmatrix}
a & b & 0 & 0 & 0 \\
b & a & b & 0 & 0 \\
0 & b & a & d & 0 \\
0 & 0 & d & a & b \\
0 & 0 & 0 & b & a
\end{bmatrix}
```

Off-diagonal sequence: $[b, b, d, b]$

**Eigenvalues** (closed-form - SAME as Pattern 1!):
```math
\left\{a + \sqrt{2b^2 + d^2}, \; a - \sqrt{2b^2 + d^2}, \; a + b, \; a - b, \; a\right\}
```

**Example**:
```julia
@variables a b d
M = [a b 0 0 0; b a b 0 0; 0 b a d 0; 0 0 d a b; 0 0 0 b a]

eigvals(M)  # {a ± √(2b² + d²), a ± b, a}  (identical to Pattern 1!)
```

**Remarkable discovery**: Despite different positions of perturbation $d$, both patterns yield IDENTICAL eigenvalues!

**Why this works**: Unknown! This is an empirical discovery. The mathematical explanation likely involves:
- Hidden symmetry in characteristic polynomial
- Isospectrality (different matrices, same spectrum)
- Possibly related to Chebyshev polynomial structure

**What doesn't work** (patterns without closed form):
- **$[d, b, b, b]$**: Boundary perturbation at position 0
- **$[b, b, b, d]$**: Boundary perturbation at position 3
- **$[d, d, b, b]$**, **$[b, d, d, b]$**: Multiple perturbations
- **$[c, d, e, f]$**: General perturbations

**Key insight from research**: Interior perturbations (positions 1-2) admit closed forms, but boundary perturbations (positions 0, 3) break the pattern.

**Detection algorithm**:
1. Check $5 \times 5$ tridiagonal structure
2. Extract super-diagonal sequence
3. Match against known patterns $[b,d,b,b]$ or $[b,b,d,b]$
4. Verify constant diagonal = $a$
5. Apply closed-form formula

**Applications**:
- Perturbed vibrating strings
- Defects in 1D lattice models
- Specific quantum systems

**Research opportunity**: 
- Why do these specific patterns work?
- Are there similar $7 \times 7$ or $9 \times 9$ patterns?
- Can we characterize all solvable perturbation patterns?

See `notes/PATTERN_DISCOVERIES.md` for more details on discovery methodology.

---

## Permutation Patterns

### Permutation Matrices

**Structure**: Exactly one 1 in each row and column, rest zeros

**Matrix Form** (example: permutation $1 \to 2$, $2 \to 3$, $3 \to 1$):
```math
\mathbf{P} = \begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{bmatrix}
```

**Eigenvalues** (closed-form for any permutation):
- Each $k$-cycle contributes $k$-th roots of unity: $\{1, \omega, \omega^2, \ldots, \omega^{k-1}\}$ where $\omega = \exp(2\pi i/k)$
- All eigenvalues have magnitude 1 (on unit circle)

**Complexity**: $O(n)$ to compute cycle decomposition and roots

**Example** ($6 \times 6$ with 3-cycle, 2-cycle, fixed point):
```julia
# Permutation: (1→2→3→1), (4↔5), (6)
P = [0 1 0 0 0 0;   # 1→2
     0 0 1 0 0 0;   # 2→3
     1 0 0 0 0 0;   # 3→1
     0 0 0 0 1 0;   # 4→5
     0 0 0 1 0 0;   # 5→4
     0 0 0 0 0 1]   # 6→6

eigvals(P)
# 3-cycle contributes: {1, ω, ω²} where ω = exp(2πi/3)
# 2-cycle contributes: {1, -1}
# Fixed point contributes: {1}
# Total: {1, 1, 1, -1, ω, ω²}
```

**Why this works**: 
1. Permutation matrices have finite order: $\mathbf{P}^k = \mathbf{I}$ for some $k$
2. If $\mathbf{P}^k\mathbf{v} = \mathbf{v}$, then $\mathbf{P}$ has eigenvalues that are $k$-th roots of unity
3. Cycle decomposition determines which roots appear

**Mathematical foundation**:
- **Cayley's theorem**: Every permutation decomposes into disjoint cycles
- **Root of unity property**: If $\mathbf{P}^k = \mathbf{I}$, then $\lambda^k = 1$, so $\lambda$ is a $k$-th root of unity
- **Cycle length = root order**: $k$-cycle gives $k$-th roots of unity

**Cycle decomposition algorithm**:
1. Find permutation: which column has 1 in each row
2. Follow cycle: start at 1, go to $P(1)$, then $P(P(1))$, until back to 1
3. Repeat for unvisited elements
4. Count cycle lengths

**Eigenvalue computation**:
For each $k$-cycle, add $k$ eigenvalues $\{\exp(2\pi ij/k) : j = 0, \ldots, k-1\}$

**Applications**:
- Sorting algorithms (permutation analysis)
- Group theory (symmetric group representations)
- Graph automorphisms
- Coding theory (permutation codes)
- Fourier analysis on symmetric group

**Mathematical properties**:
- Always orthogonal: $\mathbf{P}^\dagger\mathbf{P} = \mathbf{I}$
- Eigenvalues on unit circle: $|\lambda| = 1$
- Diagonalizable (orthogonally)
- Determinant = $\pm 1$ (sign of permutation)

**Special cases**:
- **Identity**: One $n$-cycle → eigenvalue 1 with multiplicity $n$
- **Full cycle**: One $n$-cycle → all $n$-th roots of unity $\{\exp(2\pi ik/n) : k=0\ldots n-1\}$
- **Transposition**: (2-cycle) → eigenvalues $\{1, -1\}$
- **Fixed points**: Contribute eigenvalue 1

**Detection algorithm**:
1. Check if matrix has exactly one 1 per row and column
2. All other entries must be 0
3. Extract permutation mapping
4. Compute cycle decomposition
5. Generate roots of unity for each cycle

---

## Summary Table

| Pattern | Size Limit | Complexity | Key Property |
|---------|-----------|------------|--------------|
| Diagonal | Any $n$ | $O(n)$ | Direct reading |
| Triangular | Any $n$ | $O(n)$ | Diagonal entries |
| Jordan Block | Any $n$ | $O(1)$ | Single eigenvalue |
| Block-Diagonal | Any $n$ | $\sum O(n_i)$ | Union of blocks |
| Persymmetric | Any $n$ | $O(n/2)$ | Half-size split |
| Circulant | Any $n$ | $O(n)$ | DFT formula |
| Block Circulant | Any $n$ ($k \leq 4$) | $O(n \cdot k)$ | Block DFT |
| Kronecker $\mathbf{A}\otimes\mathbf{B}$ | Any $m,n$ ($\leq 4$) | $O(m+n)$ | Product rule |
| Toeplitz Tridiag | Any $n$ | $O(n)$ | Cosine formula |
| Anti-Diagonal | Any $n$ | $O(n)$ | $\pm$ pairs |
| Permutation | Any $n$ | $O(n)$ | Roots of unity |
| Special $5 \times 5$ [b,d,b,b] | $5 \times 5$ only | $O(1)$ | Empirical |
| Special $5 \times 5$ [b,b,d,b] | $5 \times 5$ only | $O(1)$ | Empirical |

---

## Pattern Discovery Resources

Additional resources available in the repository:

- **PATTERN_DISCOVERIES.md** (`notes/PATTERN_DISCOVERIES.md`) - Detailed discovery notes and mathematical justifications
- **DISCOVERY_METHODOLOGY.md** (`notes/DISCOVERY_METHODOLOGY.md`) - How to discover new patterns
- **RESEARCH_SUMMARY.md** (`notes/RESEARCH_SUMMARY.md`) - Summary of research findings
- **explore_patterns.jl** (`examples/explore_patterns.jl`) - Interactive pattern exploration tools

---

## Contributing New Patterns

If you discover a new pattern, please:

1. Verify numerically (at least 3 test cases)
2. Verify symbolically (check characteristic polynomial)
3. Explain WHY it works (mathematical justification)
4. Implement detector and solver
5. Add comprehensive tests
6. Document in this file

See [Contributing Guide](contributing.md) for details.
