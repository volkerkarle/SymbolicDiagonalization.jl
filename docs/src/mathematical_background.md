# Mathematical Background

This document explains the mathematical foundations underlying `SymbolicDiagonalization.jl`.

## Table of Contents

- [The Eigenvalue Problem](#the-eigenvalue-problem)
- [The Abel-Ruffini Theorem](#the-abel-ruffini-theorem)
- [Why Structure Matters](#why-structure-matters)
- [Closed-Form Root Formulas](#closed-form-root-formulas)
- [Special Pattern Theory](#special-pattern-theory)

## The Eigenvalue Problem

### Definition

Given an $n \times n$ matrix $\mathbf{A}$, we seek scalars $\lambda$ (eigenvalues) and non-zero vectors $\mathbf{v}$ (eigenvectors) satisfying:

```math
\mathbf{A} \mathbf{v} = \lambda \mathbf{v}
```

Equivalently, we need:
```math
(\mathbf{A} - \lambda \mathbf{I}) \mathbf{v} = \mathbf{0}
```

For non-trivial solutions, the matrix $(\mathbf{A} - \lambda \mathbf{I})$ must be singular:
```math
\det(\mathbf{A} - \lambda \mathbf{I}) = 0
```

This determinant is a polynomial in $\lambda$ called the **characteristic polynomial**:
```math
p(\lambda) = \det(\lambda \mathbf{I} - \mathbf{A}) = \lambda^n + c_{n-1}\lambda^{n-1} + \cdots + c_1\lambda + c_0
```

**The eigenvalue problem reduces to finding roots of a polynomial.**

### Why Eigenvalues Matter

Eigenvalues reveal fundamental properties of linear transformations:

- **Dynamics**: Matrix powers $\mathbf{A}^n$ grow/decay based on eigenvalue magnitudes
- **Stability**: System stability determined by eigenvalue real parts
- **Geometry**: Eigenvalues measure stretching/compression along eigenvector directions
- **Spectral decomposition**: $\mathbf{A} = \mathbf{P}\mathbf{D}\mathbf{P}^{-1}$ where $\mathbf{D}$ is diagonal (if $\mathbf{A}$ is diagonalizable)

### The Fundamental Challenge

For general matrices:
- $2 \times 2$ → degree-2 polynomial (quadratic formula)
- $3 \times 3$ → degree-3 polynomial (cubic formula)
- $4 \times 4$ → degree-4 polynomial (quartic formula)
- $5 \times 5$ → degree-5 polynomial (**no general formula!**)
- $n \times n$ → degree-$n$ polynomial (**no general formula for $n \geq 5$**)

This limitation is not computational—it's a fundamental impossibility proven by the Abel-Ruffini theorem.

## The Abel-Ruffini Theorem

### Statement

**Theorem** (Abel-Ruffini, 1824): There is no general algebraic solution (using radicals) for polynomial equations of degree 5 or higher.

### What This Means

- **Degree $\leq 4$**: Closed-form formulas exist using $+$, $-$, $\times$, $\div$, and $n$th roots
- **Degree $\geq 5$**: No general formula exists using these operations
- **Specific polynomials**: May still have closed-form solutions (e.g., $\lambda^5 - 1 = 0$)

### Implications for Symbolic Eigenvalues

✅ **Can solve symbolically**:
- All $1 \times 1$, $2 \times 2$, $3 \times 3$, $4 \times 4$ general matrices
- Larger matrices with exploitable structure

❌ **Cannot solve symbolically** (in general):
- Generic $5 \times 5$, $6 \times 6$, ..., $n \times n$ matrices
- Matrices whose structure doesn't reduce to degree $\leq 4$

### Galois Theory Background

The Abel-Ruffini theorem is a consequence of Galois theory:

1. **Solvable groups**: Polynomial roots can be expressed in radicals iff its Galois group is solvable
2. **Symmetric groups**: The generic degree-$n$ polynomial has Galois group $S_n$ (symmetric group)
3. **Non-solvability**: $S_n$ is not solvable for $n \geq 5$
4. **Conclusion**: No general radical formula for degree $\geq 5$

**Key insight**: Specific polynomials may have special Galois groups that ARE solvable, even for degree $\geq 5$. This is why structure detection is crucial.

## Why Structure Matters

### The Structure-Exploitation Strategy

Since we can't solve general $n \times n$ matrices for $n \geq 5$, we must:

1. **Detect structure** in the matrix
2. **Exploit structure** to reduce problem complexity
3. **Solve** using available methods

### Types of Exploitable Structure

#### 1. Reducible Structure

**Block-diagonal matrices**:
```math
\mathbf{A} = \begin{bmatrix}
\mathbf{A}_1 & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{A}_2 & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{A}_3
\end{bmatrix}
```
Eigenvalues: $\lambda(\mathbf{A}) = \lambda(\mathbf{A}_1) \cup \lambda(\mathbf{A}_2) \cup \lambda(\mathbf{A}_3)$

**Why this helps**: An $8 \times 8$ block diagonal with $2 \times 2$ blocks requires solving 4 quadratic equations, not one degree-8 equation.

#### 2. Symmetry Structure

**Circulant matrices**: Each row is cyclic shift of previous
```math
\mathbf{C} = \begin{bmatrix}
c_0 & c_1 & c_2 & \cdots & c_{n-1} \\
c_{n-1} & c_0 & c_1 & \cdots & c_{n-2} \\
c_{n-2} & c_{n-1} & c_0 & \cdots & c_{n-3} \\
\vdots & \vdots & \vdots & \ddots & \vdots
\end{bmatrix}
```

Eigenvalues: $\lambda_k = \sum_{j=0}^{n-1} c_j \cdot \omega^{kj}$ where $\omega = \exp(2\pi i/n)$

**Why this helps**: The DFT matrix diagonalizes ALL circulant matrices, giving closed form for any $n$.

#### 3. Tensor Product Structure

**Kronecker products**: $\mathbf{A} \otimes \mathbf{B}$
```math
\mathbf{A} \otimes \mathbf{B} = \begin{bmatrix}
a_{11}\mathbf{B} & a_{12}\mathbf{B} & \cdots \\
a_{21}\mathbf{B} & a_{22}\mathbf{B} & \cdots \\
\vdots & \vdots & \ddots
\end{bmatrix}
```

Eigenvalues: $\lambda(\mathbf{A} \otimes \mathbf{B}) = \{\lambda_i(\mathbf{A}) \cdot \lambda_j(\mathbf{B}) : \text{all } i,j\}$

**Why this helps**: A $6 \times 6 = (2 \times 2) \otimes (3 \times 3)$ reduces to quadratic + cubic, not degree-6.

#### 4. Special Pattern Structure

**Symmetric Toeplitz tridiagonal**: Constant diagonals
```math
\mathbf{T} = \begin{bmatrix}
a & b & 0 & 0 & \cdots & 0 \\
b & a & b & 0 & \cdots & 0 \\
0 & b & a & b & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots
\end{bmatrix}
```

Eigenvalues: $\lambda_k = a + 2b \cos\left(\frac{k\pi}{n+1}\right)$ for $k = 1, 2, \ldots, n$

**Why this helps**: Known eigenvector basis (discrete sine transform) gives closed-form eigenvalues.

### The Key Principle

**Matrices with structure have special Galois groups.**

The characteristic polynomial of a structured matrix isn't generic—it has symmetries that make its Galois group smaller or more solvable, even for degree $\geq 5$.

**Examples**:
- Circulant $n \times n$: Galois group is cyclic $\mathbb{Z}_n$ (always solvable)
- Generic $n \times n$: Galois group is symmetric $S_n$ (not solvable for $n \geq 5$)

## Closed-Form Root Formulas

### Degree 1: Linear

Equation: $\lambda + c_0 = 0$

**Solution**: $\lambda = -c_0$

### Degree 2: Quadratic

Equation: $\lambda^2 + c_1\lambda + c_0 = 0$

**Solution** (Quadratic formula):
```math
\lambda = \frac{-c_1 \pm \sqrt{c_1^2 - 4c_0}}{2}
```

**Discriminant**: $\Delta = c_1^2 - 4c_0$
- $\Delta > 0$: Two real roots
- $\Delta = 0$: One repeated root
- $\Delta < 0$: Two complex conjugate roots

### Degree 3: Cubic (Cardano's Formula)

Equation: $\lambda^3 + c_2\lambda^2 + c_1\lambda + c_0 = 0$

**Strategy**: Eliminate quadratic term via substitution $\lambda = t - c_2/3$

Reduced form: $t^3 + pt + q = 0$

**Solution**:
```math
t = \sqrt[3]{-\frac{q}{2} + \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}} + \sqrt[3]{-\frac{q}{2} - \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}}
```

**Discriminant**: $\Delta = -4p^3 - 27q^2$
- $\Delta > 0$: Three distinct real roots
- $\Delta = 0$: Repeated root
- $\Delta < 0$: One real root, two complex conjugates

**Historical note**: Cardano published this in 1545, but Tartaglia discovered it earlier.

### Degree 4: Quartic (Ferrari's Method)

Equation: $\lambda^4 + c_3\lambda^3 + c_2\lambda^2 + c_1\lambda + c_0 = 0$

**Strategy**: 
1. Eliminate cubic term: $\lambda = t - c_3/4$
2. Add and subtract clever term to create perfect squares
3. Reduce to solving a cubic (Cardano) plus quadratics

**Solution** (simplified outline):
1. Solve cubic resolvent: $y^3 + p_2y^2 + p_1y + p_0 = 0$ (use Cardano)
2. Use root $y$ to decompose into two quadratics
3. Solve quadratics for final four roots

**Complexity**: Formula involves nested cube roots and square roots, producing extremely large expressions.

**Expression size**: For fully symbolic $4 \times 4$ matrix, each eigenvalue is ~13.5 MB of symbolic expressions!

### Degree ≥ 5: No General Formula

**Abel-Ruffini Theorem**: No algebraic solution using radicals.

**Our approach**: Detect structure and apply specialized methods.

## Special Pattern Theory

### Why Certain Patterns Have Closed Forms

A matrix pattern has closed-form eigenvalues if one of these conditions holds:

#### 1. Commutes with Known Transformation

If $\mathbf{A}\mathbf{U} = \mathbf{U}\mathbf{\Lambda}$ where:
- $\mathbf{U}$ is a known unitary matrix (independent of matrix entries)
- $\mathbf{\Lambda}$ is diagonal

Then $\mathbf{U}$ diagonalizes $\mathbf{A}$, and eigenvalues are diagonal entries of $\mathbf{\Lambda} = \mathbf{U}^\dagger\mathbf{A}\mathbf{U}$.

**Examples**:
- Circulant: $\mathbf{U}$ = DFT matrix
- Symmetric Toeplitz tridiagonal: $\mathbf{U}$ = discrete sine transform

#### 2. Reducible via Transformation

If a known transformation splits $\mathbf{A}$ into independent subproblems of degree $\leq 4$.

**Examples**:
- Block-diagonal: Already split
- Persymmetric: Symmetry transformation splits into half-sized problems
- Block circulant: Block DFT reduces to smaller problems

#### 3. Tensor Product Structure

If $\mathbf{A} = \mathbf{B} \otimes \mathbf{C}$, then $\lambda(\mathbf{A}) = \{\lambda_i(\mathbf{B}) \cdot \lambda_j(\mathbf{C})\}$.

**Works when**: Can detect Kronecker product structure.

#### 4. Polynomial in Simple Matrix

If $\mathbf{A} = f(\mathbf{J})$ where $\mathbf{J}$ is a simple matrix (e.g., shift operator, permutation).

**Examples**:
- Circulant: Polynomial in cyclic shift matrix
- Symmetric Toeplitz tridiagonal: Related to Chebyshev polynomials

#### 5. Special Symmetry Group

If the matrix has a symmetry group that:
- Constrains eigenvalue structure
- Reduces to solvable subproblems

**Examples**:
- Anti-diagonal: Persymmetric symmetry forces $\pm$ pairs
- Permutation: Cycle structure gives roots of unity

### The Discovery Process

Finding new patterns requires:

1. **Pattern hypothesis**: Guess a family of matrices
2. **Numerical testing**: Check if eigenvalues have simple form
3. **Symbolic verification**: Verify closed form holds symbolically
4. **Mathematical justification**: Explain WHY it works (find the structure)
5. **Generalization**: Identify what structural property makes it work

See `notes/DISCOVERY_METHODOLOGY.md` for practical exploration techniques.

### Current Limitations

**Unknown patterns**: Many patterns with closed forms remain undiscovered.

**Detection is hard**: Even when pattern has closed form, detecting it automatically is challenging.

**Partial solutions**: Some patterns have partial closed forms (e.g., some eigenvalues closed, others not).

## Group Theory and Symmetry-Based Diagonalization

Group theory provides powerful tools for symbolic eigenvalue computation, especially for matrices with $n \geq 5$ where direct algebraic methods fail. The key insight: **matrix symmetries constrain eigenvalue structure**.

### Representation Theory Fundamentals

#### Basic Principle

When a matrix $\mathbf{A}$ commutes with a group $G$ of symmetry operations:
```math
g \mathbf{A} = \mathbf{A} g \quad \forall g \in G
```

Then:
1. **Eigenvectors organize into irreducible representations** (irreps) of $G$
2. **Eigenvalues within an irrep are equal** (degeneracy)
3. **Character tables predict eigenspace dimensions** before computation

#### Key Groups for Matrix Analysis

**1. Cyclic Groups $\mathbb{Z}_n$** (Already exploited in circulant matrices)

- **Elements**: $\{e, g, g^2, \ldots, g^{n-1}\}$ where $g^n = e$
- **Irreps**: All 1-dimensional, labeled by $k = 0, 1, \ldots, n-1$
- **Characters**: $\chi_k(g^j) = \omega^{jk}$ where $\omega = e^{2\pi i/n}$
- **Application**: Circulant matrices are diagonalized by DFT (Fourier basis = irrep basis)

**2. Dihedral Groups $D_n$** (Rotation + Reflection symmetries)

- **Elements**: $n$ rotations + $n$ reflections (total: $2n$ elements)
- **Order**: $|D_n| = 2n$
- **Irreps**: 
  - Two 1-dimensional (symmetric and antisymmetric under reflection)
  - $(n-2)/2$ two-dimensional irreps (for even $n$)
  - $(n-1)/2$ two-dimensional irreps (for odd $n$)
- **Application**: Matrices with both rotational and mirror symmetry

**Example (Pentagon - $D_5$)**: 
- Matrix with pentagonal symmetry (5-fold rotation + 5 reflection axes)
- Eigenspaces: dim 1 (symmetric), dim 1 (antisymmetric), dim 2, dim 2
- Without solving: know two eigenvalues are simple, two are doubly degenerate

**3. Symmetric Group $S_n$** (All permutations)

- **Elements**: All $n!$ permutations
- **Irreps**: Labeled by Young diagrams
- **Key fact**: Generic $n \times n$ matrix has $S_n$ Galois group
- **Application**: Understanding why generic $n \geq 5$ case is impossible

**4. Point Groups** (3D Molecular symmetries)

Used in quantum chemistry for molecular orbitals:
- $C_n$: Cyclic rotation
- $D_n$: Dihedral (rotation + reflection)
- $T_d$: Tetrahedral
- $O_h$: Octahedral
- $I_h$: Icosahedral

### Galois Groups and Solvability by Radicals

#### Connection to Eigenvalue Problem

For characteristic polynomial $p(\lambda) = \det(\lambda I - A)$:

**Galois group** = Symmetries of roots that preserve polynomial relations

**Key theorem**: Polynomial solvable by radicals $\iff$ Galois group is solvable

#### Solvable Groups (Can express eigenvalues symbolically)

1. **Cyclic groups** $\mathbb{Z}_n$: Always solvable
   - Example: $\lambda^5 - 2 = 0$ has Galois group $\mathbb{Z}_5$
   - Roots: $\sqrt[5]{2} \cdot \omega^k$, $k = 0,1,2,3,4$ where $\omega = e^{2\pi i/5}$

2. **Dihedral groups** $D_n$: Always solvable
   - Example: $\lambda^5 - \lambda - 1 = 0$ often has $D_5$ Galois group
   - Solvable in radicals (though expressions may be complex)

3. **Frobenius groups**: Solvable
   - Semi-direct products of cyclic groups

4. **Nilpotent groups**: Solvable
   - Products of $p$-groups

#### Non-Solvable Groups (Cannot express in radicals)

1. **Symmetric group** $S_n$ for $n \geq 5$: NOT solvable
   - Generic degree-$n$ polynomial
   - Abel-Ruffini theorem consequence

2. **Alternating group** $A_n$ for $n \geq 5$: NOT solvable
   - Even permutations only

#### Practical Application: Galois Group Detection

**Strategy for degree 5-8 polynomials**:

1. **Compute resolvent polynomial**
   - Degree 5: Sextic resolvent
   - Degree 6: Quintic resolvent
   
2. **Check factorization patterns**:
   - Complete factorization → Cyclic group (solvable!)
   - Irreducible of specific degree → Dihedral (solvable!)
   - Stays degree $n!$ → Likely $S_n$ (not solvable)

3. **Use discriminant**:
   - Square discriminant → Galois group in $A_n$
   - Study discriminant factors

**Example (Degree 5)**:

For $\lambda^5 + c_3\lambda^3 + c_2\lambda^2 + c_1\lambda + c_0 = 0$:

- If factors as $(λ^2 + a)(λ^3 + b) = 0$ → Reducible, solvable
- If $c_3 = c_2 = 0$ → Likely cyclic Galois group
- Resolvent factorization reveals Galois group structure

### Commutant-Based Diagonalization

#### Theoretical Foundation

**Schur's Lemma**: If $\mathbf{A}$ commutes with all elements of an irreducible representation, then $\mathbf{A}$ is a scalar multiple of identity in that irrep space.

**Consequence**: 
- Find transformations $\mathbf{T}$ that commute with $\mathbf{A}$
- $\mathbf{T}$'s eigenvectors are also $\mathbf{A}$'s eigenvectors
- If $\mathbf{T}$ is "simple" (known eigenvectors), we get $\mathbf{A}$'s eigenvectors!

#### Applications Already in Code

**1. Circulant Matrices**

- Commute with cyclic shift operator $\mathbf{S}$
- $\mathbf{S}$ is diagonalized by DFT matrix
- Therefore DFT diagonalizes ALL circulant matrices

**2. Symmetric Toeplitz Tridiagonal**

- Commutes with certain reflection operators
- Eigenvectors are discrete sine transform (DST)
- Closed-form eigenvalues: $\lambda_k = a + 2b\cos(k\pi/(n+1))$

#### New Opportunities

**3. Block-Circulant with Circulant Blocks (BCCB)**

Structure:
```math
\mathbf{A} = \begin{bmatrix}
C_0 & C_1 & \cdots & C_{n-1} \\
C_{n-1} & C_0 & \cdots & C_{n-2} \\
\vdots & \vdots & \ddots & \vdots \\
C_1 & C_2 & \cdots & C_0
\end{bmatrix}
```
where each $C_i$ is itself circulant.

- **Symmetry**: Cyclic in both block structure and within blocks
- **Diagonalization**: 2D discrete Fourier transform
- **Eigenvalues**: $\lambda_{jk} = \sum_{p,q} c_{pq} \omega_n^{jp} \omega_m^{kq}$
- **Application**: Image processing, 2D signal processing

**4. Dihedral Symmetry**

Matrix commuting with:
- Rotation: $\mathbf{R}$ (order $n$)
- Reflection: $\mathbf{F}$ (order 2)

Detection criteria:
- $\mathbf{A}\mathbf{R} = \mathbf{R}\mathbf{A}$
- $\mathbf{A}\mathbf{F} = \mathbf{F}\mathbf{A}$
- $\mathbf{R}^n = \mathbf{I}$, $\mathbf{F}^2 = \mathbf{I}$, $\mathbf{F}\mathbf{R}\mathbf{F} = \mathbf{R}^{-1}$

**Consequence**: 
- Eigenvalues come in patterns dictated by irreps
- Some eigenvalues non-degenerate, others doubly degenerate
- Can block-diagonalize into smaller problems

### Graph Symmetries and Spectral Graph Theory

#### Automorphism Groups

For graph adjacency matrix $\mathbf{A}$ or Laplacian $\mathbf{L}$:

**Graph automorphism** = Permutation of vertices preserving edges

**Symmetry → Degeneracy**: 
- If vertex $i$ and $j$ are in same orbit, eigenvector components related
- Automorphism group determines eigenvalue multiplicities

#### Strongly Regular Graphs

**Definition**: Graph with parameters $(n, k, \lambda, \mu)$ where:
- $n$ vertices
- Each vertex has $k$ neighbors (regular)
- Adjacent vertices have $\lambda$ common neighbors
- Non-adjacent vertices have $\mu$ common neighbors

**Remarkable property**: **Only 3 distinct eigenvalues!**

Eigenvalues:
```math
\lambda_1 = k \text{ (multiplicity 1)}
```
```math
\lambda_2 = \frac{(\lambda - \mu) + \sqrt{(\lambda-\mu)^2 + 4(k-\mu)}}{2} \text{ (multiplicity } f \text{)}
```
```math
\lambda_3 = \frac{(\lambda - \mu) - \sqrt{(\lambda-\mu)^2 + 4(k-\mu)}}{2} \text{ (multiplicity } g \text{)}
```

where $f + g = n - 1$.

**Examples**:
- **Petersen graph** (10 vertices): Eigenvalues $\{3, 1^5, -2^4\}$
- **Paley graphs**: $(p, (p-1)/2, (p-5)/4, (p-1)/4)$ for prime $p \equiv 1 \pmod 4$
- **Complete bipartite** $K_{m,n}$: Eigenvalues $\{\pm\sqrt{mn}, 0^{m+n-2}\}$

**Detection**:
1. Check if adjacency matrix is 0-1
2. Verify regularity: all row sums equal
3. Count common neighbors for adjacent and non-adjacent pairs
4. If all pairs follow pattern → strongly regular → instant eigenvalues!

#### Distance-Regular Graphs

Generalization of strongly regular:
- Number of neighbors at distance $d$ depends only on $d$
- Eigenvalues satisfy polynomial recurrence relations
- Many have closed-form eigenvalue expressions

**Examples**:
- **Hypercubes** $Q_n$: Eigenvalues $\{n-2k : k=0,1,\ldots,n\}$ with multiplicity $\binom{n}{k}$
- **Johnson graphs** $J(n,k)$: Eigenvalues related to binomial coefficients
- **Hamming graphs**: Product structure → Kronecker product of simpler graphs

### Lie Groups and Continuous Symmetries

#### Motivation

Some matrices commute with **continuous** symmetry groups:
- Rotations in 2D or 3D
- Quantum mechanical angular momentum
- Special unitary transformations

#### SO(2) - Rotations in 2D

Matrix commuting with all 2D rotations:
```math
\mathbf{R}(\theta) = \begin{bmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{bmatrix}
```

**Consequence**: Must be proportional to identity or to rotation itself

**Application**: Angular momentum, circular symmetry problems

#### SO(3) - Rotations in 3D

Irreps labeled by angular momentum $j = 0, 1/2, 1, 3/2, 2, \ldots$

Dimension: $2j + 1$

**Example**: $j = 1$ (spin-1, dimension 3)
- Three eigenvalues: $-1, 0, +1$ (z-component of angular momentum)

**Application**: Molecular orbital theory, atomic physics

#### SU(2) - Spin Systems

Similar to SO(3) but for half-integer spins.

**Example**: Spin-$5/2$ system (6-dimensional)
- Eigenvalues: $m = -5/2, -3/2, -1/2, +1/2, +3/2, +5/2$
- **Equally spaced** by symmetry alone!
- No polynomial solving needed

**Detection**:
1. Check if matrix commutes with $\mathbf{J}_x, \mathbf{J}_y, \mathbf{J}_z$ (angular momentum operators)
2. Verify $[\mathbf{J}_i, \mathbf{J}_j] = i\epsilon_{ijk}\mathbf{J}_k$
3. If yes → eigenvalues are $\{j, j-1, j-2, \ldots, -j\}$ where $2j+1 = n$

### Practical Detection Strategies

#### Algorithm: Detect Matrix Symmetries

```
Input: Matrix A
Output: Symmetry group and eigenvalue structure

1. Check basic structure:
   - Circulant? → Use DFT
   - Block-diagonal? → Recurse on blocks
   - Triangular? → Read diagonal

2. Test group generators:
   - For each candidate generator g:
     * Compute A*g and g*A
     * If equal, g is a symmetry
   - Identify minimal generating set

3. Determine group:
   - If cyclic generator only → ℤ_n
   - If cyclic + reflection → D_n
   - If all permutations → S_n (generic case)

4. Apply representation theory:
   - Look up character table
   - Decompose into irreps
   - Determine eigenspace dimensions
   - Construct symmetry-adapted basis

5. Solve reduced problems:
   - Each irrep gives smaller eigenvalue problem
   - Often degree ≤ 4 even if original is n ≥ 5
```

#### Algorithm: Galois Group Analysis (Degree 5-8)

```
Input: Characteristic polynomial p(λ)
Output: Solvability and solution strategy

1. Compute discriminant Δ
   - If Δ is a perfect square → Galois group ⊆ A_n

2. Compute resolvent polynomial
   - Degree 5 → sextic resolvent
   - Check if it factors

3. Analyze factorization:
   - Completely factors → Cyclic group → Use roots of unity
   - Factors into degree-2 and degree-3 → Reducible → Solve separately
   - Irreducible → Check for dihedral group

4. If solvable:
   - Construct tower of field extensions
   - Express roots as nested radicals
   
5. If not solvable:
   - Return numerical approximation
   - Warn user about algebraic impossibility
```

### Implementation Roadmap

#### High Priority (Immediately Useful)

1. **Dihedral symmetry detector**
   - Input: Matrix
   - Detect rotation and reflection generators
   - Block-diagonalize by irreps
   - Solve reduced problems (often degree ≤ 2)

2. **Strongly regular graph detector**
   - Input: 0-1 adjacency matrix
   - Verify strong regularity property
   - Compute three eigenvalues from parameters
   - Works for arbitrary size $n$

3. **Galois group analyzer**
   - Input: Degree 5-8 polynomial
   - Compute resolvent and discriminant
   - Identify Galois group
   - If solvable, provide symbolic solution

#### Medium Priority

4. **BCCB matrix handler**
   - Detect block-circulant with circulant blocks
   - Apply 2D DFT
   - Common in image processing

5. **Hypercube graph detector**
   - Identify $n$-dimensional hypercube structure
   - Eigenvalues: $\{n-2k\}$ with multiplicities $\binom{n}{k}$

6. **Lie algebra symmetry checker**
   - Test commutation with SO(2), SO(3), SU(2) generators
   - If symmetric, use representation theory formulas

### Examples in Practice

#### Example 1: Pentagon Symmetry (5×5, Dihedral $D_5$)

Adjacency matrix of regular pentagon:
```math
\mathbf{A} = \begin{bmatrix}
0 & 1 & 0 & 0 & 1 \\
1 & 0 & 1 & 0 & 0 \\
0 & 1 & 0 & 1 & 0 \\
0 & 0 & 1 & 0 & 1 \\
1 & 0 & 0 & 1 & 0
\end{bmatrix}
```

**Without group theory**: Solve degree-5 polynomial (impossible in general)

**With group theory**:
- Has $D_5$ symmetry (5-fold rotation + reflections)
- $D_5$ irreps: two 1D, two 2D
- Eigenvalues: $\lambda_1 = 2$ (1D), $\lambda_2 = -2$ (1D), then two pairs
- Eigenvalues: $2, \phi, -\phi^{-1}, -\phi^{-1}, \phi$ where $\phi = (1+\sqrt{5})/2$ (golden ratio!)

**Key**: Group theory reduces degree-5 problem to quadratic!

#### Example 2: Petersen Graph (10×10, Strongly Regular)

Parameters: $(10, 3, 0, 1)$
- 10 vertices
- Each vertex has 3 neighbors
- Adjacent vertices have 0 common neighbors
- Non-adjacent vertices have 1 common neighbor

**Eigenvalues** (from strong regularity alone):
```math
\lambda_1 = 3 \text{ (multiplicity 1)}
```
```math
\lambda_2 = 1 \text{ (multiplicity 5)}
```
```math
\lambda_3 = -2 \text{ (multiplicity 4)}
```

**No polynomial solving needed!** Pure symmetry argument.

#### Example 3: 6×6 Spin-5/2 System

Matrix representing $J_z$ (angular momentum operator):
```math
\mathbf{J}_z = \text{diag}\left(-\frac{5}{2}, -\frac{3}{2}, -\frac{1}{2}, \frac{1}{2}, \frac{3}{2}, \frac{5}{2}\right)
```

More interesting: $\mathbf{J}_x$ or $\mathbf{J}_y$ matrices (non-diagonal)

**With SU(2) symmetry**:
- Eigenvalues MUST be $\{-5/2, -3/2, -1/2, 1/2, 3/2, 5/2\}$
- Equally spaced by 1
- Order may vary but values are determined by symmetry

**Application**: Quantum mechanics, no polynomial solving needed for any spin-$j$ system.

#### Example 4: Solvable Degree-5 (Cyclic Galois Group)

Polynomial: $\lambda^5 - 3 = 0$

**Galois group**: $\mathbb{Z}_5$ (cyclic)

**Solution**: 
```math
\lambda_k = \sqrt[5]{3} \cdot e^{2\pi ik/5}, \quad k = 0,1,2,3,4
```

**Explicitly**:
```math
\lambda_k = \sqrt[5]{3} \cdot \left(\cos\frac{2\pi k}{5} + i\sin\frac{2\pi k}{5}\right)
```

**Key**: Cyclic Galois group → solvable by radicals, even though degree = 5.

## Further Reading

### Classic References

1. **Galois Theory**: 
   - Stewart, I. (2015). *Galois Theory* (4th ed.). Chapman and Hall/CRC.
   
2. **Matrix Theory**:
   - Horn, R. A., & Johnson, C. R. (2012). *Matrix Analysis* (2nd ed.). Cambridge University Press.

3. **Special Matrices**:
   - Davis, P. J. (1979). *Circulant Matrices* (2nd ed.). AMS Chelsea Publishing.

4. **Group Representation Theory**:
   - Serre, J.-P. (1977). *Linear Representations of Finite Groups*. Springer.
   - Fulton, W., & Harris, J. (1991). *Representation Theory: A First Course*. Springer.

5. **Spectral Graph Theory**:
   - Godsil, C., & Royle, G. (2001). *Algebraic Graph Theory*. Springer.
   - Brouwer, A. E., & Haemers, W. H. (2012). *Spectra of Graphs*. Springer.

### Online Resources

- [Galois Theory Explained](https://en.wikipedia.org/wiki/Galois_theory)
- [Abel-Ruffini Theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem)
- [Circulant Matrices](https://en.wikipedia.org/wiki/Circulant_matrix)
- [Kronecker Product](https://en.wikipedia.org/wiki/Kronecker_product)
- [Representation Theory](https://en.wikipedia.org/wiki/Representation_theory)
- [Dihedral Group](https://en.wikipedia.org/wiki/Dihedral_group)
- [Strongly Regular Graphs](https://en.wikipedia.org/wiki/Strongly_regular_graph)
- [Character Tables](https://en.wikipedia.org/wiki/Character_table)

### Research Papers

For specific patterns implemented in this package, see references in `notes/PATTERN_DISCOVERIES.md`.
