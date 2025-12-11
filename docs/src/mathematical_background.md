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

## Further Reading

### Classic References

1. **Galois Theory**: 
   - Stewart, I. (2015). *Galois Theory* (4th ed.). Chapman and Hall/CRC.
   
2. **Matrix Theory**:
   - Horn, R. A., & Johnson, C. R. (2012). *Matrix Analysis* (2nd ed.). Cambridge University Press.

3. **Special Matrices**:
   - Davis, P. J. (1979). *Circulant Matrices* (2nd ed.). AMS Chelsea Publishing.

### Online Resources

- [Galois Theory Explained](https://en.wikipedia.org/wiki/Galois_theory)
- [Abel-Ruffini Theorem](https://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem)
- [Circulant Matrices](https://en.wikipedia.org/wiki/Circulant_matrix)
- [Kronecker Product](https://en.wikipedia.org/wiki/Kronecker_product)

### Research Papers

For specific patterns implemented in this package, see references in `notes/PATTERN_DISCOVERIES.md`.
