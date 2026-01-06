# Mathematical Background

This page covers the mathematical theory underlying SymbolicDiagonalization.jl.

## The Eigenvalue Problem

Given a matrix $A \in \mathbb{C}^{n \times n}$, the eigenvalue problem seeks values $\lambda$ and vectors $v$ satisfying:

$$Av = \lambda v$$

The eigenvalues are roots of the **characteristic polynomial**:

$$p(\lambda) = \det(A - \lambda I)$$

## The Abel-Ruffini Barrier

The Abel-Ruffini theorem (1824) states that polynomials of degree 5 or higher have no general closed-form solution in radicals. This means:

- **2×2, 3×3, 4×4**: Always solvable via quadratic, Cardano, and Ferrari formulas
- **5×5+**: No general formula exists

However, matrices with special structure can bypass this barrier.

## Closed-Form Root Formulas

### Quadratic (Degree 2)

For $x^2 + bx + c = 0$:

$$x = \frac{-b \pm \sqrt{b^2 - 4c}}{2}$$

### Cardano's Formula (Degree 3)

For $x^3 + px + q = 0$ (depressed cubic):

$$x = \sqrt[3]{-\frac{q}{2} + \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}} + \sqrt[3]{-\frac{q}{2} - \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}}$$

The three roots are obtained by multiplying cube roots by $1, \omega, \omega^2$ where $\omega = e^{2\pi i/3}$.

### Ferrari's Formula (Degree 4)

The quartic $x^4 + ax^3 + bx^2 + cx + d = 0$ is solved by:

1. Substitute $x = y - a/4$ to eliminate the cubic term
2. Find an auxiliary cubic (the resolvent)
3. Factor the quartic into two quadratics
4. Solve each quadratic

## Structure Exploitation

### Block-Diagonal Decomposition

If $A$ is block-diagonal:

$$A = \begin{pmatrix} B_1 & 0 \\ 0 & B_2 \end{pmatrix}$$

Then $\text{eigenvalues}(A) = \text{eigenvalues}(B_1) \cup \text{eigenvalues}(B_2)$.

This reduces an $n \times n$ problem to smaller subproblems.

### Persymmetric Matrices

A matrix $A$ is **persymmetric** if $JAJ = A$ where $J$ is the exchange (anti-diagonal) matrix. Such matrices decompose:

$$P^T A P = \begin{pmatrix} A_+ & 0 \\ 0 & A_- \end{pmatrix}$$

reducing an $n \times n$ problem to two $n/2 \times n/2$ problems.

## Circulant Matrices and the DFT

An $n \times n$ circulant matrix with first row $(c_0, c_1, \ldots, c_{n-1})$ is diagonalized by the DFT matrix $F$:

$$A = F^{-1} \Lambda F$$

where $\Lambda = \text{diag}(\lambda_0, \ldots, \lambda_{n-1})$ with:

$$\lambda_k = \sum_{j=0}^{n-1} c_j \omega^{jk}, \quad \omega = e^{2\pi i/n}$$

**Key insight**: The eigenvectors (columns of $F$) are independent of the specific entries—only the eigenvalues depend on $c_j$.

### Roots of Unity

For circulant matrices, eigenvalues involve exact roots of unity:

| $n$ | $\omega = e^{2\pi i/n}$ | Algebraic form |
|-----|-------------------------|----------------|
| 2 | $-1$ | $-1$ |
| 3 | $e^{2\pi i/3}$ | $-1/2 + i\sqrt{3}/2$ |
| 4 | $i$ | $i$ |
| 6 | $e^{\pi i/3}$ | $1/2 + i\sqrt{3}/2$ |

## Kronecker Products

For matrices $A \in \mathbb{C}^{m \times m}$ and $B \in \mathbb{C}^{n \times n}$:

$$\text{eigenvalues}(A \otimes B) = \{\lambda_i(A) \cdot \mu_j(B) : i=1,\ldots,m; j=1,\ldots,n\}$$

**Proof**: If $Av = \lambda v$ and $Bw = \mu w$, then $(A \otimes B)(v \otimes w) = \lambda\mu (v \otimes w)$.

### Nested Kronecker Products

For $A_1 \otimes A_2 \otimes \cdots \otimes A_k$ where each $A_i$ is $2 \times 2$:
- Size: $2^k \times 2^k$
- Eigenvalues: products $\lambda_{j_1}(A_1) \cdot \lambda_{j_2}(A_2) \cdots \lambda_{j_k}(A_k)$

A $1024 \times 1024$ matrix ($k=10$) reduces to 10 quadratic problems!

## Lie Group Eigenvalues

### SO(2): 2D Rotations

$$R(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

Eigenvalues: $e^{\pm i\theta}$

**Proof**: Characteristic polynomial is $\lambda^2 - 2\cos\theta \cdot \lambda + 1 = 0$.

### SO(3): 3D Rotations

Any $R \in SO(3)$ has eigenvalues $\{1, e^{i\theta}, e^{-i\theta}\}$ where $\theta$ is the rotation angle.

**Proof**: $\det(R) = 1$, so eigenvalues multiply to 1. Complex eigenvalues come in conjugate pairs, so one eigenvalue must be real (±1). Since $R$ is a rotation (not reflection), this eigenvalue is 1.

### SU(2): Special Unitary Group

$$U = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix}, \quad |\alpha|^2 + |\beta|^2 = 1$$

Eigenvalues: $e^{\pm i\theta/2}$ where $\cos(\theta/2) = \text{Re}(\alpha)$.

## Symmetric Toeplitz Tridiagonal

A symmetric Toeplitz tridiagonal matrix with diagonal $a$ and off-diagonal $b$:

$$T = \begin{pmatrix} a & b & & \\ b & a & b & \\ & \ddots & \ddots & \ddots \\ & & b & a \end{pmatrix}$$

has eigenvalues:

$$\lambda_k = a + 2b\cos\left(\frac{k\pi}{n+1}\right), \quad k = 1, \ldots, n$$

**Proof**: Eigenvectors are $v_k = (\sin(k\pi/(n+1)), \sin(2k\pi/(n+1)), \ldots)$, verified by Chebyshev polynomial recurrence.

## Hadamard Matrices

The Sylvester-Hadamard matrix $H_n$ of order $2^n$ satisfies:

$$H_n = \begin{pmatrix} H_{n-1} & H_{n-1} \\ H_{n-1} & -H_{n-1} \end{pmatrix}$$

**Eigenvalues**: $\pm\sqrt{2^n}$, each with multiplicity $2^{n-1}$.

**Proof**: $H_n H_n^T = 2^n I$, so eigenvalues satisfy $\lambda^2 = 2^n$.

## DFT Eigenvalues

The $n \times n$ DFT matrix $F$ with entries $F_{jk} = \omega^{jk}$ satisfies $F^4 = nI$.

**Eigenvalues**: Fourth roots of $n$, i.e., $\{\pm\sqrt{n}, \pm i\sqrt{n}\}$.

For perfect square $n = m^2$, eigenvalues are exact integers: $\{\pm m, \pm mi\}$.

## Cartan Matrices (Type Aₙ)

The Cartan matrix of type $A_n$ (associated with $SU(n+1)$):

$$C = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 2 \end{pmatrix}$$

has eigenvalues:

$$\lambda_k = 4\sin^2\left(\frac{k\pi}{2(n+1)}\right), \quad k = 1, \ldots, n$$

**Proof**: This is a symmetric Toeplitz tridiagonal with $a=2$, $b=-1$.

## Special Angle Simplification

Trigonometric values at special angles have algebraic forms:

| Angle | $\cos$ | $\sin$ |
|-------|--------|--------|
| $\pi/6$ | $\sqrt{3}/2$ | $1/2$ |
| $\pi/5$ | $(1+\sqrt{5})/4$ | $\sqrt{10-2\sqrt{5}}/4$ |
| $\pi/4$ | $\sqrt{2}/2$ | $\sqrt{2}/2$ |
| $\pi/3$ | $1/2$ | $\sqrt{3}/2$ |
| $2\pi/5$ | $(\sqrt{5}-1)/4$ | $\sqrt{10+2\sqrt{5}}/4$ |

The values at $\pi/5$ and $2\pi/5$ involve the **golden ratio** $\phi = (1+\sqrt{5})/2$.

## Complexity Analysis

| Pattern | Matrix Size | Polynomial Degree | Effective Degree |
|---------|-------------|-------------------|------------------|
| Block-diagonal | $n \times n$ | $n$ | $\max(n_1, n_2, \ldots)$ |
| Circulant | $n \times n$ | $n$ | 1 (closed form) |
| Kronecker $A \otimes B$ | $mn \times mn$ | $mn$ | $\max(m, n)$ |
| Nested $2^{\otimes k}$ | $2^k \times 2^k$ | $2^k$ | 2 |
| Hadamard | $2^n \times 2^n$ | $2^n$ | 1 (two values) |

## References

1. **Abel-Ruffini Theorem**: Ruffini (1799), Abel (1824)
2. **Circulant Matrices**: Davis, P. J. (1979). *Circulant Matrices*
3. **Lie Groups**: Hall, B. C. (2015). *Lie Groups, Lie Algebras, and Representations*
4. **Root Formulas**: Cardano (1545), Ferrari (1540)
5. **Hadamard Matrices**: Sylvester (1867)
