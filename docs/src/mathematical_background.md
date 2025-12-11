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

Given an n×n matrix **A**, we seek scalars λ (eigenvalues) and non-zero vectors **v** (eigenvectors) satisfying:

```
A v = λ v
```

Equivalently, we need:
```
(A - λI) v = 0
```

For non-trivial solutions, the matrix (A - λI) must be singular:
```
det(A - λI) = 0
```

This determinant is a polynomial in λ called the **characteristic polynomial**:
```
p(λ) = det(λI - A) = λⁿ + cₙ₋₁λⁿ⁻¹ + ... + c₁λ + c₀
```

**The eigenvalue problem reduces to finding roots of a polynomial.**

### Why Eigenvalues Matter

Eigenvalues reveal fundamental properties of linear transformations:

- **Dynamics**: Matrix powers Aⁿ grow/decay based on eigenvalue magnitudes
- **Stability**: System stability determined by eigenvalue real parts
- **Geometry**: Eigenvalues measure stretching/compression along eigenvector directions
- **Spectral decomposition**: A = PDP⁻¹ where D is diagonal (if A is diagonalizable)

### The Fundamental Challenge

For general matrices:
- 2×2 → degree-2 polynomial (quadratic formula)
- 3×3 → degree-3 polynomial (cubic formula)
- 4×4 → degree-4 polynomial (quartic formula)
- 5×5 → degree-5 polynomial (**no general formula!**)
- n×n → degree-n polynomial (**no general formula for n ≥ 5**)

This limitation is not computational—it's a fundamental impossibility proven by the Abel-Ruffini theorem.

## The Abel-Ruffini Theorem

### Statement

**Theorem** (Abel-Ruffini, 1824): There is no general algebraic solution (using radicals) for polynomial equations of degree 5 or higher.

### What This Means

- **Degree ≤ 4**: Closed-form formulas exist using +, -, ×, ÷, and nth roots
- **Degree ≥ 5**: No general formula exists using these operations
- **Specific polynomials**: May still have closed-form solutions (e.g., λ⁵ - 1 = 0)

### Implications for Symbolic Eigenvalues

✅ **Can solve symbolically**:
- All 1×1, 2×2, 3×3, 4×4 general matrices
- Larger matrices with exploitable structure

❌ **Cannot solve symbolically** (in general):
- Generic 5×5, 6×6, ..., n×n matrices
- Matrices whose structure doesn't reduce to degree ≤ 4

### Galois Theory Background

The Abel-Ruffini theorem is a consequence of Galois theory:

1. **Solvable groups**: Polynomial roots can be expressed in radicals iff its Galois group is solvable
2. **Symmetric groups**: The generic degree-n polynomial has Galois group Sₙ (symmetric group)
3. **Non-solvability**: Sₙ is not solvable for n ≥ 5
4. **Conclusion**: No general radical formula for degree ≥ 5

**Key insight**: Specific polynomials may have special Galois groups that ARE solvable, even for degree ≥ 5. This is why structure detection is crucial.

## Why Structure Matters

### The Structure-Exploitation Strategy

Since we can't solve general n×n matrices for n ≥ 5, we must:

1. **Detect structure** in the matrix
2. **Exploit structure** to reduce problem complexity
3. **Solve** using available methods

### Types of Exploitable Structure

#### 1. Reducible Structure

**Block-diagonal matrices**:
```
A = [A₁  0   0 ]
    [0   A₂  0 ]
    [0   0   A₃]
```
Eigenvalues: λ(A) = λ(A₁) ∪ λ(A₂) ∪ λ(A₃)

**Why this helps**: An 8×8 block diagonal with 2×2 blocks requires solving 4 quadratic equations, not one degree-8 equation.

#### 2. Symmetry Structure

**Circulant matrices**: Each row is cyclic shift of previous
```
C = [c₀  c₁  c₂  ...  cₙ₋₁]
    [cₙ₋₁ c₀  c₁  ...  cₙ₋₂]
    [cₙ₋₂ cₙ₋₁ c₀  ...  cₙ₋₃]
    ...
```

Eigenvalues: λₖ = Σⱼ cⱼ · ωᵏʲ where ω = exp(2πi/n)

**Why this helps**: The DFT matrix diagonalizes ALL circulant matrices, giving closed form for any n.

#### 3. Tensor Product Structure

**Kronecker products**: A ⊗ B
```
A ⊗ B = [a₁₁B  a₁₂B  ...
         a₂₁B  a₂₂B  ...
         ...]
```

Eigenvalues: λ(A ⊗ B) = {λᵢ(A) · λⱼ(B) : all i,j}

**Why this helps**: A 6×6 = (2×2) ⊗ (3×3) reduces to quadratic + cubic, not degree-6.

#### 4. Special Pattern Structure

**Symmetric Toeplitz tridiagonal**: Constant diagonals
```
T = [a  b  0  0  ...
     b  a  b  0  ...
     0  b  a  b  ...
     ...]
```

Eigenvalues: λₖ = a + 2b cos(kπ/(n+1)) for k = 1, 2, ..., n

**Why this helps**: Known eigenvector basis (discrete sine transform) gives closed-form eigenvalues.

### The Key Principle

**Matrices with structure have special Galois groups.**

The characteristic polynomial of a structured matrix isn't generic—it has symmetries that make its Galois group smaller or more solvable, even for degree ≥ 5.

**Examples**:
- Circulant n×n: Galois group is cyclic ℤₙ (always solvable)
- Generic n×n: Galois group is symmetric Sₙ (not solvable for n ≥ 5)

## Closed-Form Root Formulas

### Degree 1: Linear

Equation: λ + c₀ = 0

**Solution**: λ = -c₀

### Degree 2: Quadratic

Equation: λ² + c₁λ + c₀ = 0

**Solution** (Quadratic formula):
```
λ = (-c₁ ± √(c₁² - 4c₀)) / 2
```

**Discriminant**: Δ = c₁² - 4c₀
- Δ > 0: Two real roots
- Δ = 0: One repeated root
- Δ < 0: Two complex conjugate roots

### Degree 3: Cubic (Cardano's Formula)

Equation: λ³ + c₂λ² + c₁λ + c₀ = 0

**Strategy**: Eliminate quadratic term via substitution λ = t - c₂/3

Reduced form: t³ + pt + q = 0

**Solution**:
```
t = ∛(-q/2 + √(q²/4 + p³/27)) + ∛(-q/2 - √(q²/4 + p³/27))
```

**Discriminant**: Δ = -4p³ - 27q²
- Δ > 0: Three distinct real roots
- Δ = 0: Repeated root
- Δ < 0: One real root, two complex conjugates

**Historical note**: Cardano published this in 1545, but Tartaglia discovered it earlier.

### Degree 4: Quartic (Ferrari's Method)

Equation: λ⁴ + c₃λ³ + c₂λ² + c₁λ + c₀ = 0

**Strategy**: 
1. Eliminate cubic term: λ = t - c₃/4
2. Add and subtract clever term to create perfect squares
3. Reduce to solving a cubic (Cardano) plus quadratics

**Solution** (simplified outline):
1. Solve cubic resolvent: y³ + p₂y² + p₁y + p₀ = 0 (use Cardano)
2. Use root y to decompose into two quadratics
3. Solve quadratics for final four roots

**Complexity**: Formula involves nested cube roots and square roots, producing extremely large expressions.

**Expression size**: For fully symbolic 4×4 matrix, each eigenvalue is ~13.5 MB of symbolic expressions!

### Degree ≥ 5: No General Formula

**Abel-Ruffini Theorem**: No algebraic solution using radicals.

**Our approach**: Detect structure and apply specialized methods.

## Special Pattern Theory

### Why Certain Patterns Have Closed Forms

A matrix pattern has closed-form eigenvalues if one of these conditions holds:

#### 1. Commutes with Known Transformation

If **AU = UΛ** where:
- **U** is a known unitary matrix (independent of matrix entries)
- **Λ** is diagonal

Then U diagonalizes A, and eigenvalues are diagonal entries of Λ = U†AU.

**Examples**:
- Circulant: U = DFT matrix
- Symmetric Toeplitz tridiagonal: U = discrete sine transform

#### 2. Reducible via Transformation

If a known transformation splits A into independent subproblems of degree ≤ 4.

**Examples**:
- Block-diagonal: Already split
- Persymmetric: Symmetry transformation splits into half-sized problems
- Block circulant: Block DFT reduces to smaller problems

#### 3. Tensor Product Structure

If A = B ⊗ C, then λ(A) = {λᵢ(B) · λⱼ(C)}.

**Works when**: Can detect Kronecker product structure.

#### 4. Polynomial in Simple Matrix

If A = f(J) where J is a simple matrix (e.g., shift operator, permutation).

**Examples**:
- Circulant: Polynomial in cyclic shift matrix
- Symmetric Toeplitz tridiagonal: Related to Chebyshev polynomials

#### 5. Special Symmetry Group

If the matrix has a symmetry group that:
- Constrains eigenvalue structure
- Reduces to solvable subproblems

**Examples**:
- Anti-diagonal: Persymmetric symmetry forces ±pairs
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
