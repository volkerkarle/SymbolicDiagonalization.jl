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
```
D = [d₁  0   0  ...  0 ]
    [0   d₂  0  ...  0 ]
    [0   0   d₃ ...  0 ]
    [...            ...]
    [0   0   0  ... dₙ]
```

**Eigenvalues**: `{d₁, d₂, d₃, ..., dₙ}` (diagonal entries)

**Eigenvectors**: Standard basis vectors eᵢ

**Complexity**: O(n) - direct reading

**Example**:
```julia
@variables a b c
D = [a 0 0; 0 b 0; 0 0 c]

eigvals(D)  # [a, b, c]
```

**Why this works**: Diagonal matrices are already in eigenvalue form. The eigenvalue equation Dv = λv is satisfied by standard basis vectors.

**Mathematical properties**:
- Always diagonalizable
- Eigenvectors are orthogonal
- Simple eigenvalue structure

---

### Triangular Matrices

**Structure**: Upper or lower triangular (zeros below or above diagonal)

**Upper Triangular Form**:
```
U = [u₁₁ u₁₂ u₁₃ ... u₁ₙ]
    [0   u₂₂ u₂₃ ... u₂ₙ]
    [0   0   u₃₃ ... u₃ₙ]
    [...               ...]
    [0   0   0   ... uₙₙ]
```

**Eigenvalues**: `{u₁₁, u₂₂, u₃₃, ..., uₙₙ}` (diagonal entries)

**Eigenvectors**: Computed via back-substitution

**Complexity**: O(n) for eigenvalues, O(n²) for eigenvectors

**Example**:
```julia
@variables a b c
U = [a 1 2; 0 b 3; 0 0 c]

eigvals(U)  # [a, b, c]
```

**Why this works**: For triangular matrices, det(λI - U) is a product of diagonal terms (λ - uᵢᵢ), so eigenvalues are exactly the diagonal entries.

**Mathematical properties**:
- May not be diagonalizable (if repeated eigenvalues)
- Eigenvalues independent of off-diagonal entries
- Characteristic polynomial factorizes immediately

---

### Jordan Blocks

**Structure**: Single repeated eigenvalue on diagonal, ones on superdiagonal

**Matrix Form**:
```
J(λ, n) = [λ  1  0  ...  0]
          [0  λ  1  ...  0]
          [0  0  λ  ...  0]
          [...           ...]
          [0  0  0  ...  λ]
```

**Eigenvalues**: `{λ}` with algebraic multiplicity n

**Eigenvectors**: Only ONE eigenvector (not diagonalizable for n > 1)

**Complexity**: O(1) for eigenvalues

**Example**:
```julia
@variables λ
J = [λ 1 0; 0 λ 1; 0 0 λ]

eigvals(J)  # [λ, λ, λ]
```

**Why this works**: det(tI - J) = (t - λ)ⁿ, giving repeated root λ. Geometric multiplicity is 1 (rank deficiency).

**Mathematical properties**:
- **Not diagonalizable** (unless n = 1)
- Prototype of non-diagonalizable matrices
- Appears in Jordan normal form decomposition

**Caution**: `symbolic_diagonalize()` will throw error since Jordan blocks (n > 1) are not diagonalizable.

---

## Structure-Based Patterns

### Block-Diagonal Matrices

**Structure**: Non-zero only in diagonal blocks

**Matrix Form**:
```
B = [A₁  0   0  ...  0 ]
    [0   A₂  0  ...  0 ]
    [0   0   A₃ ...  0 ]
    [...             ...]
    [0   0   0  ... Aₖ]
```
where each Aᵢ is a square matrix (possibly different sizes).

**Eigenvalues**: λ(B) = λ(A₁) ∪ λ(A₂) ∪ ... ∪ λ(Aₖ)

**Complexity**: O(Σᵢ f(nᵢ)) where f(n) is cost to solve size-n block

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
```
det(λI - B) = det(λI - A₁) · det(λI - A₂) · ... · det(λI - Aₖ)
```

So eigenvalues are union of eigenvalues of each block.

**Detection algorithm**:
1. Find connected components of non-zero structure
2. Recursively solve each block independently
3. Concatenate eigenvalues

**Practical importance**: An 8×8 matrix with four 2×2 blocks requires solving 4 quadratics, not one degree-8 polynomial!

**Mathematical properties**:
- Diagonalizable iff all blocks are diagonalizable
- Eigenvectors are block-structured
- Extremely common in practice (physical subsystems, symmetry blocks)

---

### Persymmetric Matrices

**Structure**: Symmetric about anti-diagonal: Q[i,j] = Q[n+1-j, n+1-i]

**Matrix Form** (4×4 example):
```
P = [a  b  c  d]
    [b  e  f  c]
    [c  f  e  b]
    [d  c  b  a]
```

**Eigenvalues**: Computed by splitting into half-sized problems

**Complexity**: O(f(n/2)) where f is cost for size n/2

**Example**:
```julia
@variables a b c d e f
P = [a b c d; b e f c; c f e b; d c b a]

eigvals(P)  # Splits into two 2×2 problems
```

**Why this works**: Persymmetric structure allows transformation:
```
Q = J * P * J  (where J is anti-diagonal identity)
```

This creates block structure in eigenvector space, splitting problem in half.

**Transformation details**:
- Let J be the flip matrix (anti-diagonal identity)
- Form symmetric combinations: (P ± JPJ)/2
- These commute and can be simultaneously diagonalized
- Reduces to two half-sized eigenvalue problems

**Detection**: Check if Q[i,j] = Q[n+1-j, n+1-i] for all i,j

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
```
C = [c₀   c₁   c₂   ... cₙ₋₁]
    [cₙ₋₁  c₀   c₁   ... cₙ₋₂]
    [cₙ₋₂  cₙ₋₁  c₀   ... cₙ₋₃]
    [...                   ...]
    [c₁   c₂   c₃   ... c₀  ]
```

**Eigenvalues** (closed-form for any n):
```
λₖ = Σⱼ₌₀ⁿ⁻¹ cⱼ · ωᵏʲ  for k = 0, 1, ..., n-1
```
where ω = exp(2πi/n) is the nth root of unity.

**Eigenvectors**: Columns of DFT matrix (independent of entries!)

**Complexity**: O(n) to compute all n eigenvalues

**Example (3×3)**:
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
1. Check if C[i,j] = C[(i-1) mod n, (j-1) mod n] for all i,j
2. Extract first row [c₀, c₁, ..., cₙ₋₁]
3. Apply DFT formula

**Applications**:
- Signal processing (circular convolution)
- Time series analysis (periodic processes)
- Graph Laplacians (cycle graphs)
- Coding theory

**Practical notes**:
- Works for ANY size n (even n = 100+!)
- Generic n×n requires degree-n polynomial (impossible for n ≥ 5)
- Circulant structure bypasses Abel-Ruffini limitation

**Variants**:
- **Skew-circulant**: Last element negated on each shift
- **g-circulant**: Generalized with group operation
- **Block circulant**: See next section

---

### Block Circulant Matrices

**Structure**: Block version of circulant (blocks shift instead of elements)

**Matrix Form** (4 blocks, size k×k each):
```
BC = [A₀  A₁  A₂  A₃]
     [A₃  A₀  A₁  A₂]
     [A₂  A₃  A₀  A₁]
     [A₁  A₂  A₃  A₀]
```

**Eigenvalues**: Compute by solving n eigenvalue problems of size k×k:
```
λ(BC) = ⋃ₘ₌₀ⁿ⁻¹ λ(Σⱼ₌₀ⁿ⁻¹ Aⱼ · ωᵐʲ)
```
where ω = exp(2πi/n).

**Complexity**: O(n · f(k)) where f(k) is cost for k×k matrix

**Example (4×4 with 2×2 blocks)**:
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

**Why this works**: Block circulant matrices are diagonalized by block DFT. The reduction formula creates n different k×k matrices (linear combinations of blocks with DFT coefficients), which can be solved independently.

**Mathematical foundation**:
- Generalization of circulant to matrix entries
- Block DFT plays same role as regular DFT
- Eigenvalue problem decouples into n independent k×k problems

**Detection algorithm**:
1. Check block structure (equal-sized blocks)
2. Verify block circulant property
3. Extract blocks [A₀, A₁, ..., Aₙ₋₁]
4. Form n combinations with DFT weights
5. Solve each k×k system

**Key examples**:
- **12×12 with 3×3 blocks**: 4 blocks → solve 4 cubic equations (feasible)
- **8×8 with 4×4 blocks**: 2 blocks → solve 2 quartic equations (feasible but large)
- **16×16 with 2×2 blocks**: 8 blocks → solve 8 quadratic equations (easy)

Without block structure, these would require degree-12, degree-8, and degree-16 polynomials (impossible)!

**Applications**:
- Multi-channel signal processing
- Block Toeplitz systems
- Vectorized circulant operations
- Kronecker-structured problems

**Limitation**: Requires k ≤ 4 for closed form (unless further structure in blocks).

---

### Anti-Diagonal Matrices

**Structure**: Non-zero only on anti-diagonal, with symmetry

**Matrix Form** (n×n, symmetric about anti-diagonal):
```
A = [0    0    0   ... 0   aₙ  ]
    [0    0    0   ... aₙ₋₁ 0  ]
    [0    0    0   ... 0    0  ]
    [...                       ]
    [aₙ₋₁ 0    0   ... 0    0  ]
    [aₙ   0    0   ... 0    0  ]
```
with symmetry: A[i,j] = A[n+1-j, n+1-i]

**Eigenvalues** (closed-form for any n):
- **Odd n**: `{aₖ, -aₖ}` for k = 1, ..., (n-1)/2, plus `a₍ₙ₊₁₎/₂` (center element)
- **Even n**: `{aₖ, -aₖ}` for k = 1, ..., n/2

**Eigenvector structure**: ±pairs with specific symmetry

**Complexity**: O(n) - direct computation

**Example (5×5)**:
```julia
@variables a b c
A = [0 0 0 0 a;
     0 0 0 b 0;
     0 0 c 0 0;
     0 b 0 0 0;
     a 0 0 0 0]

eigvals(A)  # {c, ±a, ±b}
```

**Why this works**: The persymmetric symmetry combined with anti-diagonal structure forces eigenvalues into ±pairs. The transformation J*A (flip matrix) equals ±A, creating symmetric/antisymmetric subspaces.

**Mathematical foundation**:
- Anti-diagonal flip symmetry: JAJ† = ±A
- Eigenvectors come in symmetric/antisymmetric pairs
- Eigenvalue equation: If Av = λv, then A(Jv) = -λ(Jv)

**Detection algorithm**:
1. Check if zeros everywhere except anti-diagonal
2. Verify symmetry: A[i, n+1-i] = A[n+1-i, i]
3. Extract anti-diagonal elements
4. Form ±pairs (center unpaired if odd n)

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

**Structure**: Tensor product A ⊗ B = [aᵢⱼ·B]

**Matrix Form** (2×2 ⊗ 2×2 example):
```
A ⊗ B = [a₁₁B  a₁₂B]
        [a₂₁B  a₂₂B]
```

where each entry aᵢⱼ·B is a block.

**Eigenvalues** (closed-form for any sizes m, n):
```
λ(A ⊗ B) = {λᵢ(A) · λⱼ(B) : i = 1..m, j = 1..n}
```

All m·n products of eigenvalues of A and B.

**Eigenvectors**: vᵢⱼ = vᵢ(A) ⊗ vⱼ(B) (tensor products)

**Complexity**: O(f(m) + f(n)) where f is cost for each factor

**Example (4×4 = 2×2 ⊗ 2×2)**:
```julia
@variables a b c d
A = [a 0; 0 b]  # eigenvalues {a, b}
B = [c 0; 0 d]  # eigenvalues {c, d}

K = kron(A, B)  # 4×4 matrix

eigvals(K)  # {ac, ad, bc, bd}
```

**Why this works**: The eigenvalue equation for Kronecker products:
```
(A ⊗ B)(u ⊗ v) = (Au) ⊗ (Bv) = (λᵤ·u) ⊗ (λᵥ·v) = (λᵤ·λᵥ)(u ⊗ v)
```

So if u is eigenvector of A with eigenvalue λᵤ, and v is eigenvector of B with eigenvalue λᵥ, then u⊗v is eigenvector of A⊗B with eigenvalue λᵤ·λᵥ.

**Mathematical foundation**:
- Tensor product structure in linear algebra
- Factorization of eigenspaces: E(A⊗B) = E(A) ⊗ E(B)
- Multiplicative property of eigenvalues

**Detection algorithm**:
1. Check block pattern: equal-sized blocks
2. Verify scaling relationship between blocks
3. Extract factor matrices A and B
4. Solve eigenvalues of A and B independently
5. Form all products λᵢ(A) · λⱼ(B)

**Key examples**:
- **6×6 = 2×2 ⊗ 3×3**: Quadratic × cubic = feasible
- **12×12 = 3×3 ⊗ 4×4**: Cubic × quartic = feasible (but large)
- **12×12 = 4×4 ⊗ 3×3**: Same as above

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
```
T = [a  b  0  0  ... 0]
    [b  a  b  0  ... 0]
    [0  b  a  b  ... 0]
    [...            ...]
    [0  0  0  b  ... a]
```

**Eigenvalues** (closed-form for any n):
```
λₖ = a + 2b·cos(kπ/(n+1))  for k = 1, 2, ..., n
```

**Eigenvectors** (known analytically):
```
vₖ(j) = sin(jkπ/(n+1))  for j = 1, ..., n
```

**Complexity**: O(n) to compute all eigenvalues

**Example (4×4)**:
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
- Comes from finite-difference discretization of d²/dx²

**Derivation**:
The eigenvalue equation Tv = λv with v = [sin(kπ/(n+1)), sin(2kπ/(n+1)), ..., sin(nkπ/(n+1))] gives:
```
a·sin(jkπ/(n+1)) + b·sin((j-1)kπ/(n+1)) + b·sin((j+1)kπ/(n+1)) = λₖ·sin(jkπ/(n+1))
```

Using trigonometric identity sin(α)+sin(β) = 2sin((α+β)/2)cos((α-β)/2), this simplifies to λₖ = a + 2b·cos(kπ/(n+1)).

**Detection algorithm**:
1. Check tridiagonal structure (zeros beyond sub/super-diagonals)
2. Verify constant diagonals: all diagonal = a, all off-diagonal = b
3. Verify symmetry: sub-diagonal = super-diagonal
4. Apply formula for k = 1, ..., n

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

**Why it scales to any n**: The eigenvector basis is KNOWN and independent of matrix size. We don't need to solve polynomial—we directly compute eigenvalues from formula.

---

### Special 5×5 Tridiagonal Patterns

**Structure**: Tridiagonal 5×5 with one "perturbation" from constant pattern

Two specific patterns have been discovered with closed-form eigenvalues.

#### Pattern 1: [b, d, b, b]

**Matrix Form**:
```
M = [a  b  0  0  0]
    [b  a  d  0  0]
    [0  d  a  b  0]
    [0  0  b  a  b]
    [0  0  0  b  a]
```

Off-diagonal sequence: [b, d, b, b]

**Eigenvalues** (closed-form):
```
{a + √(2b² + d²), a - √(2b² + d²), a + b, a - b, a}
```

**Example**:
```julia
@variables a b d
M = [a b 0 0 0; b a d 0 0; 0 d a b 0; 0 0 b a b; 0 0 0 b a]

eigvals(M)  # {a ± √(2b² + d²), a ± b, a}
```

#### Pattern 2: [b, b, d, b]

**Matrix Form**:
```
M = [a  b  0  0  0]
    [b  a  b  0  0]
    [0  b  a  d  0]
    [0  0  d  a  b]
    [0  0  0  b  a]
```

Off-diagonal sequence: [b, b, d, b]

**Eigenvalues** (closed-form - SAME as Pattern 1!):
```
{a + √(2b² + d²), a - √(2b² + d²), a + b, a - b, a}
```

**Example**:
```julia
@variables a b d
M = [a b 0 0 0; b a b 0 0; 0 b a d 0; 0 0 d a b; 0 0 0 b a]

eigvals(M)  # {a ± √(2b² + d²), a ± b, a}  (identical to Pattern 1!)
```

**Remarkable discovery**: Despite different positions of perturbation d, both patterns yield IDENTICAL eigenvalues!

**Why this works**: Unknown! This is an empirical discovery. The mathematical explanation likely involves:
- Hidden symmetry in characteristic polynomial
- Isospectrality (different matrices, same spectrum)
- Possibly related to Chebyshev polynomial structure

**What doesn't work** (patterns without closed form):
- **[d, b, b, b]**: Boundary perturbation at position 0
- **[b, b, b, d]**: Boundary perturbation at position 3
- **[d, d, b, b]**, **[b, d, d, b]**: Multiple perturbations
- **[c, d, e, f]**: General perturbations

**Key insight from research**: Interior perturbations (positions 1-2) admit closed forms, but boundary perturbations (positions 0, 3) break the pattern.

**Detection algorithm**:
1. Check 5×5 tridiagonal structure
2. Extract super-diagonal sequence
3. Match against known patterns [b,d,b,b] or [b,b,d,b]
4. Verify constant diagonal = a
5. Apply closed-form formula

**Applications**:
- Perturbed vibrating strings
- Defects in 1D lattice models
- Specific quantum systems

**Research opportunity**: 
- Why do these specific patterns work?
- Are there similar 7×7 or 9×9 patterns?
- Can we characterize all solvable perturbation patterns?

See `notes/PATTERN_DISCOVERIES.md` for more details on discovery methodology.

---

## Permutation Patterns

### Permutation Matrices

**Structure**: Exactly one 1 in each row and column, rest zeros

**Matrix Form** (example: permutation 1→2, 2→3, 3→1):
```
P = [0 1 0]
    [0 0 1]
    [1 0 0]
```

**Eigenvalues** (closed-form for any permutation):
- Each k-cycle contributes k-th roots of unity: {1, ω, ω², ..., ωᵏ⁻¹} where ω = exp(2πi/k)
- All eigenvalues have magnitude 1 (on unit circle)

**Complexity**: O(n) to compute cycle decomposition and roots

**Example (6×6 with 3-cycle, 2-cycle, fixed point)**:
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
1. Permutation matrices have finite order: Pᵏ = I for some k
2. If Pᵏv = v, then P has eigenvalues that are k-th roots of unity
3. Cycle decomposition determines which roots appear

**Mathematical foundation**:
- **Cayley's theorem**: Every permutation decomposes into disjoint cycles
- **Root of unity property**: If Pᵏ = I, then λᵏ = 1, so λ is a k-th root of unity
- **Cycle length = root order**: k-cycle gives k-th roots of unity

**Cycle decomposition algorithm**:
1. Find permutation: which column has 1 in each row
2. Follow cycle: start at 1, go to P(1), then P(P(1)), until back to 1
3. Repeat for unvisited elements
4. Count cycle lengths

**Eigenvalue computation**:
For each k-cycle, add k eigenvalues {exp(2πij/k) : j = 0, ..., k-1}

**Applications**:
- Sorting algorithms (permutation analysis)
- Group theory (symmetric group representations)
- Graph automorphisms
- Coding theory (permutation codes)
- Fourier analysis on symmetric group

**Mathematical properties**:
- Always orthogonal: P†P = I
- Eigenvalues on unit circle: |λ| = 1
- Diagonalizable (orthogonally)
- Determinant = ±1 (sign of permutation)

**Special cases**:
- **Identity**: One n-cycle → eigenvalue 1 with multiplicity n
- **Full cycle**: One n-cycle → all n-th roots of unity {exp(2πik/n) : k=0..n-1}
- **Transposition**: (2-cycle) → eigenvalues {1, -1}
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
| Diagonal | Any n | O(n) | Direct reading |
| Triangular | Any n | O(n) | Diagonal entries |
| Jordan Block | Any n | O(1) | Single eigenvalue |
| Block-Diagonal | Any n | Σ O(nᵢ) | Union of blocks |
| Persymmetric | Any n | O(n/2) | Half-size split |
| Circulant | Any n | O(n) | DFT formula |
| Block Circulant | Any n (k ≤ 4) | O(n·k) | Block DFT |
| Kronecker A⊗B | Any m,n (≤ 4) | O(m+n) | Product rule |
| Toeplitz Tridiag | Any n | O(n) | Cosine formula |
| Anti-Diagonal | Any n | O(n) | ±pairs |
| Permutation | Any n | O(n) | Roots of unity |
| Special 5×5 [b,d,b,b] | 5×5 only | O(1) | Empirical |
| Special 5×5 [b,b,d,b] | 5×5 only | O(1) | Empirical |

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
