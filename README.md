# SymbolicDiagonalization.jl

**Status: Functional Prototype with 9 Special Patterns**

A Julia package for symbolic matrix diagonalization using closed-form root solvers and structure detection.

## The Challenge

The Abel-Ruffini theorem proves that no general closed-form solution exists for polynomials of degree ≥5. This means we can't solve general 5×5+ matrices symbolically. However, many real-world matrices have exploitable structure (block-diagonal, circulant, Toeplitz, etc.) that allows closed-form eigenvalue computation regardless of size.

## The Solution

**SymbolicDiagonalization.jl** automatically detects and exploits matrix structure to solve larger symbolic problems:

- **Closed-form root solvers** for degrees 1-4 (linear, quadratic, cubic, quartic)
- **Automatic structure detection** (block-diagonal, persymmetric, Hermitian)
- **9 special pattern solvers** that work for any matrix size n (circulant, block circulant, Kronecker products, symmetric Toeplitz tridiagonal, anti-diagonal, permutation, special 5×5 patterns)
- **Clean LinearAlgebra.jl interface** (eigen(), eigvals())

## Current State

**What works well**:
- ✅ All matrices up to 4×4 (via quartic formula)
- ✅ Block-diagonal matrices (recursive decomposition)
- ✅ 9 special patterns for arbitrary-sized matrices (3×3 to n×n)
- ✅ 172 passing tests, comprehensive test suite

**Limitations**:
- ❌ General 5×5+ matrices (requires structure detection)
- ❌ Expression explosion for fully symbolic 4×4 (~13.5 MB per eigenvalue)
- ⚠️  Basic structure detection (misses many patterns)
- ⚠️  Minimal expression simplification

## Quick Start

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## What's Implemented

### 1. Closed-Form Root Solvers (degrees 1-4)

Solves characteristic polynomials up to degree 4:
- **Linear** (1×1 matrices): direct solution
- **Quadratic** (2×2): standard quadratic formula
- **Cubic** (3×3): Cardano's method
- **Quartic** (4×4): Ferrari's method (produces large expressions)

The characteristic polynomial is computed via **Bareiss fraction-free determinant** on (λI - A).

### 2. Structure Detection

Automatically detects and exploits these matrix structures:
- **Diagonal/Triangular**: Direct eigenvalue extraction O(n)
- **Block-diagonal**: Recursive decomposition into independent subproblems
- **Persymmetric**: Splits matrices with persymmetric symmetry (Q[i,j] = Q[n+1-j,n+1-i])
- **Hermitian**: Optimization paths for Hermitian matrices

### 3. Special Pattern Solvers (any size n ≥ 3)

Nine specialized solvers that work for arbitrary-sized matrices:

1. **Circulant matrices**: DFT-based closed-form eigenvalues (O(n) complexity)
2. **Block circulant matrices**: Block DFT reduction to smaller eigenvalue problems
3. **Kronecker products** (A ⊗ B): Eigenvalues as products λᵢ(A) · λⱼ(B)
4. **Symmetric Toeplitz tridiagonal**: Closed-form via trigonometric formula
5. **Anti-diagonal matrices**: Eigenvalues in ±pairs
6. **Permutation matrices**: Roots of unity from cycle decomposition
7. **Special 5×5 pattern [b,d,b,b]**: Closed-form eigenvalues
8. **Special 5×5 pattern [b,b,d,b]**: Closed-form eigenvalues (same as [b,d,b,b])
9. **Jordan blocks**: Repeated eigenvalue on diagonal

## Installation

From package root:

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Examples

All examples below are tested and working with the current implementation.

### Basic Usage

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]  # Upper triangular

E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster) → [a, b, c]
```

### Block-Diagonal Matrices

Block-diagonal matrices are automatically detected and decomposed into independent subproblems:

```julia
@variables a b c d
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat)  # [a+b, a-b, c+d, c-d]
```

**How it works**: Detects zero blocks, recursively solves each diagonal block independently, then combines results.

### Circulant Matrices (any size n)

Circulant matrices have each row as a cyclic shift of the previous row. Eigenvalues are given by the DFT of the first row:

```julia
@variables a b c
# 3×3 circulant: each row is cyclic shift of previous
C = [a b c;
     c a b;
     b c a]

vals = eigvals(C)  # Uses ω = exp(2πi/3): [a+b+c, a+b*ω+c*ω², a+b*ω²+c*ω]
```

**Why this works**: Circulant matrices are diagonalized by the DFT matrix, giving eigenvalues λₖ = Σⱼ cⱼ·ωᵏʲ where ω = exp(2πi/n) and c is the first row.

**Performance**: Works for any size n, even n=100. A general 100×100 symbolic matrix would be impossible (degree-100 polynomial), but circulant structure reduces it to O(n) complexity.

### Block Circulant Matrices (any size n)

Block circulant matrices generalize circulant matrices to blocks. An n-block matrix with k×k blocks reduces to n eigenvalue problems of size k×k:

```julia
@variables a b c d
# 4×4 block circulant with 2×2 blocks
A = [a b; c d]
B = [1 0; 0 1]

# Create 4×4 block circulant: [A B; B A]
M = [a b 1 0;
     c d 0 1;
     1 0 a b;
     0 1 c d]

vals = eigvals(M)  # Reduces to 2 problems of size 2×2
```

**Why this works**: Block circulant matrices are diagonalized by block DFT, reducing an (n·k) × (n·k) problem to n problems of size k×k.

**Example**: 
- 12×12 with 3×3 blocks → 4 problems of size 3×3 (solvable via cubic)
- 8×8 with 4×4 blocks → 2 problems of size 4×4 (solvable via quartic)

Without this structure, these would require degree-12 and degree-8 polynomials (no closed form).

### Kronecker Products (any size m×n)

If A is m×m with eigenvalues {λ₁, ..., λₘ} and B is n×n with eigenvalues {μ₁, ..., μₙ}, then A⊗B has eigenvalues {λᵢ·μⱼ : all i,j}:

```julia
@variables a b c d
A = [a 0; 0 b]  # 2×2 diagonal, eigenvalues {a, b}
B = [c 0; 0 d]  # 2×2 diagonal, eigenvalues {c, d}

M = kron(A, B)  # 4×4 Kronecker product

vals = eigvals(M)  # {ac, ad, bc, bd}
```

**Why this works**: The eigenvectors of A⊗B are Kronecker products of eigenvectors of A and B, and eigenvalues multiply.

**Example**: 6×6 = (2×2) ⊗ (3×3) would normally require degree-6 polynomial (no closed form), but Kronecker structure reduces it to quadratic and cubic problems.

### Symmetric Toeplitz Tridiagonal (any size n)

Symmetric tridiagonal matrices with constant diagonals have closed-form eigenvalues via trigonometric formula:

```julia
@variables a b
# 4×4 symmetric Toeplitz tridiagonal
T = [a b 0 0;
     b a b 0;
     0 b a b;
     0 0 b a]

vals = eigvals(T)  # λₖ = a + 2b·cos(kπ/(n+1)) for k=1,2,3,4
```

**Why this works**: These matrices have known eigenvectors (discrete sine transform basis), leading to closed-form eigenvalues.

**Example**: Works for any size n, even n=100. The general n×n case would require degree-n polynomial, but this pattern gives immediate closed form.

### Anti-Diagonal Matrices (any size n)

Symmetric anti-diagonal matrices (non-zero only on the anti-diagonal) have eigenvalues in ±pairs:

```julia
@variables a b c
# 5×5 symmetric anti-diagonal
A = [0 0 0 0 a;
     0 0 0 b 0;
     0 0 c 0 0;
     0 b 0 0 0;
     a 0 0 0 0]

vals = eigvals(A)  # {c, ±a, ±b}  (center element c is unpaired for odd n)
```

**Why this works**: Anti-diagonal symmetry Q[i,j] = Q[n+1-j,n+1-i] combined with persymmetric structure forces eigenvalues into ±pairs.

**Note**: For even n, all eigenvalues come in ±pairs. For odd n, the center element is unpaired.

### Permutation Matrices (any size n)

Permutation matrices have eigenvalues that are roots of unity, determined by cycle decomposition:

```julia
# 6×6 permutation: (1→2→3→1) is a 3-cycle, (4↔5) is a 2-cycle, (6) is fixed
P = [0 1 0 0 0 0;  # 1→2
     0 0 1 0 0 0;  # 2→3
     1 0 0 0 0 0;  # 3→1
     0 0 0 0 1 0;  # 4→5
     0 0 0 1 0 0;  # 5→4
     0 0 0 0 0 1]  # 6→6

vals = eigvals(P)  
# 3-cycle contributes: {1, ω, ω²} where ω = exp(2πi/3)
# 2-cycle contributes: {1, -1}
# Fixed point contributes: {1}
# Total: {1, 1, 1, -1, ω, ω²}
```

**Why this works**: Each k-cycle contributes the k-th roots of unity as eigenvalues. Permutation matrices are diagonalizable with eigenvalues on the unit circle.

**Example**: Any permutation matrix, regardless of size, has eigenvalues that are roots of unity (magnitude 1).

### Special 5×5 Tridiagonal Patterns

Two specific 5×5 tridiagonal patterns have been discovered with closed-form eigenvalues:

**Pattern 1: [b, d, b, b]** (off-diagonals reading left-to-right on super-diagonal)
```julia
@variables a b d
mat = [a  b  0  0  0;
       b  a  d  0  0;
       0  d  a  b  0;
       0  0  b  a  b;
       0  0  0  b  a]

vals = eigvals(mat)  # {a ± √(2b² + d²), a ± b, a}
```

**Pattern 2: [b, b, d, b]**
```julia
@variables a b d
mat = [a  b  0  0  0;
       b  a  b  0  0;
       0  b  a  d  0;
       0  0  d  a  b;
       0  0  0  b  a]

vals = eigvals(mat)  # {a ± √(2b² + d²), a ± b, a}  (same as Pattern 1!)
```

**Remarkable insight**: Both patterns produce **identical eigenvalues** despite different structures. This suggests a deeper symmetry principle at work.

**Limitation**: These are the only known 5×5 tridiagonal patterns with closed forms. Other perturbations like [d, b, b, b] or [b, b, b, d] do not have closed-form solutions.

## Performance Characteristics

### Complexity by Matrix Type

| Matrix Type | Size | Polynomial Degree | Solution Method | Complexity |
|-------------|------|-------------------|-----------------|------------|
| General symbolic | ≤ 4×4 | 1-4 | Root formulas | O(1) |
| General symbolic | ≥ 5×5 | ≥ 5 | ❌ No closed form | Impossible |
| Block-diagonal | n×n | Depends on blocks | Recursive | O(k₁ + k₂ + ...) |
| Circulant | n×n | n | DFT formula | O(n) |
| Block circulant | nk×nk | n blocks of k×k | Block DFT | O(n · k) |
| Kronecker A⊗B | mn×mn | m·n | Product rule | O(m + n) |
| Toeplitz tridiag | n×n | n | Trig formula | O(n) |
| Anti-diagonal | n×n | n | ±pair rule | O(n) |
| Permutation | n×n | n | Cycle decomp | O(n) |

### Expression Size

⚠️ **Warning**: Fully symbolic 4×4 matrices produce extremely large expressions:
- Each eigenvalue: ~13.5 MB of symbolic expressions
- Cause: Ferrari's quartic formula involves nested radicals
- Recommendation: Use numerical methods or add structure

**Example of explosion**:
```julia
@variables a b c d e f g h i j k l m n o p
M = [a b c d; e f g h; i j k l; m n o p]  # Generic 4×4

# This will work but produce huge expressions:
vals = eigvals(M)  # Each of 4 eigenvalues is ~13.5 MB
```

**Better approach**: Add structure or substitute numeric values:
```julia
# Option 1: Block structure
M_block = [a b 0 0; b a 0 0; 0 0 c d; 0 0 d c]  # Manageable

# Option 2: Partial substitution
M_partial = substitute(M, Dict(e=>0, g=>0, i=>0, j=>0, k=>0, m=>0, n=>0, o=>0))
```

## API Reference

### LinearAlgebra Interface (Recommended)

The package extends `LinearAlgebra.jl` with symbolic support:

**`eigen(A; kwargs...)`**
- Returns `Eigen` object with `.values` and `.vectors` fields
- Computes both eigenvalues and eigenvectors
- Example: `E = eigen(mat); E.values, E.vectors`

**`eigvals(A; kwargs...)`** 
- Returns eigenvalues only (faster, skips eigenvector computation)
- Recommended when you only need eigenvalues
- Example: `λ = eigvals(mat)`

### Direct API (Advanced)

For more control over the computation:

**`symbolic_eigenvalues(A; kwargs...)`**
- Returns `(values, poly, λ)` tuple
- `values`: Vector of eigenvalue expressions
- `poly`: Characteristic polynomial
- `λ`: Symbolic variable used in polynomial

**`symbolic_eigenpairs(A; kwargs...)`**
- Returns vector of `(eigenvalue, eigenvectors)` tuples
- Each eigenvalue may have multiple eigenvectors (geometric multiplicity)
- Example: `[(λ₁, [v₁]), (λ₂, [v₂, v₃]), ...]`

**`symbolic_diagonalize(A; kwargs...)`**
- Returns `(P, D, pairs)` where `A = P * D * inv(P)`
- `P`: Matrix of eigenvectors (columns)
- `D`: Diagonal matrix of eigenvalues
- `pairs`: Same as `symbolic_eigenpairs` output
- Throws error if matrix is not diagonalizable

**`characteristic_polynomial(A; var=nothing)`**
- Returns `(poly, coeffs, λ)` tuple
- `poly`: Characteristic polynomial det(λI - A)
- `coeffs`: Polynomial coefficients
- `λ`: Symbolic variable (auto-generated if not provided)

### Keyword Arguments

All functions accept these optional keyword arguments:

**`structure::Symbol`** (default: `:auto`)
- Matrix structure hint for optimization
- Options: `:auto`, `:hermitian`, `:symmetric`, `:unitary`, `:general`
- Example: `eigvals(H, structure=:hermitian)`

**`expand::Bool`** (default: `true`)
- Whether to expand polynomial expressions
- Set to `false` for factored form (if available)
- Example: `eigvals(mat, expand=false)`

**`complexity_threshold::Int`** (default: `5`)
- Warn if symbolic variable count exceeds this threshold
- Helps prevent accidentally running huge symbolic computations
- Set to higher value to suppress warnings

**`timeout::Int`** (default: `300`)
- Maximum computation time in seconds
- Throws `ComputationTimeoutError` if exceeded
- Example: `eigvals(mat, timeout=60)`

**`max_terms::Int`** (default: `10000`)
- Maximum number of terms in symbolic expressions
- Throws `ExpressionComplexityError` if exceeded
- Prevents runaway expression growth

## Current Limitations & Known Issues

### What Doesn't Work

❌ **General 5×5+ matrices**: No closed-form solution exists (Abel-Ruffini theorem)
- Requires structure detection or numeric fallback
- Will throw error if no structure detected

❌ **Expression explosion on symbolic 4×4**: Ferrari's quartic produces ~13.5 MB per eigenvalue
- Use block structure, partial substitution, or numerical methods
- Consider `expand=false` option

❌ **Limited structure detection**: Only detects obvious patterns
- Misses: Hankel, general Toeplitz, rank-1 updates, arrow matrices
- May fall back to degree-n polynomial when closed form exists

⚠️ **Minimal simplification**: Results often not in simplest form
- `a + b - b` might not simplify to `a`
- May need manual simplification with `Symbolics.simplify()`

⚠️ **Eigenvector computation**: Can be slow for complex symbolic expressions
- Use `eigvals()` instead of `eigen()` when only eigenvalues needed
- Eigenvectors computed via nullspace (RREF-based)

⚠️ **No multiplicity handling**: Repeated eigenvalues may not report geometric multiplicity correctly
- Works for simple cases (diagonal, Jordan blocks)
- Complex cases may need manual verification

### Error Messages

**`ComputationTimeoutError`**: Computation exceeded `timeout` seconds
- Try: smaller matrix, add structure, increase timeout, or use numeric methods

**`ExpressionComplexityError`**: Expression exceeded `max_terms` limit
- Try: simplify structure, partial substitution, or increase `max_terms`

**`ArgumentError: Cannot solve degree n ≥ 5 without structure`**: No structure detected
- Try: verify matrix has structure, add structure hint, or use numeric methods

## Roadmap & Future Work

### High Priority (Next Steps)

- [ ] **More special patterns**
  - [ ] Hankel matrices (anti-diagonal symmetry)
  - [ ] General Toeplitz families beyond tridiagonal
  - [ ] Arrow matrices (tridiagonal + rank-1 update)
  - [ ] Companion matrices (standard basis)
  - [ ] Rank-1 updates to diagonal (Sherman-Morrison)
  
- [ ] **Improved structure detection**
  - [ ] Automatic circulant detection
  - [ ] Block structure inference
  - [ ] Pattern matching against known families
  - [ ] Symmetry group analysis

- [ ] **Expression simplification**
  - [ ] Automatic radical simplification
  - [ ] Trigonometric simplification for Toeplitz patterns
  - [ ] Common subexpression elimination
  - [ ] Integration with Symbolics.jl simplification

- [ ] **Comprehensive testing**
  - [ ] Edge cases (zero eigenvalues, repeated roots)
  - [ ] Numerical stability verification
  - [ ] Large-scale pattern testing (n=100+)
  - [ ] Defective matrix handling

### Medium Priority

- [ ] **Better eigenvector computation**
  - [ ] Symbolic normalization options
  - [ ] Orthogonality verification for Hermitian cases
  - [ ] Basis selection heuristics

- [ ] **Eigenvalue multiplicity**
  - [ ] Automatic geometric multiplicity detection
  - [ ] Jordan normal form for defective matrices
  - [ ] Repeated root handling in root solvers

- [ ] **Numerical integration**
  - [ ] Automatic fallback to numerical methods
  - [ ] Mixed symbolic-numeric computation
  - [ ] Condition number estimation
  - [ ] Perturbation analysis

- [ ] **Documentation**
  - [ ] Complete API reference
  - [ ] Mathematical background guide
  - [ ] Pattern library catalog
  - [ ] Tutorial notebooks

### Future Ideas

- [ ] **Machine learning for pattern recognition**
  - [ ] Train classifier on known solvable patterns
  - [ ] Suggest similar patterns from database
  - [ ] Automated pattern discovery

- [ ] **User-defined pattern libraries**
  - [ ] Plugin system for custom patterns
  - [ ] Pattern DSL for specification
  - [ ] Community pattern repository

- [ ] **Advanced features**
  - [ ] Symbolic perturbation theory
  - [ ] Parallel processing for block systems
  - [ ] Lazy evaluation for large expressions
  - [ ] Integration with other CAS systems

## Testing

The package has a comprehensive test suite with **172 tests, all passing** (execution time ~38 seconds):

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

### Test Coverage

- **test/test_basic.jl** (97 lines): Diagonal, triangular, Jordan blocks, quartic, simple 2×2
- **test/test_structure.jl** (119 lines): Hermitian, persymmetric, block-diagonal detection
- **test/test_patterns.jl** (624 lines): All 9 special patterns (108 tests)
  - Circulant (3×3, 4×4, 5×5)
  - Block circulant (4×4, 6×6, 8×8)
  - Kronecker products (4×4, 6×6)
  - Symmetric Toeplitz tridiagonal (3×3 to 6×6)
  - Anti-diagonal (3×3, 4×4, 5×5)
  - Permutation (various cycle structures)
  - Special 5×5 patterns [b,d,b,b] and [b,b,d,b]
- **test/test_interface.jl** (93 lines): LinearAlgebra interface (eigen, eigvals)

All test matrices verified numerically for correctness.

## Pattern Discovery & Documentation

### Research Notes

Comprehensive documentation of discovered patterns and methodology:

- **[docs/notes/RESEARCH_SUMMARY.md](docs/notes/RESEARCH_SUMMARY.md)** - Summary of all pattern discoveries
- **[docs/notes/PATTERN_DISCOVERIES.md](docs/notes/PATTERN_DISCOVERIES.md)** - Detailed pattern catalog with mathematical details
- **[docs/notes/DISCOVERY_METHODOLOGY.md](docs/notes/DISCOVERY_METHODOLOGY.md)** - How to discover new patterns

### Exploration Tools

**[examples/explore_patterns.jl](examples/explore_patterns.jl)** provides interactive tools to discover new solvable patterns:

```julia
include("examples/explore_patterns.jl")

# Systematically test perturbations of known patterns
explore_tridiagonal_perturbations(n=5, base_pattern=[b,b,b,b])

# Test if a specific pattern has closed form
test_pattern(your_matrix)

# Analyze eigenvalue structure
analyze_eigenvalue_structure(mat)
```

**Key insights from research**:
- 5×5 patterns [b,d,b,b] and [b,b,d,b] produce **identical** eigenvalues
- Symmetric perturbations in 7×7 give **partial** closed-forms
- Boundary perturbations (position 0 or n-1) **break** closed-form structure
- Interior perturbations (positions 1-2) work best for small matrices

## Contributing

Contributions are very welcome! This is experimental research software with many opportunities for improvement.

### Ways to Contribute

**1. Implement known special patterns**
- Hankel matrices, arrow matrices, companion matrices
- See linear algebra literature for solvable families
- Add pattern detector + solver + tests

**2. Improve structure detection**
- Better algorithms for recognizing patterns
- Symmetry group analysis
- Sparsity pattern matching

**3. Find new solvable patterns**
- Use `examples/explore_patterns.jl` to discover patterns
- Document discoveries in `docs/notes/`
- Submit patterns with mathematical justification

**4. Performance optimization**
- Expression simplification improvements
- Parallel processing for block matrices
- Caching and memoization

**5. Testing and documentation**
- Edge case testing
- Tutorial notebooks
- API documentation
- Mathematical background guides

### Development Setup

```bash
# Clone repository
git clone https://github.com/yourusername/SymbolicDiagonalization.jl.git
cd SymbolicDiagonalization.jl

# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Build documentation
julia --project=docs docs/make.jl
```

### Code Style

- Follow Julia style guide
- Add docstrings for public functions
- Include tests for new features
- Document mathematical foundations

### Submitting Patterns

When submitting a new pattern, please include:
1. **Mathematical justification**: Why does this pattern have closed-form eigenvalues?
2. **Pattern detector**: Function to detect when a matrix matches the pattern
3. **Solver implementation**: Function to compute eigenvalues
4. **Test cases**: At least 3 test matrices with verified eigenvalues
5. **Documentation**: Add to pattern library with examples

## License

See LICENSE file for details.
