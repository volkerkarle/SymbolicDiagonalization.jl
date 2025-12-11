# SymbolicDiagonalization.jl

**Status: Early Work in Progress**

A Julia package for symbolic matrix diagonalization using closed-form root solvers.

The ultimate goal here is to build a practical symbolic eigenvalue solver that goes beyond simple 3x3 matrices. Closed-form solutions only exist for polynomials up to degree 4 (Abel-Ruffini theorem); this means we can't solve general 5×5+ matrices symbolically. However, many real-world matrices have exploitable structure. If we can automatically detect and exploit these structures, we can potentially solve much larger symbolic problems.

**What we're building**:

- Automatic structure detection (block-diagonal, persymmetric, special patterns)
- Library of special solvable patterns (tridiagonal, circulant, Toeplitz, etc.)
- Tools to discover new solvable patterns
- Clean integration with Symbolics.jl ecosystem

**Current state**: Early but functional. We have block-diagonal detection, circulant matrix support, symmetric Toeplitz tridiagonal support, and some special patterns. Many ideas from the vision are still unimplemented.

## Quick Start

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## What Works

### Closed-Form Root Solvers (1-4 variables)

- Linear, quadratic, cubic (Cardano), quartic (Ferrari)
- Characteristic polynomial via Bareiss determinant
- Works but produces huge expressions for fully symbolic 4×4

### Basic Structure Detection (Work in Progress)

- **Block-diagonal**: Detects and splits into independent subproblems
- **Persymmetric**: Splits certain symmetric-about-antidiagonal matrices
- **Circulant matrices** (any size n): Uses DFT-based closed-form eigenvalue formula
- **Block circulant matrices** (any size n): Reduces to smaller eigenvalue problems using block DFT
- **Kronecker products** (A ⊗ B): Computes eigenvalues as products of factor eigenvalues
- **Symmetric Toeplitz tridiagonal** (any size n): Closed-form eigenvalues using trigonometric formula
- **Symmetric anti-diagonal** (any size n): Eigenvalues come in ±pairs
- **Permutation matrices** (any size n): Eigenvalues are roots of unity based on cycle structure
- **Special 5×5 tridiagonal**: Specific patterns `[b,d,b,b]` with known eigenvalues

## Installation

From package root:

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Examples

### Block-Diagonal (Works Well)

```julia
@variables a b c d
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat)  # [a+b, a-b, c+d, c-d]
```

### Circulant Matrices (Any Size)

Circulant matrices have eigenvalues given by the DFT of the first row. Works for any size n, even beyond degree 4.

```julia
@variables a b c
# 3×3 circulant: each row is cyclic shift of previous
C = [a b c;
     c a b;
     b c a]

vals = eigvals(C)  # Uses ω = exp(2πi/3): [a+b+c, a+b*ω+c*ω², a+b*ω²+c*ω]
```

### Block Circulant Matrices (Any Size)

Block circulant matrices generalize circulant matrices to blocks. An n-block matrix with block size k reduces to n eigenvalue problems of size k×k.

```julia
# 8×8 block circulant with 2×2 blocks (reduces to 4 problems of size 2×2)
A = [1 2; 3 4]
B = [5 6; 7 8]
C = [9 10; 11 12]
D = [13 14; 15 16]

# Each row of blocks is cyclic shift: [A B C D; D A B C; C D A B; B C D A]
M = [A B C D;
     D A B C;
     C D A B;
     B C D A]

vals = eigvals(M)  # Computes eigenvalues by solving 4 smaller 2×2 problems

# Example: 12×12 with 3×3 blocks → 4 problems of size 3×3
# Example: 8×8 with 4×4 blocks → 2 problems of size 4×4 (solvable via quartic)
```

**Performance**: A 12×12 general matrix would require degree-12 polynomial (no closed form), but block circulant structure reduces it to solvable subproblems.

### Kronecker Products (Any Size)

Kronecker products (A ⊗ B) have eigenvalues that are products of the eigenvalues of the factors. If A has eigenvalues {λ₁, ..., λₘ} and B has eigenvalues {μ₁, ..., μₙ}, then A ⊗ B has eigenvalues {λᵢμⱼ : i=1..m, j=1..n}.

```julia
# 4×4 = 2×2 ⊗ 2×2
A = [1 2; 3 4]
B = [5 6; 7 8]
M = kron(A, B)  # Creates the 4×4 Kronecker product

vals = eigvals(M)  # Eigenvalues are all products λᵢ(A) * λⱼ(B)

# Works with symbolic matrices too
@variables a b c d
A_sym = [a 0; 0 b]  # diagonal, eigenvalues {a, b}
B_sym = [c 0; 0 d]  # diagonal, eigenvalues {c, d}
M_sym = kron(A_sym, B_sym)  # 4×4 matrix

vals_sym = eigvals(M_sym)  # {ac, ad, bc, bd}

# Also works with complex numeric matrices
A_complex = ComplexF64[1.0+0.0im 2.0-1.0im; 2.0+1.0im -0.5+0.0im]
B_complex = ComplexF64[-1.5+0.0im 1.0+1.0im; 1.0-1.0im 1.0+0.0im]
M_complex = kron(A_complex, B_complex)

vals_complex = eigvals(M_complex)  # Complex eigenvalues computed efficiently
```

**Performance**: A 6×6 = 2×2 ⊗ 3×3 matrix would normally require degree-6 polynomial (no closed form), but Kronecker structure reduces it to quadratic and cubic problems.

### Symmetric Toeplitz Tridiagonal (Any Size)

Symmetric tridiagonal matrices with constant diagonals have closed-form eigenvalues using cosines.

```julia
@variables a b
# 4×4 symmetric Toeplitz tridiagonal
T = [a b 0 0;
     b a b 0;
     0 b a b;
     0 0 b a]

vals = eigvals(T)  # λₖ = a + 2b*cos(kπ/5) for k=1,2,3,4
```

### Anti-Diagonal Matrices (Any Size)

Symmetric anti-diagonal matrices (non-zero only on the anti-diagonal) have eigenvalues that come in ±pairs.

```julia
@variables a b c
# 5×5 symmetric anti-diagonal
A = [0 0 0 0 a;
     0 0 0 b 0;
     0 0 c 0 0;
     0 b 0 0 0;
     a 0 0 0 0]

vals = eigvals(A)  # {c, ±a, ±b} for odd dimension
```

### Permutation Matrices (Any Size)

Permutation matrices have eigenvalues that are roots of unity, determined by their cycle structure.

```julia
# 6×6 permutation: 1→2→3→1 (3-cycle), 4↔5 (2-cycle), 6 fixed
P = [0 0 1 0 0 0;  # 1→3
     1 0 0 0 0 0;  # 2→1
     0 1 0 0 0 0;  # 3→2
     0 0 0 0 1 0;  # 4→5
     0 0 0 1 0 0;  # 5→4
     0 0 0 0 0 1]  # 6→6

vals = eigvals(P)  # 3rd roots of unity from 3-cycle: {1, ω, ω²}
                   # 2nd roots from 2-cycle: {1, -1}
                   # Fixed point: {1}
                   # Total: {1, 1, 1, -1, ω, ω²} where ω = exp(2πi/3)
```

### Special 5×5 Pattern (Experimental)

```julia
@variables a b d
mat = [a  b  0  0  0;
       b  a  d  0  0;
       0  d  a  b  0;
       0  0  b  a  b;
       0  0  0  b  a]

vals = eigvals(mat)  # a ± √(2b² + d²), a ± b, a
```

## API

### LinearAlgebra Interface

- `eigen(A; kwargs...)` - Returns `Eigen` object with `.values` and `.vectors`
- `eigvals(A; kwargs...)` - Eigenvalues only (faster)

### Direct API

- `symbolic_eigenvalues(A; kwargs...)` - Returns `(values, poly, λ)`
- `symbolic_eigenpairs(A; kwargs...)` - Returns `[(λ, [v₁, v₂, ...]), ...]`
- `symbolic_diagonalize(A; kwargs...)` - Returns `(P, D, pairs)`

### Options

- `structure` - Hint: `:auto`, `:hermitian`, `:symmetric`, `:unitary`
- `expand` - Expand polynomials (default: `true`)
- `complexity_threshold` - Warn if too many variables (default: 5)
- `timeout` - Max computation time in seconds (default: 300)
- `max_terms` - Expression complexity limit (default: 10000)

## Current Limitations

- **No general 5×5+ solver**: Only works if structure detected
- **Structure detection is basic**: Misses most exploitable patterns
- **Huge expressions**: Fully symbolic 4×4 produces ~13.5 MB symbols per eigenvalue
- **No simplification**: Results are often not in simplest form
- **Fragile**: Many edge cases not handled

## What Needs Work

### High Priority

- [ ] Robust structure detection algorithms
- [ ] Group theory infrastructure (character tables, projection operators)
- [x] Circulant matrix pattern (closed-form via DFT)
- [x] Block circulant matrix pattern (reduction to smaller problems via block DFT)
- [x] Kronecker product pattern (eigenvalues as products of factor eigenvalues)
- [x] Toeplitz tridiagonal pattern (symmetric, closed-form via cosines)
- [x] Anti-diagonal matrix pattern (symmetric, closed-form via ±pairs)
- [x] Permutation matrix pattern (closed-form via cycle decomposition and roots of unity)
- [ ] More special patterns (Hankel, general tridiagonal families, rank-1 updates, block Toeplitz)
- [ ] Better expression simplification
- [ ] Comprehensive testing of edge cases

### Medium Priority

- [ ] Eigenvalue multiplicity handling
- [ ] Symbolic condition number estimation
- [ ] Integration with numerical fallback
- [ ] Documentation of all special patterns

### Future Ideas

- [ ] Machine learning for pattern recognition
- [ ] User-defined pattern libraries
- [ ] Symbolic perturbation theory
- [ ] Parallel processing for large block systems

## Development

### Testing

```julia
julia --project -e 'using Pkg; Pkg.test()'
```

### Exploring Patterns

See [examples/explore_patterns.jl](examples/explore_patterns.jl) for tools to discover new solvable patterns.

See [docs/notes/](docs/notes/) for documentation of discovered patterns.

### Documentation

```julia
julia --project=docs docs/make.jl
```

## Contributing

This is experimental software. Contributions, ideas, and feedback very welcome:

- Implementing known special patterns from linear algebra literature
- Improving structure detection algorithms
- Finding new solvable patterns
- Performance optimization
- Better testing

## License

See LICENSE file for details.
