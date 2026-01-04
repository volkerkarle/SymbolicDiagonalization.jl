# API Reference

Complete API documentation for `SymbolicDiagonalization.jl`.

## Table of Contents

- [User-Facing API](#user-facing-api)
  - [LinearAlgebra Interface](#linearalgebra-interface)
  - [Direct Eigenvalue API](#direct-eigenvalue-api)
  - [Characteristic Polynomial](#characteristic-polynomial)
  - [Root Solvers](#root-solvers)
- [Keyword Arguments](#keyword-arguments)
- [Exception Types](#exception-types)
- [Internal API](#internal-api)

---

## User-Facing API

### LinearAlgebra Interface

The recommended interface for most users. Extends `LinearAlgebra.jl` with symbolic matrix support.

#### `LinearAlgebra.eigen`

```julia
eigen(A; kwargs...) → Eigen
```

Computes eigenvalues and eigenvectors of symbolic matrix A.

**Arguments**:
- `A::Matrix{Union{Num, Complex{Num}}}` - Symbolic matrix to diagonalize

**Keyword Arguments**: See [Keyword Arguments](#keyword-arguments) section

**Returns**: `Eigen` object with fields:
- `.values` - Vector of eigenvalues
- `.vectors` - Matrix of eigenvectors (as columns)

**Example**:
```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)
E.values   # [a, b, c]
E.vectors  # 3×3 eigenvector matrix
```

**Notes**:
- Use `eigvals()` instead if you only need eigenvalues (faster)
- Throws error if matrix is not diagonalizable
- Eigenvectors are columns of the returned matrix

---

#### `LinearAlgebra.eigvals`

```julia
eigvals(A; kwargs...) → Vector
```

Computes eigenvalues of symbolic matrix A (faster than `eigen` since it skips eigenvectors).

**Arguments**:
- `A::Matrix{Union{Num, Complex{Num}}}` - Symbolic matrix

**Keyword Arguments**: See [Keyword Arguments](#keyword-arguments) section

**Returns**: Vector of eigenvalue expressions

**Example**:
```julia
@variables a b
mat = [a b; b a]

λ = eigvals(mat)  # [a + b, a - b]
```

**Notes**:
- Recommended when eigenvectors not needed
- Much faster than `eigen()` for large matrices
- Same structure detection and special pattern solvers

---

### Direct Eigenvalue API

Lower-level API for advanced users who need more control.

#### `symbolic_eigenvalues`

```julia
symbolic_eigenvalues(A; kwargs...) → (values, poly, λ)
```

Computes eigenvalues along with characteristic polynomial.

**Arguments**:
- `A::Matrix` - Symbolic matrix

**Keyword Arguments**: See [Keyword Arguments](#keyword-arguments) section

**Returns**: Tuple of:
- `values` - Vector of eigenvalue expressions
- `poly` - Characteristic polynomial det(λI - A)
- `λ` - Symbolic variable used in polynomial

**Example**:
```julia
@variables a b
mat = [a b; b a]

vals, poly, λ = symbolic_eigenvalues(mat)
# vals: [a+b, a-b]
# poly: (λ - a)² - b²
# λ: symbolic variable
```

**Use cases**:
- Need both eigenvalues and characteristic polynomial
- Want to manipulate polynomial directly
- Need the eigenvalue variable for further computation

---

#### `symbolic_eigenpairs`

```julia
symbolic_eigenpairs(A; kwargs...) → Vector{Tuple{Num, Vector{Vector{Num}}}}
```

Computes eigenvalue-eigenvector pairs with multiplicity handling.

**Arguments**:
- `A::Matrix` - Symbolic matrix

**Keyword Arguments**: 
- `compute_vectors::Bool = true` - Whether to compute eigenvectors
- Plus all [standard kwargs](#keyword-arguments)

**Returns**: Vector of tuples `(eigenvalue, eigenvectors)`:
- Each eigenvalue may have multiple eigenvectors (geometric multiplicity)
- Eigenvectors are vectors, not matrix columns

**Example**:
```julia
@variables a b
mat = [a b; b a]

pairs = symbolic_eigenpairs(mat)
# [(a+b, [v₁]), (a-b, [v₂])]

# Access first eigenvalue and its eigenvectors
λ₁, vecs₁ = pairs[1]
```

**Use cases**:
- Need explicit multiplicity information
- Want eigenvalues grouped with their eigenvectors
- Need to handle degenerate cases

---

#### `symbolic_diagonalize`

```julia
symbolic_diagonalize(A; kwargs...) → (P, D, pairs)
```

Computes full diagonalization A = P D P⁻¹.

**Arguments**:
- `A::Matrix` - Symbolic matrix

**Keyword Arguments**: See [Keyword Arguments](#keyword-arguments) section

**Returns**: Tuple of:
- `P` - Matrix of eigenvectors (columns)
- `D` - Diagonal matrix of eigenvalues
- `pairs` - Same as `symbolic_eigenpairs()` output

**Example**:
```julia
@variables a b
mat = [a b; b a]

P, D, pairs = symbolic_diagonalize(mat)
# Verify: mat ≈ P * D * inv(P)
```

**Throws**:
- `ArgumentError` if matrix is not diagonalizable
- `ArgumentError` if insufficient eigenvectors found

**Use cases**:
- Need explicit matrix factorization
- Want to verify diagonalizability
- Need P and D matrices for further computation

---

### Characteristic Polynomial

#### `characteristic_polynomial`

```julia
characteristic_polynomial(A; var=nothing) → (poly, coeffs, λ)
```

Computes characteristic polynomial det(λI - A) using Bareiss algorithm.

**Arguments**:
- `A::Matrix` - Symbolic matrix
- `var` - Optional symbolic variable to use (auto-generated if not provided)

**Returns**: Tuple of:
- `poly` - Polynomial expression det(λI - A)
- `coeffs` - Vector of polynomial coefficients [c₀, c₁, ..., cₙ]
- `λ` - Symbolic variable

**Example**:
```julia
@variables a b
mat = [a b; b a]

poly, coeffs, λ = characteristic_polynomial(mat)
# poly: λ² - 2a·λ + (a² - b²)
# coeffs: [a² - b², -2a, 1]
```

**Implementation details**:
- Uses **Bareiss fraction-free determinant** for efficiency
- Avoids expression explosion from cofactor expansion
- Coefficients extracted by differentiation at λ = 0

**Use cases**:
- Need characteristic polynomial for analysis
- Want polynomial coefficients explicitly
- Implementing custom root-finding methods

---

### Root Solvers

#### `symbolic_roots`

```julia
symbolic_roots(poly, λ; expand=true, max_terms=10000) → Vector
```

Finds symbolic roots of polynomial using closed-form formulas (degrees 1-4).

**Arguments**:
- `poly` - Polynomial expression in variable λ
- `λ` - Symbolic variable
- `expand::Bool = true` - Whether to expand expressions
- `max_terms::Int = 10000` - Maximum expression complexity

**Returns**: Vector of symbolic root expressions

**Example**:
```julia
@variables λ a b
poly = λ^2 - 2a*λ + (a^2 - b^2)

roots = symbolic_roots(poly, λ)
# [a + b, a - b]
```

**Supported degrees**:
- **Degree 1**: Linear formula
- **Degree 2**: Quadratic formula
- **Degree 3**: Cardano's cubic formula
- **Degree 4**: Ferrari's quartic formula
- **Degree ≥ 5**: Throws error (Abel-Ruffini theorem)

**Throws**:
- `ArgumentError` if degree ≥ 5 and no structure detected
- `ExpressionComplexityError` if expression exceeds `max_terms`

**Notes**:
- Quartic formula produces very large expressions
- Consider structure detection for degree ≥ 4
- For degree 4, expressions can be ~13.5 MB each

---

### Rotation Matrix Constructors

Convenience constructors for rotation matrices with clean symbolic eigenvalues.

#### `SO2_rotation`

```julia
SO2_rotation(θ) → Matrix{Num}
```

Construct a 2×2 rotation matrix (SO(2)).

**Arguments**:
- `θ` - Rotation angle (symbolic or numeric)

**Returns**: 2×2 rotation matrix `[cos(θ) -sin(θ); sin(θ) cos(θ)]`

**Example**:
```julia
@variables θ
R = SO2_rotation(θ)
eigvals(R)  # [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
```

---

#### `SO3_Rx`, `SO3_Ry`, `SO3_Rz`

```julia
SO3_Rx(θ) → Matrix{Num}  # Rotation around x-axis
SO3_Ry(θ) → Matrix{Num}  # Rotation around y-axis
SO3_Rz(θ) → Matrix{Num}  # Rotation around z-axis
```

Construct 3×3 elementary rotation matrices in SO(3).

**Arguments**:
- `θ` - Rotation angle (symbolic or numeric)

**Returns**: 3×3 rotation matrix

**Example**:
```julia
@variables θ
eigvals(SO3_Rx(θ))  # [1, cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
eigvals(SO3_Ry(θ))  # [1, cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
eigvals(SO3_Rz(θ))  # [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ), 1]
```

---

#### `SO2_kron`

```julia
SO2_kron(angles::Vector) → Matrix{Num}
```

Construct the Kronecker product of SO(2) rotation matrices.

**Arguments**:
- `angles` - Vector of rotation angles `[θ₁, θ₂, ..., θₖ]`

**Returns**: 2ᵏ × 2ᵏ matrix `SO2_rotation(θ₁) ⊗ SO2_rotation(θ₂) ⊗ ... ⊗ SO2_rotation(θₖ)`

**Example**:
```julia
@variables α β γ
K = SO2_kron([α, β, γ])  # 8×8 matrix
size(K)  # (8, 8)
```

---

#### `SO2_kron_eigenvalues`

```julia
SO2_kron_eigenvalues(angles::Vector) → Vector
```

Compute eigenvalues of SO(2) Kronecker products directly in clean trigonometric form.

**Arguments**:
- `angles` - Vector of rotation angles `[θ₁, θ₂, ..., θₖ]`

**Returns**: Vector of 2ᵏ eigenvalues `cos(±θ₁±θ₂±...±θₖ) + i·sin(±θ₁±θ₂±...±θₖ)`

**Example**:
```julia
@variables α β
vals = SO2_kron_eigenvalues([α, β])
# [cos(-α-β) + im*sin(-α-β),
#  cos(α-β) + im*sin(α-β),
#  cos(-α+β) + im*sin(-α+β),
#  cos(α+β) + im*sin(α+β)]
```

**Notes**:
- Faster than building the matrix and computing eigenvalues
- Always produces clean form without messy symbolic expressions
- The same clean form is automatically detected when using `eigvals(SO2_kron(angles))`

---

### SU(2) Constructors

#### `pauli_x`, `pauli_y`, `pauli_z`

```julia
pauli_x() → Matrix{Complex{Int}}  # Pauli σx matrix
pauli_y() → Matrix{Complex{Int}}  # Pauli σy matrix
pauli_z() → Matrix{Int}           # Pauli σz matrix
```

Return the three Pauli matrices.

**Example**:
```julia
σx = pauli_x()  # [0 1; 1 0]
σy = pauli_y()  # [0 -im; im 0]
σz = pauli_z()  # [1 0; 0 -1]
```

---

#### `SU2_Ux`, `SU2_Uy`, `SU2_Uz`

```julia
SU2_Ux(θ) → Matrix{Complex{Num}}  # exp(-i θ σx/2)
SU2_Uy(θ) → Matrix{Num}           # exp(-i θ σy/2)
SU2_Uz(θ) → Matrix{Complex{Num}}  # exp(-i θ σz/2)
```

Construct SU(2) rotation matrices (spin-1/2 representation).

**Arguments**:
- `θ` - Rotation angle (symbolic or numeric)

**Example**:
```julia
@variables θ
U = SU2_Uz(θ)  # diagonal: [exp(-iθ/2), exp(iθ/2)]
```

---

#### `SU2_kron`

```julia
SU2_kron(angles::Vector) → Matrix{Complex{Num}}
```

Construct the Kronecker product of SU(2) rotation matrices.

**Arguments**:
- `angles` - Vector of rotation angles `[α, β, ...]`

**Returns**: 2ᵏ × 2ᵏ matrix `SU2_Uz(α) ⊗ SU2_Uz(β) ⊗ ...`

**Example**:
```julia
@variables α β
K = SU2_kron([α, β])  # 4×4 matrix
```

---

#### `SU2_kron_eigenvalues`

```julia
SU2_kron_eigenvalues(angles::Vector) → Vector
```

Compute eigenvalues of SU(2) Kronecker products directly.

**Returns**: Eigenvalues with half-angle formulas: `cos((±α±β±...)/2) + i·sin((±α±β±...)/2)`

---

### SU(3) Constructors

#### `gellmann_1` through `gellmann_8`

```julia
gellmann_1() → Matrix  # First Gell-Mann matrix
...
gellmann_8() → Matrix  # Eighth Gell-Mann matrix
```

Return the eight Gell-Mann matrices (generators of SU(3)).

---

#### `gellmann_matrices`

```julia
gellmann_matrices() → Vector{Matrix}
```

Return all eight Gell-Mann matrices as a vector.

---

#### `SU3_diagonal_trig`

```julia
SU3_diagonal_trig(α₁, α₂) → Matrix{Complex{Num}}
```

Construct a diagonal SU(3) matrix (Cartan subalgebra element).

**Returns**: `diag(exp(iα₁), exp(iα₂), exp(-i(α₁+α₂)))`

---

#### `SU3_kron`

```julia
SU3_kron((α₁, α₂), (β₁, β₂)) → Matrix{Complex{Num}}
```

Construct the Kronecker product of two diagonal SU(3) matrices.

**Returns**: 9×9 diagonal matrix

---

#### `SU3_kron_eigenvalues`

```julia
SU3_kron_eigenvalues((α₁, α₂), (β₁, β₂)) → Vector
```

Compute eigenvalues of SU(3) Kronecker products directly.

**Returns**: 9 eigenvalues of the form `cos(θ) + i·sin(θ)` where θ is a phase sum.

---

### Lie Algebra Generators

#### `spin_j_generators`

```julia
spin_j_generators(j) → (Jx, Jy, Jz)
```

Construct the standard spin-j generators for the (2j+1)-dimensional irreducible representation of su(2)/so(3).

**Arguments**:
- `j` - Spin value (half-integer: 1//2, 1, 3//2, etc.)

**Returns**: Tuple of three matrices (Jx, Jy, Jz)

**Example**:
```julia
Jx, Jy, Jz = spin_j_generators(1//2)  # 2×2 matrices
Jx, Jy, Jz = spin_j_generators(1)      # 3×3 matrices
```

---

#### `so3_generators`

```julia
so3_generators() → (Lx, Ly, Lz)
```

Return the three so(3) generators (3×3 skew-symmetric matrices).

These satisfy `[Lᵢ, Lⱼ] = εᵢⱼₖ Lₖ`.

---

#### `su2_generators`

```julia
su2_generators() → (τ1, τ2, τ3)
```

Return the three su(2) generators: `τₖ = i·σₖ/2`.

These are skew-Hermitian and satisfy `[τᵢ, τⱼ] = εᵢⱼₖ τₖ`.

---

### Trigonometric Simplification

#### `trig_simplify`

```julia
trig_simplify(expr) → Num
```

Apply trigonometric simplification rules to a symbolic expression using SymbolicUtils.jl
rule-based rewriting.

**Arguments**:
- `expr` - Symbolic expression containing trigonometric functions

**Returns**: Simplified expression

**Example**:
```julia
@variables θ φ
trig_simplify(sin(θ)^2 + cos(θ)^2)           # → 1
trig_simplify(cos(θ)*cos(φ) - sin(θ)*sin(φ)) # → cos(θ + φ)
trig_simplify(2*sin(θ)*cos(θ))               # → sin(2θ)
```

**Supported identities**:
- Pythagorean: $\sin^2(x) + \cos^2(x) = 1$
- Squared Pythagorean: $\sin^4(x) + 2\sin^2(x)\cos^2(x) + \cos^4(x) = 1$
- Angle addition: $\cos(a)\cos(b) - \sin(a)\sin(b) = \cos(a+b)$
- Angle subtraction: $\cos(a)\cos(b) + \sin(a)\sin(b) = \cos(a-b)$
- Double angle: $2\sin(x)\cos(x) = \sin(2x)$, $\cos^2(x) - \sin^2(x) = \cos(2x)$

**Use cases**:
- Simplifying eigenvalues of rotation matrices
- Cleaning up trigonometric expressions in physics problems
- Post-processing SO(2) Kronecker product eigenvalues

---

#### `aggressive_simplify`

```julia
aggressive_simplify(expr) → Num
```

Apply aggressive simplification combining trigonometric, algebraic, and sqrt-trig rules.

**Key feature**: Transforms `sqrt(1 - cos(θ)^2)` → `sin(θ)` and `sqrt(1 - sin(θ)^2)` → `cos(θ)`.

**Example**:
```julia
@variables θ
aggressive_simplify(sqrt(1 - cos(θ)^2))  # → sin(θ)
aggressive_simplify(sqrt(1 - sin(θ)^2))  # → cos(θ)
aggressive_simplify(sin(θ)^2 + cos(θ)^2) # → 1
```

---

#### `simplify_eigenvalue`

```julia
simplify_eigenvalue(expr) → Complex{Num}
```

Specialized simplification for eigenvalue expressions. Handles real and imaginary parts separately.

**Example**:
```julia
@variables θ
expr = cos(θ) + im*sqrt(1 - cos(θ)^2)
simplify_eigenvalue(expr)  # → cos(θ) + im*sin(θ)
```

---

#### `simplify_eigenvalues`

```julia
simplify_eigenvalues(vals::Vector) → Vector
```

Apply `simplify_eigenvalue` to a vector of eigenvalues.

**Example**:
```julia
@variables θ
vals = [cos(θ) + im*sqrt(1 - cos(θ)^2), cos(θ) - im*sqrt(1 - cos(θ)^2)]
simplify_eigenvalues(vals)  # → [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
```

---

## Keyword Arguments

All main functions (`eigen`, `eigvals`, `symbolic_eigenvalues`, etc.) accept these keyword arguments:

### `var::Union{Nothing, Num}`
**Default**: `nothing` (auto-generate)

Symbolic variable to use for eigenvalue. If `nothing`, a fresh variable is created.

**Example**:
```julia
@variables λ
vals = eigvals(mat, var=λ)
```

**Use cases**: 
- Need specific variable name
- Integrating with existing symbolic expressions
- Want to control variable scope

---

### `structure::Symbol`
**Default**: `:auto`

Hint about matrix structure for optimization.

**Options**:
- `:auto` - Automatic structure detection (default)
- `:hermitian` - Hermitian matrix (A† = A)
- `:symmetric` - Real symmetric (Aᵀ = A)
- `:unitary` - Unitary matrix (A† = A⁻¹)
- `:general` - No structure assumptions
- `:diagonal` - Diagonal matrix (hint to skip detection)
- `:triangular` - Triangular matrix (hint to skip detection)
- `:none` - Skip all structure detection entirely

**Example**:
```julia
# Tell solver that matrix is Hermitian
vals = eigvals(hermitian_mat, structure=:hermitian)
```

**Notes**:
- Providing correct hint can speed up computation
- Incorrect hint may give wrong results
- `:auto` is safe but may be slower
- `:none` skips detection entirely (fastest, but no pattern optimizations)

---

### `expand::Bool`
**Default**: `true`

Whether to expand polynomial expressions.

**Options**:
- `true` - Expand products and powers
- `false` - Keep factored form when possible

**Example**:
```julia
# Factored form (if available)
vals = eigvals(mat, expand=false)
```

**Use cases**:
- Factored form may be simpler
- Expanded form needed for numerical evaluation
- Debugging expression structure

---

### `complexity_threshold::Int`
**Default**: `5`

Warn if matrix contains more than this many symbolic variables.

**Example**:
```julia
# Suppress warning for 10-variable matrix
vals = eigvals(large_mat, complexity_threshold=10)
```

**Purpose**: Prevent accidentally running huge symbolic computations that may timeout or produce enormous expressions.

**Recommendation**: 
- Keep ≤ 5 for general matrices
- Can increase for structured matrices
- Set to `Inf` to disable warnings

---

### `timeout::Int`
**Default**: `300` (5 minutes)

Maximum computation time in seconds.

**Example**:
```julia
# Allow only 60 seconds
vals = eigvals(mat, timeout=60)
```

**Throws**: `ComputationTimeoutError` if exceeded

**Use cases**:
- Prevent runaway computations
- Testing with time limits
- Interactive use with quick feedback

---

### `max_terms::Int`
**Default**: `10000`

Maximum number of terms allowed in symbolic expressions.

**Example**:
```julia
# Allow more complex expressions
vals = eigvals(mat, max_terms=50000)
```

**Throws**: `ExpressionComplexityError` if exceeded

**Purpose**: Prevent expression explosion (especially in quartic formula)

**Recommendation**:
- Keep default for most cases
- Increase for structured matrices
- Watch memory usage if increasing

---

## Exception Types

### `ExpressionComplexityError`

Thrown when symbolic expression exceeds complexity limit.

**Fields**:
- `message::String` - Error message with details and suggestions for resolution

**Example**:
```julia
try
    vals = eigvals(huge_mat, max_terms=1000)
catch e
    if e isa ExpressionComplexityError
        println("Expression too complex: $(e.message)")
    end
end
```

**Common causes**:
- Fully symbolic 4×4 matrices (quartic formula)
- Lack of detected structure
- Nested radical expressions

**Solutions**:
- Add structure to matrix
- Use partial numeric substitution
- Increase `max_terms` (watch memory!)
- Use numerical methods instead

---

### `ComputationTimeoutError`

Thrown when computation exceeds time limit.

**Fields**:
- `message::String` - Error message with details and suggestions for resolution

**Example**:
```julia
try
    vals = eigvals(mat, timeout=30)
catch e
    if e isa ComputationTimeoutError
        println("Timed out: $(e.message)")
    end
end
```

**Common causes**:
- Very large symbolic matrices
- Complex expression simplification
- Expensive structure detection

**Solutions**:
- Increase timeout
- Simplify matrix structure
- Use numerical methods
- Pre-substitute some variables

---

## Internal API

These functions are not exported but may be useful for advanced users or contributors.

### Structure Detection

#### `_detect_structure(mat)`
Automatically detects matrix structure (diagonal, triangular, block-diagonal, etc.)

#### `_is_diagonal(mat)`
Tests if matrix is diagonal

#### `_is_triangular(mat)` 
Tests if matrix is upper or lower triangular

#### `_is_hermitian(mat)`
Tests if matrix is Hermitian (A† = A)

#### `_is_symmetric(mat)`
Tests if matrix is real symmetric

#### `_is_persymmetric(mat)`
Tests if matrix has persymmetric structure Q[i,j] = Q[n+1-j,n+1-i]

---

### Special Pattern Detectors

#### `_is_circulant(mat)`
Tests if matrix is circulant (each row cyclic shift of previous)

#### `_is_block_circulant(mat)`
Tests if matrix is block circulant, returns `(true, n_blocks, block_size, blocks)` or `(false, ...)`

#### `_is_kronecker_product(mat)`
Tests if matrix is Kronecker product, returns `(true, A, B, m, n)` or `(false, ...)`

#### `_is_toeplitz_tridiagonal(mat)`
Tests if matrix is symmetric Toeplitz tridiagonal

#### `_is_antidiagonal(mat)`
Tests if matrix is symmetric anti-diagonal

#### `_is_permutation_matrix(mat)`
Tests if matrix is a permutation matrix

#### `_detect_special_5x5_tridiagonal(mat)`
Tests for special 5×5 patterns with known eigenvalues

#### `_is_symbolic_orthogonal(mat)`
Tests if a symbolic matrix is orthogonal (M'*M = I) using trig simplification

#### `_try_so2_kronecker_decomposition(mat, k)`
Attempts to detect and decompose SO(2) Kronecker products R(θ₁) ⊗ R(θ₂) ⊗ ... ⊗ R(θₖ)

---

### Special Pattern Solvers

#### `_circulant_eigenvalues(mat)`
Computes eigenvalues of circulant matrix using DFT

#### `_block_circulant_eigenvalues(mat, n_blocks, block_size, blocks; ...)`
Computes eigenvalues of block circulant matrix

#### `_kronecker_eigenvalues(A, B, m, n; ...)`
Computes eigenvalues of Kronecker product A ⊗ B

#### `_toeplitz_tridiagonal_eigenvalues(n, a, b, c)`
Computes eigenvalues of symmetric Toeplitz tridiagonal

#### `_antidiagonal_eigenvalues(mat)`
Computes eigenvalues of symmetric anti-diagonal matrix

#### `_compute_permutation_eigenvalues(A)`
Computes eigenvalues of permutation matrix via cycle decomposition

#### `_issymzero_trig(expr)`
Tests if expression is symbolically zero, using trigonometric simplification

---

### Utility Functions

#### `_block_split(mat)`
Detects and splits block-diagonal structure

#### `_persymmetric_split(mat)`
Splits persymmetric matrix into half-sized problems

#### `_count_symbolic_vars(A)`
Counts number of unique symbolic variables in matrix

#### `_check_complexity(A; threshold, quiet)`
Checks if matrix exceeds complexity threshold, optionally warns

---

## Usage Examples

### Example 1: Basic Eigenvalues

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 0 0; 0 b 0; 0 0 c]

λ = eigvals(mat)  # [a, b, c]
```

### Example 2: With Structure Hint

```julia
@variables a b
H = [a b; conj(b) a]  # Hermitian

E = eigen(H, structure=:hermitian)
```

### Example 3: Error Handling

```julia
@variables a b c d e
mat = [a b c d e;
       b a c d e;
       c c a d e;
       d d d a e;
       e e e e a]

try
    λ = eigvals(mat, timeout=60, max_terms=5000)
catch e
    if e isa ComputationTimeoutError
        println("Computation timed out")
    elseif e isa ExpressionComplexityError
        println("Expression too complex")
    else
        rethrow(e)
    end
end
```

### Example 4: Characteristic Polynomial

```julia
@variables a b
mat = [a b; b a]

poly, coeffs, λ = characteristic_polynomial(mat)
# poly = λ² - 2a·λ + (a² - b²)

# Now solve manually if desired
roots = symbolic_roots(poly, λ)
```

### Example 5: Full Diagonalization

```julia
@variables a b
mat = [a b; b a]

P, D, pairs = symbolic_diagonalize(mat)

# Verify: mat ≈ P * D * inv(P)
# P contains eigenvectors as columns
# D is diagonal with eigenvalues
```

---

## See Also

- [User Guide](user_guide.md) - Practical examples and workflows
- [Pattern Library](pattern_library.md) - Special patterns and their eigenvalues
- [Implementation](implementation.md) - Algorithm details
- [Mathematical Background](mathematical_background.md) - Theory and proofs
