# Contributing to SymbolicDiagonalization.jl

Thank you for your interest in contributing! This guide covers development setup, coding standards, testing requirements, and how to add new patterns.

## Development Setup

### Prerequisites

- Julia ≥1.6
- Git (if the project becomes a repository)
- Recommended: VSCode with Julia extension

### Installation for Development

```bash
# Clone or navigate to the package directory
cd SymbolicDiagonalization.jl

# Activate the project environment
julia --project=.

# Install dependencies
julia> using Pkg
julia> Pkg.instantiate()

# Run tests to verify setup
julia> Pkg.test()
```

### Project Structure

```
SymbolicDiagonalization.jl/
├── src/
│   ├── SymbolicDiagonalization.jl  # Main module
│   ├── characteristic_poly.jl       # Bareiss determinant
│   ├── root_solvers.jl              # Closed-form root formulas (1-4)
│   ├── structure_detection.jl       # Pattern detection algorithms
│   ├── special_patterns.jl          # Pattern-specific eigensolvers
│   ├── eigenvectors.jl              # Nullspace computation
│   └── utils.jl                     # Helper functions
├── test/
│   ├── runtests.jl                  # Test suite entry point
│   ├── test_basic.jl                # Small matrix tests
│   ├── test_structures.jl           # Pattern-specific tests
│   └── test_edge_cases.jl           # Edge case coverage
├── docs/
│   ├── src/                         # Documentation source
│   └── make.jl                      # Documentation builder
└── examples/
    └── explore_patterns.jl          # Interactive pattern exploration
```

## Code Style Guidelines

### General Principles

- **Clarity over cleverness**: Readable code is maintainable code
- **Type stability**: Avoid type instabilities for performance
- **Document intent**: Comments should explain *why*, not *what*
- **Fail fast**: Use early error checks with clear messages

### Naming Conventions

- **Functions**: `snake_case` (e.g., `detect_block_diagonal`)
- **Types**: `PascalCase` (e.g., `EigenPair`)
- **Constants**: `SCREAMING_SNAKE_CASE` (e.g., `MAX_ITERATIONS`)
- **Private functions**: Prefix with `_` (e.g., `_compute_helper`)

### Code Formatting

```julia
# Good: Clear structure, well-documented
function detect_circulant(A::AbstractMatrix{T}) where T
    n = size(A, 1)
    n == size(A, 2) || return false
    
    # Check if all rows are cyclic shifts of first row
    for i in 2:n
        for j in 1:n
            if A[i, j] != A[1, mod1(j - i + 1, n)]
                return false
            end
        end
    end
    return true
end

# Bad: Dense, unclear, no documentation
function detect_circulant(A::AbstractMatrix{T}) where T;n=size(A,1);n==size(A,2)||return false;for i in 2:n;for j in 1:n;A[i,j]!=A[1,mod1(j-i+1,n)]&&return false;end;end;true;end
```

### Documentation Strings

Every exported function must have a docstring following this template:

```julia
"""
    function_name(arg1, arg2; kwarg1=default)

Brief one-line description.

Extended description with usage details, algorithm notes, and complexity.

# Arguments
- `arg1`: Description of arg1
- `arg2`: Description of arg2

# Keywords
- `kwarg1`: Description and default value

# Returns
Description of return value(s)

# Examples
```julia
julia> using SymbolicDiagonalization, Symbolics
julia> @variables a b
julia> result = function_name([a b; b a])
[a + b, a - b]
```

# Complexity
Time and space complexity analysis

# References
- [Author Year] Paper title, Journal/Conference
"""
function function_name(arg1, arg2; kwarg1=default)
    # Implementation
end
```

## Testing Requirements

### Test Coverage Standards

- **All new functions**: Must have at least one test
- **Pattern detectors**: Test true positives AND false negatives
- **Eigensolvers**: Test correctness (A*v = λ*v) and completeness (all eigenvalues found)
- **Edge cases**: Empty matrices, 1×1 matrices, repeated eigenvalues

### Writing Tests

Use the `@testset` macro to organize tests:

```julia
@testset "Circulant Matrices" begin
    @testset "Detection" begin
        # Test positive detection
        @variables a b c
        circ = [a b c; c a b; b c a]
        @test detect_circulant(circ) == true
        
        # Test negative detection
        non_circ = [a b c; c a b; a c b]
        @test detect_circulant(non_circ) == false
    end
    
    @testset "Eigenvalues" begin
        # Test known eigenvalue formula
        circ = [1 2 3; 3 1 2; 2 3 1]
        λ = eigvals(circ)
        
        # Verify A*v = λ*v for each eigenvalue
        for (i, val) in enumerate(λ)
            v = eigvecs(circ)[:, i]
            @test circ * v ≈ val * v atol=1e-10
        end
    end
    
    @testset "Size Scaling" begin
        # Test pattern works for various sizes
        for n in [3, 5, 7, 10]
            A = circulant_matrix(n)
            λ = eigvals(A)
            @test length(λ) == n
        end
    end
end
```

### Running Tests

```julia
# Run full test suite
julia> using Pkg; Pkg.test()

# Run specific test file
julia> include("test/test_structures.jl")

# Run with coverage
julia --code-coverage=user -e 'using Pkg; Pkg.test()'
```

## Adding New Patterns

Follow this 5-step process to add a new matrix pattern:

### Step 1: Detection Function

Create a detector in `src/structure_detection.jl`:

```julia
"""
    detect_my_pattern(A::AbstractMatrix) -> Bool

Detect if matrix A has MyPattern structure.

MyPattern is defined as [mathematical definition].

# Complexity
O(n²) - Checks all matrix entries
"""
function detect_my_pattern(A::AbstractMatrix{T}) where T
    n = size(A, 1)
    n == size(A, 2) || return false
    
    # Check pattern-specific properties
    for i in 1:n, j in 1:n
        # Pattern condition
        if !condition(A[i,j])
            return false
        end
    end
    return true
end
```

### Step 2: Eigensolver Function

Create a solver in `src/special_patterns.jl`:

```julia
"""
    solve_my_pattern(A::AbstractMatrix; kwargs...) -> Vector

Compute eigenvalues of MyPattern matrices using [algorithm name].

# Algorithm
Brief description of the algorithm and its basis (e.g., DFT, closed-form formula).

# Complexity
- Time: O(f(n))
- Space: O(g(n))

# References
- [Author Year] Paper title
"""
function solve_my_pattern(A::AbstractMatrix{T}; kwargs...) where T
    n = size(A, 1)
    
    # Extract pattern parameters
    params = extract_parameters(A)
    
    # Compute eigenvalues using pattern-specific algorithm
    eigenvalues = compute_eigenvalues(params)
    
    return eigenvalues
end
```

### Step 3: Integration

Add pattern to detection chain in `src/structure_detection.jl`:

```julia
function detect_and_solve(A::AbstractMatrix; kwargs...)
    # ... existing patterns ...
    
    # Add new pattern
    if detect_my_pattern(A)
        return solve_my_pattern(A; kwargs...)
    end
    
    # ... fallback ...
end
```

### Step 4: Testing

Create comprehensive tests in `test/test_structures.jl`:

```julia
@testset "MyPattern Matrices" begin
    @testset "Detection" begin
        # True positive
        @variables a b
        mat = construct_my_pattern(a, b)
        @test detect_my_pattern(mat) == true
        
        # False negative
        non_mat = [a b; b a]
        @test detect_my_pattern(non_mat) == false
    end
    
    @testset "Known Eigenvalues" begin
        # Test with known theoretical eigenvalues
        A = [concrete example]
        λ = eigvals(A)
        expected = [known eigenvalues]
        @test sort(λ) ≈ sort(expected)
    end
    
    @testset "Eigenvector Verification" begin
        # Verify A*v = λ*v
        A = construct_my_pattern(...)
        E = eigen(A)
        for i in 1:length(E.values)
            @test A * E.vectors[:, i] ≈ E.values[i] * E.vectors[:, i]
        end
    end
    
    @testset "Size Scaling" begin
        # Test multiple sizes
        for n in [3, 5, 10, 20]
            A = construct_my_pattern(n)
            λ = eigvals(A)
            @test length(λ) == n
        end
    end
end
```

### Step 5: Documentation

Add pattern to `docs/src/pattern_library.md`:

```markdown
## MyPattern

### Definition
Mathematical definition with notation.

### Properties
- Property 1
- Property 2

### Eigenvalue Formula
Closed-form expression or algorithm description.

### Example
\```julia
@variables a b
A = [matrix construction]
λ = eigvals(A)
# Result: [eigenvalues]
\```

### Complexity
- Detection: O(n²)
- Eigenvalues: O(f(n))
- Eigenvectors: O(g(n))

### Implementation Details
Algorithm notes, numerical considerations, limitations.

### References
- [Author Year] Paper title
```

## Pattern Submission Template

When proposing a new pattern, open an issue with:

```markdown
## Pattern Name: [Name]

### Mathematical Definition
LaTeX or clear description of the pattern

### Eigenvalue Formula
Known closed-form formula or algorithm

### Complexity
Expected time/space complexity

### References
Papers or textbooks describing the pattern

### Example
Concrete example matrix and its eigenvalues

### Implementation Notes
Any special considerations (numerical stability, edge cases, etc.)
```

## Git Workflow (When Applicable)

If this project becomes a Git repository:

### Branch Naming

- `feature/pattern-name` - New pattern implementations
- `fix/issue-description` - Bug fixes
- `docs/topic` - Documentation updates
- `perf/optimization-area` - Performance improvements

### Commit Messages

Follow conventional commits:

```
type(scope): Brief description

Detailed explanation of changes and motivation.

Fixes #issue-number
```

Types: `feat`, `fix`, `docs`, `test`, `perf`, `refactor`, `style`

### Pull Request Process

1. **Fork and branch** from `main`
2. **Implement** with tests and documentation
3. **Run tests** - All must pass
4. **Update docs** - Add examples and API references
5. **Open PR** with clear description
6. **Address review** comments
7. **Squash and merge** when approved

## Code Review Guidelines

### For Reviewers

- **Check correctness**: Verify algorithm implementation against references
- **Test coverage**: Ensure all code paths are tested
- **Documentation**: Confirm docstrings and examples are clear
- **Performance**: Look for obvious inefficiencies
- **Style**: Check adherence to guidelines

### For Authors

- **Self-review first**: Read your own changes critically
- **Respond promptly**: Address feedback quickly
- **Be open**: Accept constructive criticism gracefully
- **Test thoroughly**: Don't rely on reviewers to find bugs

## Performance Optimization Guidelines

### Profiling

```julia
using Profile, ProfileView

# Profile a function
@profile eigvals(large_matrix)
ProfileView.view()

# Benchmark with BenchmarkTools
using BenchmarkTools
@benchmark eigvals($matrix)
```

### Common Optimizations

1. **Avoid allocations**: Use in-place operations where possible
2. **Type stability**: Ensure functions return consistent types
3. **Loop order**: Access arrays in column-major order
4. **SIMD**: Use `@simd` for vectorizable loops
5. **Preallocation**: Allocate output arrays before loops

### Example Optimization

```julia
# Slow: Allocates in loop
function slow_transform(A)
    result = []
    for i in 1:size(A, 1)
        push!(result, transform(A[i, :]))
    end
    return result
end

# Fast: Preallocated, type-stable
function fast_transform(A::Matrix{T}) where T
    n = size(A, 1)
    result = Vector{T}(undef, n)
    for i in 1:n
        result[i] = transform(view(A, i, :))
    end
    return result
end
```

## Documentation Building

### Local Documentation

```bash
# Build documentation
julia --project=docs docs/make.jl

# View in browser with HTTP server (required for navigation links to work)
cd docs
./serve.sh
# Then open http://localhost:8000
```

### Documentation Style

- **Be concise**: Short sentences, clear structure
- **Use examples**: Show, don't just tell
- **Link references**: Cross-reference related functions
- **Update regularly**: Keep docs in sync with code

## CI/CD Setup for GitHub

### Automatic Documentation Deployment

The package includes GitHub Actions workflow for automatic documentation deployment to GitHub Pages.

#### Initial Setup (One-Time)

**1. Generate SSH Deploy Key**

```bash
# Install DocumenterTools if needed
julia -e 'using Pkg; Pkg.add("DocumenterTools")'

# Generate deploy key
julia -e 'using DocumenterTools; DocumenterTools.genkeys()'
```

This creates:
- `docs/src/.documenter` (private key - do NOT commit)
- Public key displayed in terminal

**2. Add Deploy Key to GitHub**

1. Go to your repository on GitHub
2. Navigate to: **Settings** → **Deploy keys** → **Add deploy key**
3. Title: `documenter-key`
4. Key: Paste the **public key** from step 1
5. ✅ Check "Allow write access"
6. Click **Add key**

**3. Add Private Key as GitHub Secret**

1. Copy the contents of `docs/src/.documenter` (the private key file)
2. Go to: **Settings** → **Secrets and variables** → **Actions**
3. Click **New repository secret**
4. Name: `DOCUMENTER_KEY`
5. Value: Paste the entire private key contents
6. Click **Add secret**

**4. Update Repository URLs in `docs/make.jl`**

Replace `YOUR_USERNAME` with your GitHub username:

```julia
repo = "https://github.com/YOUR_USERNAME/SymbolicDiagonalization.jl/blob/{commit}{path}#{line}",
canonical = "https://YOUR_USERNAME.github.io/SymbolicDiagonalization.jl",
```

And:

```julia
deploydocs(;
    repo = "github.com/YOUR_USERNAME/SymbolicDiagonalization.jl",
```

**5. Enable GitHub Pages**

1. Go to: **Settings** → **Pages**
2. Source: **Deploy from a branch**
3. Branch: `gh-pages` / `root`
4. Click **Save**

#### How It Works

- **On push to `main`**: Builds and deploys docs to `gh-pages` branch
- **On pull requests**: Builds docs to verify they work (doesn't deploy)
- **On tags**: Deploys versioned documentation

Your documentation will be available at:
```
https://YOUR_USERNAME.github.io/SymbolicDiagonalization.jl/
```

#### Testing CI Locally

You can test that the CI build will work:

```bash
# Simulate CI environment
CI=true julia --project=docs -e '
    using Pkg
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
'
CI=true julia --project=docs docs/make.jl
```

#### Troubleshooting

**"SSH key authentication failed"**
- Ensure `DOCUMENTER_KEY` secret contains the entire private key
- Verify deploy key has write access enabled

**"gh-pages branch not found"**
- First deployment creates the branch automatically
- Wait a few minutes after first push

**"Documentation build failed"**
- Check the Actions tab for error logs
- Ensure all dependencies are in `docs/Project.toml`
- Test locally with `CI=true` as shown above

**"Pages not updating"**
- Check GitHub Pages is enabled and pointing to `gh-pages`
- May take 5-10 minutes for changes to appear
- Check Pages build status in Actions tab

## Getting Help

### Resources

- **Julia Discourse**: https://discourse.julialang.org/
- **Julia Slack**: Get invited at https://julialang.org/slack/
- **Symbolics.jl Docs**: https://symbolics.juliasymbolics.org/

### Questions?

Open an issue with:
- Clear description of the problem
- Minimal reproducible example
- Julia version and package versions (`Pkg.status()`)
- Expected vs actual behavior

## License

Contributions are assumed to be licensed under the same license as the package (MIT).

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to SymbolicDiagonalization.jl!
