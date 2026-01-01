# Contributing

Guide for contributing to SymbolicDiagonalization.jl.

## Setup

```bash
cd SymbolicDiagonalization.jl
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()
Pkg.test()
```

## Code Style

- **Functions**: `snake_case`, prefix private with `_`
- **Types**: `PascalCase`
- **Clarity over cleverness**: Readable code is maintainable
- **Document intent**: Comments explain *why*, not *what*

### Docstrings

```julia
"""
    function_name(arg1; kwarg1=default)

Brief description.

# Arguments
- `arg1`: Description

# Returns
Description of return value

# Examples
```julia
result = function_name(input)
```
"""
```

## Testing

```julia
# Run all tests
Pkg.test()

# Run specific file
include("test/test_patterns.jl")
```

### Test Requirements

- All new functions must have tests
- Test both detection (true/false positives) and correctness (A*v = λ*v)
- Test edge cases: 1×1, repeated eigenvalues, etc.

## Adding a Pattern

### 1. Detection

Add `_is_pattern(mat)` in appropriate `src/patterns/*.jl`:

```julia
function _is_my_pattern(mat)
    # Return nothing if not detected
    # Return detection info if detected
end
```

### 2. Solver

Add `_my_pattern_eigenvalues(mat, info)`:

```julia
function _my_pattern_eigenvalues(mat, info)
    # Return vector of eigenvalues
end
```

### 3. Integration

Add to detection cascade in `src/diagonalize.jl`:

```julia
# In symbolic_eigenvalues_core:
info = _is_my_pattern(mat)
if !isnothing(info)
    return _my_pattern_eigenvalues(mat, info), nothing, nothing
end
```

### 4. Tests

Add to `test/test_patterns.jl`:

```julia
@testset "MyPattern" begin
    @test _is_my_pattern(valid_matrix) !== nothing
    @test _is_my_pattern(invalid_matrix) === nothing
    
    vals = eigvals(valid_matrix)
    @test length(vals) == n
    # Verify eigenvalue equation
end
```

### 5. Documentation

Add section to `docs/src/pattern_library.md`.

## Building Docs

```bash
julia --project=docs docs/make.jl

# Serve locally (required for working navigation)
cd docs/build && python3 -m http.server 8000
```

## Git Workflow

- Branch: `feature/pattern-name`, `fix/description`
- Commit: `type(scope): description` (feat, fix, docs, test)
- All tests must pass before merge

## Questions

Open an issue with:
- Clear problem description
- Minimal reproducible example
- Julia version (`versioninfo()`)

## License

Contributions are MIT licensed.
