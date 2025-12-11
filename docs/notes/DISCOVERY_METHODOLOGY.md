# Methodology for Discovering Closed-Form Eigenvalue Patterns

## Executive Summary

We developed a systematic approach to discover tridiagonal matrix patterns that have closed-form eigenvalue expressions, despite exceeding the degree limit imposed by the Abel-Ruffini theorem. This document describes the methodology, findings, and future directions.

## The Challenge

For matrices larger than 4×4, computing symbolic eigenvalues requires solving polynomial equations of degree ≥5, which generally have no closed-form solutions. However, certain *structured* matrices have special properties that allow closed-form eigenvalues even at large sizes.

## Methodology: A Systematic Discovery Process

### Step 1: Choose a Starting Point

Begin with **well-known patterns** that have closed-form eigenvalues:

1. **Constant Tridiagonal (Toeplitz)**: All diagonals = a, all off-diagonals = b
   - Eigenvalues: λₖ = a + 2b·cos(kπ/(n+1))
   - Well-studied in linear algebra literature

2. **Diagonal matrices**: Trivial closed-form (the diagonal elements)

3. **Block-diagonal matrices**: Eigenvalues of each block

### Step 2: Apply Systematic Perturbations

Starting from known patterns, apply small, structured perturbations:

**Single Perturbation:**
- Change ONE off-diagonal element from the constant value
- Test all positions: 1st, 2nd, ..., (n-1)th off-diagonal

**Symmetric Perturbations:**
- Apply perturbations symmetrically (e.g., positions 1 and n-1)
- This preserves certain symmetries in the characteristic polynomial

**Double Perturbations:**
- Two adjacent or non-adjacent changes
- Look for patterns in placement

### Step 3: Numerical Testing

For each pattern, test with **concrete numeric values**:

```julia
# Example: Test pattern [b, d, b, b] with a=1, b=2, d=3
mat = [1 2 0 0 0
       2 1 3 0 0
       0 3 1 2 0
       0 0 2 1 2
       0 0 0 2 1]

eigs = eigvals(mat)  # Compute numerically
# Result: [-3.123106, -1.0, 1.0, 3.0, 5.123106]
```

### Step 4: Pattern Recognition

Look for eigenvalues that match **simple formulas**:

**Level 1: Exact algebraic combinations**
- λ = a (exactly)
- λ = a ± b
- λ = a ± 2b
- λ = a ± d

**Level 2: Square root expressions**
- λ = a ± √(αb² + βd²) for small integers α, β
- Common patterns:
  - √(b² + d²)
  - √(2b² + d²)
  - √(b² + 2d²)
  - √(2b² + 2d²) = √(2(b² + d²))

**Level 3: More complex expressions**
- λ = a ± √(b² + d² + 2bd)
- λ = a ± √((b±d)²) = a ± |b±d|

### Step 5: Formula Verification

Once a pattern is suspected:

1. **Test with multiple parameter sets**:
   ```julia
   test_params = [
       (a=1, b=2, d=3),
       (a=0, b=1, d=1),
       (a=5, b=3, d=7),
       (a=-1, b=4, d=2)
   ]
   ```

2. **Verify all eigenvalues** are accounted for

3. **Check edge cases**: b=d, b=0, d=0, etc.

### Step 6: Symbolic Verification

Implement detection and verify symbolically:

```julia
@variables a b d
mat = construct_pattern(a, b, d)
vals = symbolic_eigenvalues(mat)

# Check if formulas match
expected = [a - √(2b² + d²), a - b, a, a + b, a + √(2b² + d²)]
verify_match(vals, expected)
```

## Discovered Patterns

### Pattern 1: [b, d, b, b] in 5×5 (Position 1)

**Matrix:**
```
[a  b  0  0  0]
[b  a  d  0  0]
[0  d  a  b  0]
[0  0  b  a  b]
[0  0  0  b  a]
```

**Eigenvalues:**
- λ₁ = a - √(2b² + d²)
- λ₂ = a - b
- λ₃ = a
- λ₄ = a + b
- λ₅ = a + √(2b² + d²)

**Status:** ✅ Implemented and tested

### Pattern 2: [b, b, d, b] in 5×5 (Position 2)

**Matrix:**
```
[a  b  0  0  0]
[b  a  b  0  0]
[0  b  a  d  0]
[0  0  d  a  b]
[0  0  0  b  a]
```

**Eigenvalues:** Same as Pattern 1!

**Status:** ✅ Implemented and tested

**Key Insight:** Perturbation positions 1 and 2 give identical eigenvalues. This suggests an underlying symmetry or equivalence class.

### Pattern 3: [b, d, b, b, d, b] in 6×6 (Symmetric)

**Matrix:**
```
[a  b  0  0  0  0  0]
[b  a  d  0  0  0  0]
[0  d  a  b  0  0  0]
[0  0  b  a  d  0  0]
[0  0  0  d  a  b  0]
[0  0  0  0  b  a  b]
[0  0  0  0  0  b  a]
```

**Eigenvalues (partial):**
- λ = a (exactly!)
- λ = a ± √(b² + d²)
- Plus 4 more eigenvalues (currently unknown closed-form)

**Status:** ⚠️ Partially characterized - needs more investigation

## Patterns That FAILED

Not all perturbations work. These patterns do NOT have simple closed-forms:

### Failed Pattern: [d, b, b, b] (Perturbation at start)

**Matrix:**
```
[a  d  0  0  0]
[d  a  b  0  0]
[0  b  a  b  0]
[0  0  b  a  b]
[0  0  0  b  a]
```

**Result:** Only λ = a has closed-form; others are complex

**Lesson:** Perturbations at the boundaries (position 0 or n-1) tend to break the structure

### Failed Pattern: [b, d, d, b, b] (Adjacent double perturbation)

**Result:** Only λ = a has closed-form

**Lesson:** Adjacent perturbations seem to create interference that prevents closed-form solutions

## Key Insights and Heuristics

### What Makes a Pattern Work?

1. **Interior Placement**: Perturbations in the "middle" (positions 1-3 for 5×5) work better than edges

2. **Symmetry Preservation**: Patterns that preserve certain symmetries are more likely to have closed-forms

3. **Single Perturbation**: One deviation from constant is easier than multiple

4. **Specific Positions**: Not all interior positions work equally - positions 1 and 2 in 5×5 work, but not position 0 or 4

### Discovery Heuristics

**Priority 1: Test these first**
- Single perturbation at position 1 or 2
- Symmetric perturbations (mirror positions)
- Constant tridiagonal (known to work)

**Priority 2: Worth investigating**
- Single perturbation at other interior positions
- Patterns with perturbations separated by ≥2 positions

**Priority 3: Unlikely to work**
- Perturbations at boundaries (position 0 or n-1)
- Adjacent multiple perturbations
- More than 2 unique off-diagonal values

## Implementation Strategy

### Current Implementation

The package now includes:

1. Detection of [b, d, b, b] pattern (position 1)
2. Detection of [b, b, d, b] pattern (position 2)
3. Closed-form formula: a ± √(2b² + d²), a ± b, a

### Code Structure

```julia
function _detect_special_5x5_tridiagonal(mat)
    # 1. Check size
    # 2. Check tridiagonal structure
    # 3. Check constant diagonal
    # 4. Check symmetry
    # 5. Check off-diagonal pattern
    # 6. Return (a, b, d) or nothing
end
```

### Testing Approach

For each pattern:
1. Symbolic test with variables a, b, d
2. Numeric verification with concrete values
3. Edge case testing (b=d, b=0, etc.)
4. Comparison with LinearAlgebra.eigen()

## Future Directions

### Immediate Next Steps

1. **Investigate 7×7 patterns**: Test [b, d, b, b, b, b] and [b, b, d, b, b, b]
   - Hypothesis: May have similar structure with a ± √(αb² + d²)

2. **Complete 6×6 symmetric pattern**: Find closed-forms for the remaining 4 eigenvalues

3. **Test symmetric double perturbations**: [b, d, b, b, d, b] in various sizes

### Theoretical Questions

1. **Why does [b, d, b, b] work?**
   - Can we prove the eigenvalue formula analytically?
   - Is there a connection to orthogonal polynomials?

2. **What's the pattern in α coefficient?**
   - 5×5 has √(2b² + d²)
   - Does 7×7 have √(3b² + d²)?
   - General formula for (2k+1)×(2k+1)?

3. **Equivalence classes:**
   - [b, d, b, b] ≡ [b, b, d, b]
   - Are there other equivalent patterns?

### Long-term Research

1. **General theory**: Develop criteria for "eigenvalue-solvable" tridiagonal matrices

2. **Connection to known mathematics**:
   - Sturm-Liouville theory
   - Orthogonal polynomials (Chebyshev, Legendre, etc.)
   - Quantum mechanics (tight-binding models)

3. **Extension to other structures**:
   - Pentadiagonal matrices
   - Block-tridiagonal matrices
   - Non-symmetric patterns

## Tools and Scripts

### exploration_patterns.jl

Systematic testing script that:
- Tests constant tridiagonal patterns (sizes 3-7)
- Tests single perturbations at all positions
- Tests symmetric perturbations
- Performs pattern recognition on eigenvalues
- Outputs matches for simple formulas

**Usage:**
```bash
julia --project=. explore_patterns.jl
```

### Pattern Detection Function

```julia
function detect_pattern_type(off_diagonals)
    unique_vals = unique(off_diagonals)
    
    if length(unique_vals) == 1
        return :constant
    elseif length(unique_vals) == 2
        # Analyze positions of perturbation
        counts = [count(==(v), off_diagonals) for v in unique_vals]
        rare_val_index = argmin(counts)
        rare_val = unique_vals[rare_val_index]
        positions = findall(==(rare_val), off_diagonals)
        
        if length(positions) == 1
            return (:single_perturbation, positions[1]-1)
        end
    end
    
    return :unknown
end
```

## Conclusion

We've developed a robust methodology for discovering closed-form eigenvalue patterns:

1. **Start from known patterns** (constant tridiagonal)
2. **Apply systematic perturbations** (single, symmetric, double)
3. **Test numerically** with concrete values
4. **Recognize patterns** in eigenvalues (exact values, square roots)
5. **Verify symbolically** and implement detection
6. **Test thoroughly** with edge cases

**Key Success:** Discovered that 5×5 patterns [b, d, b, b] and [b, b, d, b] have identical closed-form eigenvalues, and successfully implemented detection for both.

**Next Challenge:** Extend to 7×7 and understand the theoretical basis for these patterns.

The methodology is general and can be applied to discover other special matrix structures with closed-form eigenvalues.
