# ============================================================================
# LaTeX Display Support
# ============================================================================
#
# Provides a LaTeX wrapper type for displaying symbolic expressions as
# beautifully formatted LaTeX/MathJax output in notebooks, documentation,
# and other environments that support HTML rendering.
#
# Requires Latexify.jl to be loaded (weak dependency).
# ============================================================================

"""
    LaTeX(expr)

Wrapper type that displays symbolic expressions as rendered LaTeX.

When displayed in an environment that supports HTML (Jupyter notebooks,
Documenter.jl, Pluto, etc.), the expression will be rendered as a
beautifully formatted mathematical equation using MathJax.

# Examples

```julia
using SymbolicDiagonalization, Symbolics, LinearAlgebra

@variables a b c
M = [a b; b c]

# Display eigenvalues as LaTeX
LaTeX(eigvals(M))

# Display matrix as LaTeX
LaTeX(M)

# Display single expression
LaTeX(a^2 + sqrt(b))
```

# Requirements

Requires `Latexify.jl` to be installed. The package will attempt to load
it automatically when `LaTeX` is first used.
"""
struct LaTeX{T}
    expr::T
end

# Flag to track if Latexify is available
const _latexify_loaded = Ref(false)
const _latexify_func = Ref{Any}(nothing)

"""
    _ensure_latexify()

Attempt to load Latexify.jl if not already loaded.
Returns true if successful, false otherwise.
"""
function _ensure_latexify()
    if _latexify_loaded[]
        return true
    end
    
    try
        # Try to load Latexify
        @eval Main begin
            using Latexify
        end
        _latexify_func[] = Main.Latexify.latexify
        _latexify_loaded[] = true
        return true
    catch e
        return false
    end
end

"""
    _to_latex(expr; env=:raw)

Convert an expression to LaTeX string using Latexify.
Falls back to string representation if Latexify is not available.
"""
function _to_latex(expr; env=:raw)
    if _ensure_latexify()
        try
            return string(_latexify_func[](expr, env=env))
        catch
            return string(expr)
        end
    else
        return string(expr)
    end
end

# HTML display for LaTeX wrapper (used by notebooks, Documenter, etc.)
function Base.show(io::IO, ::MIME"text/html", r::LaTeX)
    if !_ensure_latexify()
        # Fallback: display as code block
        print(io, "<pre>", r.expr, "</pre>")
        print(io, "<p><em>Note: Install Latexify.jl for LaTeX rendering</em></p>")
        return
    end
    
    if r.expr isa AbstractVector
        # For vectors (like eigenvalue lists), display each on its own line
        print(io, "<p><strong>Eigenvalues:</strong></p>\n")
        print(io, "\\[\n\\begin{aligned}\n")
        for (i, v) in enumerate(r.expr)
            latex_str = _to_latex(v)
            print(io, "\\lambda_{$i} &= ", latex_str)
            if i < length(r.expr)
                print(io, " \\\\[0.3em]\n")
            else
                print(io, "\n")
            end
        end
        print(io, "\\end{aligned}\n\\]")
    elseif r.expr isa AbstractMatrix
        # For matrices, use array environment
        print(io, "\\[\n")
        latex_str = _to_latex(r.expr)
        print(io, latex_str)
        print(io, "\n\\]")
    else
        # Single expression
        latex_str = _to_latex(r.expr)
        print(io, "\\[", latex_str, "\\]")
    end
end

# Plain text fallback
function Base.show(io::IO, ::MIME"text/plain", r::LaTeX)
    show(io, MIME("text/plain"), r.expr)
end

# Regular show
function Base.show(io::IO, r::LaTeX)
    show(io, r.expr)
end

# LaTeX MIME for environments that support it directly
function Base.show(io::IO, ::MIME"text/latex", r::LaTeX)
    if !_ensure_latexify()
        print(io, string(r.expr))
        return
    end
    
    if r.expr isa AbstractVector
        print(io, "\\begin{aligned}\n")
        for (i, v) in enumerate(r.expr)
            latex_str = _to_latex(v)
            print(io, "\\lambda_{$i} &= ", latex_str)
            if i < length(r.expr)
                print(io, " \\\\\n")
            else
                print(io, "\n")
            end
        end
        print(io, "\\end{aligned}")
    elseif r.expr isa AbstractMatrix
        print(io, _to_latex(r.expr))
    else
        print(io, _to_latex(r.expr))
    end
end
