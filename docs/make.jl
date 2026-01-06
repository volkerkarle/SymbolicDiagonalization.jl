using Documenter
using Latexify
using Symbolics

push!(LOAD_PATH, "..")
using SymbolicDiagonalization

# ============================================================================
# LaTeX Display Wrapper for Documentation
# ============================================================================

"""
    LaTeX(expr)

Wrapper type that displays symbolic expressions as rendered LaTeX in documentation.
Use with Documenter.jl's `@example` blocks for beautiful mathematical output.

# Example
```julia
@variables a b
vals = eigvals([a b; b a])
LaTeX(vals)  # Renders as formatted LaTeX
```
"""
struct LaTeX{T}
    expr::T
end

# Display as LaTeX math block in HTML (Documenter uses this)
function Base.show(io::IO, ::MIME"text/html", r::LaTeX)
    if r.expr isa AbstractVector
        # For vectors (like eigenvalue lists), display each on its own line
        print(io, "<p><strong>Eigenvalues:</strong></p>\n")
        print(io, "\\[\n\\begin{aligned}\n")
        for (i, v) in enumerate(r.expr)
            latex_str = latexify(v, env=:raw)
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
        latex_str = latexify(r.expr, env=:raw)
        print(io, latex_str)
        print(io, "\n\\]")
    else
        # Single expression
        latex_str = latexify(r.expr, env=:raw)
        print(io, "\\[", latex_str, "\\]")
    end
end

# Plain text fallback
function Base.show(io::IO, ::MIME"text/plain", r::LaTeX)
    show(io, MIME("text/plain"), r.expr)
end

# Also allow regular show
function Base.show(io::IO, r::LaTeX)
    show(io, r.expr)
end

makedocs(;
    sitename = "SymbolicDiagonalization",
    modules = [SymbolicDiagonalization],
    pages = [
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "Pattern Library" => "pattern_library.md",
        "Mathematical Background" => "mathematical_background.md",
    ],
    checkdocs = :none,
    repo = "https://github.com/volkerkarle/SymbolicDiagonalization.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://volkerkarle.github.io/SymbolicDiagonalization.jl",
        edit_link = "main",
        assets = String[],
        repolink = "https://github.com/volkerkarle/SymbolicDiagonalization.jl",
        mathengine = Documenter.MathJax3(),
    ),
)

deploydocs(;
    repo = "github.com/volkerkarle/SymbolicDiagonalization.jl",
    devbranch = "main",
    push_preview = true,
)
