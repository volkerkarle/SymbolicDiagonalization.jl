using Documenter

# Make the package available when building docs without needing installation.
push!(LOAD_PATH, "..")
using SymbolicDiagonalization

makedocs(;
    sitename = "SymbolicDiagonalization",
    modules = [SymbolicDiagonalization],
    pages = [
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "API Reference" => "api_reference.md",
        "Pattern Library" => "pattern_library.md",
        "Group Theory Examples" => "group_theory_examples.md",
        "Implementation Details" => "implementation.md",
        "Mathematical Background" => "mathematical_background.md",
        "Contributing" => "contributing.md",
    ],
    checkdocs = :none,  # Don't error on missing docstrings for internal functions
    repo = "https://github.com/volkerkarle/SymbolicDiagonalization.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://volkerkarle.github.io/SymbolicDiagonalization.jl",
        edit_link = "main",
        assets = String[],
        repolink = "https://github.com/volkerkarle/SymbolicDiagonalization.jl",
    ),
)

# Deploy documentation to GitHub Pages
deploydocs(;
    repo = "github.com/volkerkarle/SymbolicDiagonalization.jl",
    devbranch = "main",
    push_preview = true,
)
