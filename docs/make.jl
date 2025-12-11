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
        "Implementation Details" => "implementation.md",
        "Mathematical Background" => "mathematical_background.md",
        "Contributing" => "contributing.md",
    ],
    checkdocs = :none,  # Don't error on missing docstrings for internal functions
)

# No deploy step: build locally via `julia --project=docs docs/make.jl`.
