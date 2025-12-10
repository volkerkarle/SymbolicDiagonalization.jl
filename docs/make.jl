using Documenter

# Make the package available when building docs without needing installation.
push!(LOAD_PATH, "..")
using SymbolicDiagonalization

makedocs(;
    sitename = "SymbolicDiagonalization",
    modules = [SymbolicDiagonalization],
    pages = [
        "Home" => "index.md",
    ],
)

# No deploy step: build locally via `julia --project=docs docs/make.jl`.
