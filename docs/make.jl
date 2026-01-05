using Documenter

push!(LOAD_PATH, "..")
using SymbolicDiagonalization

makedocs(;
    sitename = "SymbolicDiagonalization",
    modules = [SymbolicDiagonalization],
    pages = ["Documentation" => "index.md"],
    checkdocs = :none,
    repo = "https://github.com/volkerkarle/SymbolicDiagonalization.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://volkerkarle.github.io/SymbolicDiagonalization.jl",
        edit_link = "main",
        assets = String[],
        repolink = "https://github.com/volkerkarle/SymbolicDiagonalization.jl",
    ),
)

deploydocs(;
    repo = "github.com/volkerkarle/SymbolicDiagonalization.jl",
    devbranch = "main",
    push_preview = true,
)
