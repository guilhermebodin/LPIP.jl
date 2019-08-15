using Documenter, LPIP

makedocs(
    modules = [LPIP],
    doctest  = false,
    clean = true,
    # assets = ["assets/favicon.ico"],
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "LPIP.jl",
    authors = "Guilherme Bodin, and contributors",
    pages = [
        "Home" => "index.md",
        "manual.md",
        "examples.md",
        "reference.md"
    ]
)

deploydocs(
    repo = "github.com/guilhermebodin/LPIP.jl.git",
)