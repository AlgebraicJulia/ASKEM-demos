using ASKEM
using Documenter

DocMeta.setdocmeta!(ASKEM, :DocTestSetup, :(using ASKEM); recursive=true)

makedocs(;
    modules=[ASKEM],
    authors="AlgebraicJulia team",
    repo="https://github.com/AlgebraicJulia/ASKEM.jl/blob/{commit}{path}#{line}",
    sitename="ASKEM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://AlgebraicJulia.github.io/ASKEM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AlgebraicJulia/ASKEM.jl",
    devbranch="main",
)
