using Fractal
using Documenter

DocMeta.setdocmeta!(Fractal, :DocTestSetup, :(using Fractal); recursive=true)

makedocs(;
    modules=[Fractal],
    authors="Rochas-Martin",
    repo="https://github.com/MartiRoc/Fractal.jl/blob/{commit}{path}#{line}",
    sitename="Fractal.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MartiRoc.github.io/Fractal.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MartiRoc/Fractal.jl",
    devbranch="master",
)
