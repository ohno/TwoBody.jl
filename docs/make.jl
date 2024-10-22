using TwoBody
using Documenter

DocMeta.setdocmeta!(TwoBody, :DocTestSetup, :(using TwoBody); recursive=true)

makedocs(;
    modules=[TwoBody],
    authors="Shuhei Ohno",
    repo="https://github.com/ohno/TwoBody.jl/blob/{commit}{path}#{line}",
    sitename="TwoBody.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ohno.github.io/TwoBody.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ohno/TwoBody.jl",
)
