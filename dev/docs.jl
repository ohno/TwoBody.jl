# Please run `include("./dev/docs.jl")` on RELP.

run(`julia --project=docs/ -e 'using Pkg; Pkg.activate("./"); cd("docs"); include("make.jl")'`)
