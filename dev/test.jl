# Please run `include("./dev/test.jl")` on RELP.

dir = dirname(@__FILE__) * "/../"
cd(dir)
@show pwd()

using Pkg
Pkg.activate(dir)
Pkg.test()
