# Please run `include("./dev/revice.jl")` on RELP

# import Pkg; Pkg.add("Revise")
using Revise

dir = dirname(@__FILE__) * "/../"
cd(dir)
@show pwd()
import Pkg; Pkg.activate(".")
using TwoBody
