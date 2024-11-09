# This is just a reminder. You don't use it.

# import Pkg; Pkg.add("PkgTemplates")
using PkgTemplates

t = Template(;
  dir = pwd(),
  user = "ohno",
  authors = ["Shuhei Ohno"],
  julia = v"1.7",
  plugins = [
    ProjectFile(; version = v"0.0.1"),
    License(; name = "MIT"),
    Git(; manifest = false),
    GitHubActions(; extra_versions = ["1.10"]),
    Documenter{GitHubActions}(),
    Tests(; project = true),
    Readme(; inline_badges = true),
  ],
)

generate("TwoBody.jl", t)
