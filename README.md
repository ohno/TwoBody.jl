# TwoBody.jl [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ohno.github.io/TwoBody.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ohno.github.io/TwoBody.jl/dev) [![Build Status](https://github.com/ohno/TwoBody.jl/workflows/CI/badge.svg)](https://github.com/ohno/TwoBody.jl/actions)

TwoBody.jl: a Julia package for quantum mechanical two-body problems

## Documentation 

https://ohno.github.io/TwoBody.jl/dev/

## Dependency

```mermaid
---
config:
  layout: elk
  theme: mc
---
flowchart TD
  A["Hamiltonians.jl"]
  B["Basis.jl"]
  C["Rayleigh-Ritz.jl"]
  D["GEM.jl"]
  E["FiniteDifferenceMatrices.jl"]
  F["FDM.jl"]
  G["VMC.jl"]
  Z["TwoBody.jl"]
  A --> B & F & G
  B --> C
  B --> D
  E --> F
  C & D & F & G --> Z
```

## Developer's Guide

There are several tools for developers.

```sh
git clone https://github.com/ohno/TwoBody.jl.git
cd TwoBody.jl
julia
julia> include("dev/revice.jl")
julia> include("dev/test.jl")
julia> include("dev/docs.jl")
```
