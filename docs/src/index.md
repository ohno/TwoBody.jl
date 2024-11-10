```@meta
CurrentModule = TwoBody
```

# TwoBody.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ohno.github.io/TwoBody.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ohno.github.io/TwoBody.jl/dev) [![Build Status](https://github.com/ohno/TwoBody.jl/workflows/CI/badge.svg)](https://github.com/ohno/TwoBody.jl/actions)

[TwoBody.jl](https://github.com/ohno/TwoBody.jl): a Julia package for two-body problems

## Install

Run the following code on the REPL or Jupyter Notebook to install this package.

```julia
import Pkg; Pkg.add(url="https://github.com/ohno/TwoBody.jl.git")
```

## Usage

Run the following code before each use.
```@example index
using TwoBody
```

Define the [Hamiltoninan](@ref Hamiltonian). This is an example for the non-relativistic Hamiltonian of hydrogen atom in atomic units:
```math
\hat{H} = 
- \frac{1}{2} \nabla^2
- \frac{1}{r}
```
```@example index
hamiltonian = Hamiltonian(
  NonRelativisticKinetic(ℏ = 1 , m = 1),
  CoulombPotential(coefficient = -1),
)
nothing # hide
```

The usage depends on the method. Define the basis set for the [Rayleigh–Ritz method](@ref Rayleigh–Ritz-method):
```math
\begin{aligned}
  \phi_1(r) &= \exp(-13.00773 ~r^2), \\
  \phi_2(r) &= \exp(-1.962079 ~r^2), \\
  \phi_3(r) &= \exp(-0.444529 ~r^2), \\
  \phi_4(r) &= \exp(-0.1219492 ~r^2).
\end{aligned}
```
```@example index
basisset = BasisSet(
  SimpleGaussianBasis(13.00773),
  SimpleGaussianBasis(1.962079),
  SimpleGaussianBasis(0.444529),
  SimpleGaussianBasis(0.1219492),
)
nothing # hide
```

You should find
```math
E_{n=1} = -0.499278~E_\mathrm{h},
```
which is amazingly good for only four basis functions according to [Thijssen(2007)](https://doi.org/10.1017/CBO9781139171397). The exact ground-state energy is ``-0.5~E_\mathrm{h}``.

```@repl index
solve(hamiltonian, basisset)
```

## API reference

```@index
```
