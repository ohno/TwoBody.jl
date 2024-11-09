```@meta
CurrentModule = TwoBody
```

# Rayleigh–Ritz method

The solver finds the best $c_i$ for the trial wavefunction
```math
\varPsi(r) = \sum_i c_i \phi_i(r),
```
to minmize the expectation value of the energy
```math
E = \frac{\langle\psi|\hat{H}|\psi\rangle}{\langle\psi|\psi\rangle}.
```
Here, the energy by trial wavefunction is the upper bound for the exact energy.

## Exaples

```@setup example
using TwoBody
```

Define the [Hamiltoninan](@ref Hamiltonian). This is an example for the non-relativistic Hamiltonian of hydrogen atom in atomic units:
```math
\hat{H} = 
- \frac{1}{2} \nabla^2
- \frac{1}{r}
```
```@repl example
hamiltonian = Hamiltonian(
  NonRelativisticKinetic(ℏ = 1 , m = 1),
  CoulombPotential(coefficient = -1),
)
```

Define the basis set:
```math
\begin{aligned}
  \phi_1(r) &= \mathrm{e}^{-13.00773 r^2}, \\
  \phi_2(r) &= \mathrm{e}^{-1.962079 r^2}, \\
  \phi_3(r) &= \mathrm{e}^{-0.444529 r^2}, \\
  \phi_4(r) &= \mathrm{e}^{-0.1219492 r^2}.
\end{aligned}
```
```@repl example
basisset = BasisSet(
  SimpleGaussianBasis(13.00773),
  SimpleGaussianBasis(1.962079),
  SimpleGaussianBasis(0.444529),
  SimpleGaussianBasis(0.1219492),
)
```

You should find
```math
E_{n=1} = -0.499278~E_\mathrm{h},
```
which is amazingly good for only four basis functions according to [Thijssen(2007)](https://doi.org/10.1017/CBO9781139171397). The exact ground-state energy is ``-0.5~E_\mathrm{h}``.

```@repl example
solve(hamiltonian, basisset)
```

## Solver

```@docs; canonical=false
TwoBody.solve
```

## Basis Set

```@docs; canonical=false
TwoBody.BasisSet
```

## Basis Functions

```@docs; canonical=false
TwoBody.SimpleGaussianBasis
TwoBody.GaussianBasis
TwoBody.ContractedBasis
```

## Matrix Elements

```@docs; canonical=false
TwoBody.element
```
