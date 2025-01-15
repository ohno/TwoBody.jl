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
```@example example
H = Hamiltonian(
  NonRelativisticKinetic(ℏ = 1 , m = 1),
  CoulombPotential(coefficient = -1),
)
nothing # hide
```

Define the basis set:
```math
\begin{aligned}
  \phi_1(r) &= \exp(-13.00773 ~r^2), \\
  \phi_2(r) &= \exp(-1.962079 ~r^2), \\
  \phi_3(r) &= \exp(-0.444529 ~r^2), \\
  \phi_4(r) &= \exp(-0.1219492 ~r^2).
\end{aligned}
```
```@example example
BS = BasisSet(
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

```@repl example
solve(H, BS)
```

## Solver

```@docs; canonical=false
TwoBody.solve
TwoBody.optimize
```

## Basis Set

```@docs; canonical=false
TwoBody.BasisSet
TwoBody.GeometricBasisSet
TwoBody.geometric
```

## Basis Functions

```@docs; canonical=false
TwoBody.SimpleGaussianBasis
TwoBody.ContractedBasis
```

## Matrix Elements

```@docs; canonical=false
TwoBody.element
```
