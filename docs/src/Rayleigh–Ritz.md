```@meta
CurrentModule = TwoBody
```

# Rayleigh–Ritz Method

This method is one of the variational method. It solves the generalized eigenvalue problem,
```math
\pmb{H} \pmb{c} = E \pmb{S} \pmb{c}.
```
The Hamiltonian matrix is defined as ``H_{ij} = \langle \phi_{i} | \hat{H} | \phi_{j} \rangle`` and the overlap matrix is defined as ``S_{ij} = \langle \phi_{i} | \phi_{j} \rangle``. The eigenvector ``\pmb{c}`` is the column of the optimal coefficients $c_i$ for the linear combination,
```math
\psi(r) = \sum_i c_i \phi_i(r),
```
to minmize the expectation value of the energy,
```math
E = \frac{\langle\psi|\hat{H}|\psi\rangle}{\langle\psi|\psi\rangle}.
```
Note that the nonlinear parameters (e.g., exponents of the Gaussian basis functions) are not optimized. The expectation by trial wavefunction is the upper bound for the exact energy.

## Examples

Run the following code before each use.

```@example example
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

Solve the eigenvalue problem. You should find
```math
E_{n=1} = -0.499278~E_\mathrm{h},
```
which is amazingly good for only four basis functions according to [Thijssen(2007)](https://doi.org/10.1017/CBO9781139171397). The exact ground-state energy is ``-0.5~E_\mathrm{h}``.

```@repl example
solve(H, BS)
```

Here is a comprehensive example including calculations up to excited states.

```@example example
# solve
using TwoBody
H = Hamiltonian(NonRelativisticKinetic(1,1), CoulombPotential(-1))
BS = GeometricBasisSet(SimpleGaussianBasis, 0.1, 80.0, 20)
res = solve(H, BS, info=0)

# benchmark
import Antique
HA = Antique.HydrogenAtom(Z=1, Eₕ=1.0, a₀=1.0, mₑ=1.0, ℏ=1.0)

# plot
using CairoMakie
fig = Figure(
  size = (840,600),
  fontsize = 11.5,
  backgroundcolor = :transparent
)
for n in 1:4
  axis = Axis(
    fig[div(n-1,2)+1,rem(n-1,2)+1],
    xlabel = L"$r~/~a_0$",
    ylabel = L"$4\pi r^2|\psi(r)|^2~ /~{a_0}^{-1}$",
    xlabelsize = 16.5,
    ylabelsize = 16.5,
    limits=(
      0, [5, 15, 30, 50][n],
      0, [0.6, 0.2, 0.11, 0.07][n],
    )
  )
  lines!(axis, 0..50, r -> 4π * r^2 * abs(res.ψ[n](r))^2, label="TwoBody.jl")
  lines!(axis, 0..50, r -> 4π * r^2 * abs(Antique.ψ(HA,r,0,0,n=n))^2, label="Antique.jl", color=:black, linestyle=:dash)
  axislegend(axis, "n = $n", position=:rt, framevisible=false)
end
fig
```

## Solver

```@docs; canonical=false
solve(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4)
solve(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=4)
solve(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4)
optimize(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)
optimize(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=1, progress=true, optimizer=Optim.NelderMead(), options...)
optimize(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)
```

## Basis Set

```@docs; canonical=false
TwoBody.BasisSet
TwoBody.GeometricBasisSet
TwoBody.geometric(r₁, rₙ, n::Int; nₘₐₓ::Int=n, nₘᵢₙ::Int=1)
```

## Basis Functions

```@docs; canonical=false
TwoBody.SimpleGaussianBasis
TwoBody.ContractedBasis
```

## Matrix

```@docs; canonical=false
matrix(basisset::BasisSet)
matrix(hamiltonian::Hamiltonian, basisset::BasisSet)
```

## Matrix Elements

```@docs; canonical=false
element(SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::Hamiltonian, B1::Basis, B2::Basis)
element(o::RestEnergy, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::Laplacian, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::NonRelativisticKinetic, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::ConstantPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::LinearPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::CoulombPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::PowerLawPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
element(o::GaussianPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
```
