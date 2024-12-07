export solve, optimize

using LinearAlgebra
using Optim
using Printf
using SpecialFunctions
using Subscripts

@doc raw"""
`solve(hamiltonian::Hamiltonian, basisset::BasisSet)`

This function returns the eigenvalues ``E``  and eigenvectors ``\pmb{c}`` for
```math
\pmb{H} \pmb{c} = E \pmb{S} \pmb{c}.
```
The Hamiltonian matrix is defined as ``H_{ij} = \langle \phi_{i} | \hat{H} | \phi_{j} \rangle``. The overlap matrix is defined as ``S_{ij} = \langle \phi_{i} | \phi_{j} \rangle``.
"""
function solve(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4)

  # Initialization
  nₘₐₓ = length(basisset.basis)

  # Matrix Elements
  S = [element(basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
  H = [element(hamiltonian, basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]

  # Calculations
  E, C = eigen(Hermitian(H), Hermitian(S))

  # Output
  if 0 < info
    println("\nn \tbasis function φₙ")
    for n in 1:nₘₐₓ
    println("$n\t", basisset.basis[n])
    end
    println("\nn \twavefuntion ψₙ")
    for n in 1:min(nₘₐₓ, info)
      print("$n\t")
      for i in 1:nₘₐₓ
        @printf(" %s %.6f φ%s", C[i,n]<0 ? "-" : "+", abs(C[i,n]), sub("$i"))
      end
      println()
    end
    println("\nn \tnorm, <ψ|ψ> = c' * S * c")
    for n in 1:min(nₘₐₓ, info)
      println("$n\t", C[:,n]' * S * C[:,n])
    end
    println("\nn \teigenvalue, E")
    for n in 1:min(nₘₐₓ, info)
      println("$n\t", E[n]) 
    end
    if !isempty(perturbation.terms)
      println("\nn \tperturbation")
      M = [element(perturbation, basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
      for n in 1:min(nₘₐₓ, info)
      println("$n\t", C[:,n]' * M * C[:,n])
      end
      println("\nn \teigenvalue + perturbation")
      for n in 1:min(nₘₐₓ, info)
          println("$n\t", E[n] + C[:,n]' * M * C[:,n])
      end
    end
    println("\nn \texpectation value of the Hamiltonian, <ψ|H|ψ> = c' * H * c")
    for n in 1:min(nₘₐₓ, info)
      println("$n\t", C[:,n]' * H * C[:,n])
    end
    for term in [hamiltonian.terms..., perturbation.terms...]
      println("\nn \texpectation value of $(term)")
      M = [element(term, basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
      for n in 1:min(nₘₐₓ, info)
        println("$n\t", C[:,n]' * M * C[:,n])
      end
    end
    println()
    return (hamiltonian=hamiltonian, basisset=basisset, E=E, C=C, S=S, H=H)
  end
  return E
end

@doc raw"""
`solve(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4)`
This function is a wrapper for `solve(hamiltonian::Hamiltonian, basisset::BasisSet, ...)`.
"""
function solve(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4)
  if 0 < info
    println("\n  \tgeometric progression")
    println("type \t$(basisset.basistype)")
    println("range\tr", sub("$(basisset.nₘᵢₙ)"), " - r", sub("$(basisset.nₘₐₓ)"))
    println("r", sub("$(basisset.nₘᵢₙ)"), " \t", basisset.r₁)
    println("r", sub("$(basisset.n)"  ), " \t", basisset.rₙ)
  end
  solve(hamiltonian, BasisSet(basisset.basis...); perturbation=perturbation, info=info)
end

@doc raw"""
`optimize(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4, optimizer=Optim.NelderMead())`
This function minimizes the energy by optimizing $r_1$ and $r_n$ using Optim.jl.
```math
\frac{\partial E}{\partial r_1} = \frac{\partial E}{\partial r_n} = 0
```
"""
function optimize(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4, optimizer=Optim.NelderMead())
  res = Optim.optimize(x -> solve(hamiltonian, GeometricBasisSet(basisset.basistype, x..., basisset.n, nₘₐₓ=basisset.nₘₐₓ, nₘᵢₙ=basisset.nₘᵢₙ), perturbation=perturbation, info=0)[1], [basisset.r₁, basisset.rₙ], method=optimizer)
  solve(hamiltonian, GeometricBasisSet(basisset.basistype, res.minimizer..., basisset.n, nₘₐₓ=basisset.nₘₐₓ, nₘᵢₙ=basisset.nₘᵢₙ); perturbation=perturbation, info=info)
end

# SGEM

@doc raw"""
`element(SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  S_{ij}
   = \langle \phi_{i} | \phi_{j} \rangle
  &= \iiint
     \phi_{i}^*(r)
     \phi_{j}(r)
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \int_0^{2\pi} \mathrm{d}\varphi
     \int_0^\pi \sin\theta ~\mathrm{d}\theta
     \int_0^\infty r^{2} \mathrm{e}^{-(\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= 2\pi \times 2 \times \frac{1!!}{2^{2}} \sqrt{\frac{\pi}{a^{3}}} \\
  &= \underline{\left( \frac{\pi}{\alpha_i + \alpha_j} \right)^{3/2}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{2n} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{(2n-1)!!}{2^{n+1}} \sqrt{\frac{\pi}{a^{2n+1}}}
```
"""
function element(SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return (π/(SGB1.a+SGB2.a))^(3/2)
end

@doc raw"""
`element(o::RestEnergy, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | mc^2 | \phi_{j} \rangle
  &= mc^2 \langle \phi_{i} | \phi_{j} \rangle
  &= mc^2 \iiint
     \phi_{i}^*(r)
     \phi_{j}(r)
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= mc^2
     \int_0^{2\pi} \mathrm{d}\varphi
     \int_0^\pi \sin\theta ~\mathrm{d}\theta
     \int_0^\infty r^{2} \mathrm{e}^{-(\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= mc^2 \times 2\pi \times 2 \times \frac{1!!}{2^{2}} \sqrt{\frac{\pi}{a^{3}}} \\
  &= \underline{mc^2 \left( \frac{\pi}{\alpha_i + \alpha_j} \right)^{3/2}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{2n} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{(2n-1)!!}{2^{n+1}} \sqrt{\frac{\pi}{a^{2n+1}}}
```
"""
function element(o::RestEnergy, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.m * o.c^2 * (π/(SGB1.a+SGB2.a))^(3/2)
end

@doc raw"""
`element(o::NonRelativisticKinetic, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  T_{ij} = \langle \phi_{i} | \hat{T} | \phi_{j} \rangle
  &= \iiint
     \mathrm{e}^{-\alpha_i r^2}
     \left[ -\frac{\hbar^2}{2\mu} \nabla^2 \right]
     \mathrm{e}^{-\alpha_j r^2}
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= -\frac{\hbar^2}{2\mu} \iiint
     \mathrm{e}^{-\alpha_i r^2}
     \left[ \nabla^2 \right]
     \mathrm{e}^{-\alpha_j r^2}
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= -\frac{\hbar^2}{2\mu} \iiint
     \mathrm{e}^{-\alpha_i r^2}
     \left[ -6\alpha_j + 4\alpha_j^2 r^2 \right]
     \mathrm{e}^{-\alpha_j r^2}
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= -\frac{\hbar^2}{2\mu} \iint
     \sin\theta ~\mathrm{d}\theta \mathrm{d}\varphi
     \int
     \left[ -6\alpha_j + 4\alpha_j^2 r^2 \right]
     r^2 \mathrm{e}^{-(\alpha_i + \alpha_j) r^2}
     ~\mathrm{d}r \\
  &= -\frac{\hbar^2}{2\mu} \cdot 4\pi
     \left[
        -6\alpha_j   \mathrm{GGI}(2, \alpha_i + \alpha_j)
        +4\alpha_j^2 \mathrm{GGI}(4, \alpha_i + \alpha_j)
     \right]
     \\
  &= -\frac{\hbar^2}{2\mu} \cdot 4\pi
     \left[
        -6\alpha_j   \frac{\Gamma\left( \frac{3}{2} \right)}{2 (\alpha_i + \alpha_j)^{\frac{3}{2}}}
        +4\alpha_j^2 \frac{\Gamma\left( \frac{5}{2} \right)}{2 (\alpha_i + \alpha_j)^{\frac{5}{2}}}
     \right] \\
  &= -\frac{\hbar^2}{2\mu} \cdot 4\pi
     \left[
        -6\alpha_j   \frac{ \sqrt{\pi}/2}{2 (\alpha_i + \alpha_j)^{\frac{3}{2}}}
        +4\alpha_j^2 \frac{3\sqrt{\pi}/4}{2 (\alpha_i + \alpha_j)^{\frac{5}{2}}}
     \right] \\
  &= -\frac{\hbar^2}{2\mu} \cdot 4\pi
     \left[
        \frac{\alpha_j}{\alpha_i + \alpha_j} - 1
     \right]
     \cdot 6 \alpha_j \cdot \frac{\sqrt{\pi}/2}{2 (\alpha_i + \alpha_j)^{\frac{3}{2}}}
     \\
  &= -\frac{\hbar^2}{2\mu} \cdot 4\pi
     \left[
        - \frac{\alpha_i}{\alpha_i + \alpha_j}
     \right]
     \cdot 6 \alpha_j \cdot \frac{\sqrt{\pi}/2}{2 (\alpha_i + \alpha_j)^{\frac{3}{2}}}
     \\
  &= \underline{
        -\frac{\hbar^2}{2\mu}
        \cdot 6
        \cdot \frac{\alpha_i \alpha_j \pi^{\frac{3}{2}}}{(\alpha_i + \alpha_j)^{\frac{5}{2}}}
     }
\end{aligned}
```

or

```math
\begin{aligned}
  T_{ij} = \langle \phi_{i} | \hat{T} | \phi_{j} \rangle
  &= \iiint
     \mathrm{e}^{-\alpha_i r^2}
     \left[ -\frac{\hbar^2}{2\mu} \nabla^2 \right]
     \mathrm{e}^{-\alpha_j r^2}
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= -\frac{\hbar^2}{2\mu} \iiint
     \mathrm{e}^{-\alpha_i r^2}
     \nabla^2
     \mathrm{e}^{-\alpha_j r^2}
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \frac{\hbar^2}{2\mu} \iiint
     \left[ \nabla \mathrm{e}^{-\alpha_i r^2} \right]
     \left[ \nabla \mathrm{e}^{-\alpha_j r^2} \right]
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \frac{\hbar^2}{2\mu} \iiint
     \left[ -2 \alpha_i r \mathrm{e}^{-\alpha_i r^2} \right]
     \left[ -2 \alpha_j r \mathrm{e}^{-\alpha_j r^2} \right]
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \frac{\hbar^2}{2\mu} \cdot 4 \alpha_i \alpha_j \iiint
     \left[ r \mathrm{e}^{-\alpha_i r^2} \right]
     \left[ r \mathrm{e}^{-\alpha_j r^2} \right]
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \frac{\hbar^2}{2\mu}
     \cdot 4 \alpha_i \alpha_j
     \iint \sin\theta ~\mathrm{d}\theta \mathrm{d}\varphi
     \int r^4
     \mathrm{e}^{- (\alpha_i + \alpha_j) r^2}
     ~\mathrm{d}r \\
  &= \frac{\hbar^2}{2\mu}
     \cdot 4 \alpha_i \alpha_j
     \cdot 4 \pi
     \cdot \mathrm{GGI}(4, \alpha_i + \alpha_j) \\
  &= \frac{\hbar^2}{2\mu}
     \cdot 4 \alpha_i \alpha_j
     \cdot 4 \pi
     \cdot \frac{\Gamma\left( \frac{5}{2} \right)}{2 (\alpha_i + \alpha_j)^{\frac{5}{2}}} \\
  &= \frac{\hbar^2}{2\mu}
     \cdot 4 \alpha_i \alpha_j
     \cdot 4 \pi
     \cdot \frac{3\sqrt{\pi}/4}{2 (\alpha_i + \alpha_j)^{\frac{5}{2}}} \\
  &= \underline{
        \frac{\hbar^2}{2\mu}
        \cdot 6
        \cdot \frac{\alpha_i \alpha_j \pi^{\frac{3}{2}}}{(\alpha_i + \alpha_j)^{\frac{5}{2}}}
     }
\end{aligned}
```
"""
function element(o::NonRelativisticKinetic, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.ℏ^2/(2*o.m) * 6*π^(3/2)*SGB1.a*SGB2.a/(SGB1.a+SGB2.a)^(5/2)
end

@doc raw"""
`element(o::ConstantPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | c | \phi_{j} \rangle
  &= c \langle \phi_{i} | \phi_{j} \rangle
  &= c \iiint
     \phi_{i}^*(r)
     \phi_{j}(r)
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= c
     \int_0^{2\pi} \mathrm{d}\varphi
     \int_0^\pi \sin\theta ~\mathrm{d}\theta
     \int_0^\infty r^{2} \mathrm{e}^{-(\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= c \times 2\pi \times 2 \times \frac{1!!}{2^{2}} \sqrt{\frac{\pi}{a^{3}}} \\
  &= \underline{c \left( \frac{\pi}{\alpha_i + \alpha_j} \right)^{3/2}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{2n} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{(2n-1)!!}{2^{n+1}} \sqrt{\frac{\pi}{a^{2n+1}}}
```
"""
function element(o::ConstantPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.constant * (π/(SGB1.a+SGB2.a))^(3/2)
end

@doc raw"""
`element(o::LinearPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | r | \phi_{j} \rangle
  &= \iiint
     \phi_{i}^*(r)
     \times r \times
     \phi_{j}(r)
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \int_0^{2\pi} \mathrm{d}\varphi
     \int_0^\pi \sin\theta ~\mathrm{d}\theta
     \int_0^\infty r^3 \mathrm{e}^{-(\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= 2\pi \times 2 \times \frac{1!}{2 (\alpha_i + \alpha_j)^{2}} \\
  &= \underline{\frac{2\pi}{(\alpha_i + \alpha_j)^2}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{2n+1} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{n!}{2 a^{n+1}}
```
"""
function element(o::LinearPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * 2*π/(SGB1.a+SGB2.a)^2
end

@doc raw"""
`element(o::CoulombPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | \frac{1}{r} | \phi_{j} \rangle
  &= \iiint
    \phi_{i}^*(r)
    \times \frac{1}{r} \times
    \phi_{j}(r)
    ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \int_0^{2\pi} \mathrm{d}\varphi
    \int_0^\pi \sin\theta ~\mathrm{d}\theta
    \int_0^\infty r \mathrm{e}^{-(\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= 2\pi \times 2 \times \frac{0!}{2 (\alpha_i + \alpha_j)} \\
  &= \underline{\frac{2\pi}{\alpha_i + \alpha_j}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{2n+1} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{n!}{2 a^{n+1}}
```
"""
function element(o::CoulombPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * 2*π/(SGB1.a+SGB2.a)
end

@doc raw"""
`element(o::PowerLawPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | r^n | \phi_{j} \rangle
  &= \iiint
     \phi_{i}^*(r)
     \times r^n \times
     \phi_{j}(r)
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \int_0^{2\pi} \mathrm{d}\varphi
     \int_0^\pi \sin\theta ~\mathrm{d}\theta
     \int_0^\infty r^{n+2} \mathrm{e}^{-(\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= 2\pi \times 2 \times \frac{\Gamma\left( \frac{n+3}{2} \right)}{2 (\alpha_i + \alpha_j)^{\frac{n+3}{2}}} \\
  &= \underline{2\pi\frac{\Gamma\left( \frac{n+3}{2} \right)}{(\alpha_i + \alpha_j)^{\frac{n+3}{2}}}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{n} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{\Gamma\left( \frac{n+1}{2} \right)}{2 a^{\frac{n+1}{2}}}
```
"""
function element(o::PowerLawPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * 2 * π * gamma((o.exponent+3)/2) / (SGB1.a+SGB2.a)^((o.exponent+3)/2)
end

@doc raw"""
`element(o::GaussianPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | \exp(-br^2) | \phi_{j} \rangle
  &= \iiint
     \phi_{i}^*(r)
     \times \exp(-br^2) \times
     \phi_{j}(r)
     ~r^2 \sin\theta ~\mathrm{d}r \mathrm{d}\theta \mathrm{d}\varphi \\
  &= \int_0^{2\pi} \mathrm{d}\varphi
     \int_0^\pi \sin\theta ~\mathrm{d}\theta
     \int_0^\infty r^2 \mathrm{e}^{-(b+\alpha_i + \alpha_j) r^2} ~\mathrm{d}r \\
  &= 2\pi \times 2 \times \frac{1!!}{2^{2}} \sqrt{\frac{\pi}{(b + \alpha_i + \alpha_j)^{2\cdot1+1}}} \\
  &= \underline{\left( \frac{\pi}{b + \alpha_i + \alpha_j} \right)^{3/2}}
\end{aligned}
```

Integral Formula:
```math
\int_0^{\infty} r^{2n} \exp \left(-a r^2\right) ~\mathrm{d}r = \frac{(2n-1)!!}{2^{n+1}} \sqrt{\frac{\pi}{a^{2n+1}}}
```
"""
function element(o::GaussianPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * (π/(o.exponent+SGB1.a+SGB2.a))^(3/2)
end

@doc raw"""
`element(o::Hamiltonian, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  H_{ij}
  &= \langle \phi_{i} | \hat{H} | \phi_{j} \rangle \\
  &= \langle \phi_{i} | \sum_k \hat{o}_k | \phi_{j} \rangle \\
  &= \sum_k \langle \phi_{i} | \hat{o}_k | \phi_{j} \rangle \\
\end{aligned}
```
"""
function element(o::Hamiltonian, B1::Basis, B2::Basis)
  return sum(element(term, B1, B2) for term in o.terms)
end
