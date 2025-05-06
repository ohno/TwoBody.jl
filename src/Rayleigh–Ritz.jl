export solve, optimize

import LinearAlgebra
import Optim
import Printf
import SpecialFunctions
import Subscripts

# solver

function solve(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4)

  # Initialization
  nₘₐₓ = length(basisset.basis)

  # Matrix Elements
  S = matrix(basisset)
  H = matrix(hamiltonian, basisset)

  # Calculations
  E, C = LinearAlgebra.eigen(H, S)

  # Print
  if 0 < info
    println("\n# method\n")
    println("Rayleigh–Ritz method with $(typeof(basisset.basis[1]))")
    println("J. Thijssen, Computational Physics 2nd Edition (2013)")
    println("https://doi.org/10.1017/CBO9781139171397")
    println("\n# basis function\n")
    for n in 1:nₘₐₓ
      Printf.@printf("φ%s(r) = TwoBody.φ(%s, r)\n", Subscripts.sub("$n"), basisset.basis[n])
    end
    println("\n# eigenfunction\n")
    for n in 1:min(nₘₐₓ, info)
      Printf.@printf("ψ%s(r) = ", Subscripts.sub("$n"))
      for i in 1:nₘₐₓ
        Printf.@printf("%s %.6fφ%s(r) ", C[i,n]<0 ? "-" : "+", abs(C[i,n]), Subscripts.sub("$i"))
      end
      println()
    end
    println("\n# eigenvalue\n")
    for n in 1:min(nₘₐₓ, info)
      Printf.@printf("E%s = %s\n", Subscripts.sub("$n"), "$(E[n])")
    end
    println("\n# others")
    println("\nn \tnorm, <ψₙ|ψₙ> = cₙ' * S * cₙ")
    for n in 1:min(nₘₐₓ, info)
      println("$n\t", C[:,n]' * S * C[:,n])
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
    println("\nn \terror check, |<ψₙ|H|ψₙ> - E| = |cₙ' * H * cₙ - E| = 0")
    for n in 1:min(nₘₐₓ, info)
      println("$n\t", abs(C[:,n]' * H * C[:,n] - E[n]))
    end
    for term in [hamiltonian.terms..., perturbation.terms...]
      println("\nn \texpectation value of $(term)")
      M = [element(term, basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
      for n in 1:min(nₘₐₓ, info)
        println("$n\t", C[:,n]' * M * C[:,n])
      end
    end
    println()
  end

  # Return
  if 0 ≤ info
    return (
      hamiltonian = hamiltonian, 
      perturbation = perturbation,
      basisset = basisset,
      nₘₐₓ = nₘₐₓ,
      H = H,
      S = S,
      E = E,
      C = C,
      φ = [r -> TwoBody.φ(basisset.basis[n], r) for n in 1:nₘₐₓ],
      ψ = [r -> sum(C[i,n]*TwoBody.φ(basisset.basis[i], r) for i in 1:nₘₐₓ) for n in 1:nₘₐₓ],
    )
  else
    return (
      E = E,
    )
  end
end

function optimize(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)

  # optimizer & initial values
  if 0 < info
    println("\n# optimizer\n")
    println(optimizer)
    println("P. K. Mogensen, A. N. Riseth, J. Open Source Softw., 3(24), 615 (2018)")
    println("https://doi.org/10.21105/joss.00615")
    println("\n# initial basis function\n")
    nₘₐₓ = length(basisset.basis)
    for n in 1:nₘₐₓ
      Printf.@printf("φ%s(r) = TwoBody.φ(%s, r)\n", Subscripts.sub("$n"), basisset.basis[n])
    end
  end
  if 0 < info && progress
    println("\n# optimization log\n")
  end

  # optimize
  history = []
  res = Optim.optimize(
    x -> try
      E = solve(
        hamiltonian,
        BasisSet([typeof(basisset.basis[i])(x[i]) for i in keys(basisset.basis)]...),
        perturbation = perturbation,
        info = -1
      ).E[1]
      if 0 ≤ info
        push!(history, (energy=E, parameters=x))
      end
      if 0 < info && progress
        Printf.@printf("%.9e  %s\n", E, "[" * join([Printf.@sprintf("%.3e", x[i]) for i in keys(x)], ", ") *"]")
      end
    E
    catch
      if 0 ≤ info
        push!(history, (energy=Inf, parameters=x))
      end
      if 0 < info && progress
        Printf.@printf("%.9e  %s\n", Inf, "[" * join([Printf.@sprintf("%.3e", x[i]) for i in keys(x)], ", ") *"]")
      end
      Inf
    end,
    [basisset.basis[i].a for i in keys(basisset.basis)],
    method = optimizer,
    options...
  )

  # final results
  res = solve(hamiltonian, BasisSet([typeof(basisset.basis[i])(res.minimizer[i]) for i in keys(basisset.basis)]...), perturbation=perturbation, info=info)
  if 0 ≤ info
    return (res..., history=history)
  else
    return res
  end

end

function solve(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=4)
  return solve(hamiltonian, BasisSet(basis); perturbation=perturbation, info=info)
end

function optimize(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=1, progress=true, optimizer=Optim.NelderMead(), options...)
  return optimize(hamiltonian, BasisSet(basis); perturbation=perturbation, info=info, progress=progress, optimizer=optimizer, options...)
end

function solve(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4)
  if 0 < info
    println("\n# geometric progression\n")
    println("type \t$(basisset.basistype)")
    println("range\tr", Subscripts.sub("$(basisset.nₘᵢₙ)"), " - r", Subscripts.sub("$(basisset.nₘₐₓ)"))
    println("r", Subscripts.sub("$(basisset.nₘᵢₙ)"), " \t", basisset.r₁)
    println("r", Subscripts.sub("$(basisset.n)"  ), " \t", basisset.rₙ)
  end
  solve(hamiltonian, BasisSet(basisset.basis...); perturbation=perturbation, info=info)
end

function optimize(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)

  # optimizer & initial values
  if 0 < info && progress
    println("\n# optimizer\n")
    println(optimizer)
    println("P. K. Mogensen, A. N. Riseth, J. Open Source Softw., 3(24), 615 (2018)")
    println("https://doi.org/10.21105/joss.00615")
    println("\n# initial geometric progression\n")
    println("type \t$(basisset.basistype)")
    println("range\tr", Subscripts.sub("$(basisset.nₘᵢₙ)"), " - r", Subscripts.sub("$(basisset.nₘₐₓ)"))
    println("r", Subscripts.sub("$(basisset.nₘᵢₙ)"), " \t", basisset.r₁)
    println("r", Subscripts.sub("$(basisset.n)"  ), " \t", basisset.rₙ)
    println("\n# optimization log\n")
  end

  # optimize
  history = []
  res = Optim.optimize(
    x -> try
      E = solve(
        hamiltonian,
        GeometricBasisSet(basisset.basistype, x..., basisset.n, nₘₐₓ=basisset.nₘₐₓ, nₘᵢₙ=basisset.nₘᵢₙ),
        perturbation = perturbation,
        info = -1
      ).E[1]
      if 0 ≤ info
        push!(history, (energy=E, parameters=x))
      end
      if 0 < info && progress
        Printf.@printf("%.9e  %s\n", E, "[" * join([Printf.@sprintf("%+.3e", x[i]) for i in keys(x)], ", ") *"]")
      end
    E
    catch
      if 0 ≤ info
        push!(history, (energy=Inf, parameters=x))
      end
      if 0 < info && progress
        Printf.@printf("%.9e  %s\n", Inf, "[" * join([Printf.@sprintf("%+.3e", x[i]) for i in keys(x)], ", ") *"]")
      end
      Inf
    end,
    [basisset.r₁, basisset.rₙ],
    method = optimizer,
    options...
  )

  # final results
  res = solve(hamiltonian, GeometricBasisSet(basisset.basistype, res.minimizer..., basisset.n, nₘₐₓ=basisset.nₘₐₓ, nₘᵢₙ=basisset.nₘᵢₙ); perturbation=perturbation, info=info)
  if 0 ≤ info
    return (res..., history=history)
  else
    return res
  end

end

# element

function element(SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return (π/(SGB1.a+SGB2.a))^(3/2)
end

function element(o::RestEnergy, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.m * o.c^2 * (π/(SGB1.a+SGB2.a))^(3/2)
end

function element(o::Laplacian, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return - 6*π^(3/2)*SGB1.a*SGB2.a/(SGB1.a+SGB2.a)^(5/2)
end

function element(o::NonRelativisticKinetic, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.ℏ^2/(2*o.m) * 6*π^(3/2)*SGB1.a*SGB2.a/(SGB1.a+SGB2.a)^(5/2)
end

function element(o::ConstantPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.constant * (π/(SGB1.a+SGB2.a))^(3/2)
end

function element(o::LinearPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * 2*π/(SGB1.a+SGB2.a)^2
end

function element(o::CoulombPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * 2*π/(SGB1.a+SGB2.a)
end

function element(o::PowerLawPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * 2 * π * SpecialFunctions.gamma((o.exponent+3)/2) / (SGB1.a+SGB2.a)^((o.exponent+3)/2)
end

function element(o::GaussianPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)
  return o.coefficient * (π/(o.exponent+SGB1.a+SGB2.a))^(3/2)
end

function element(o::Hamiltonian, B1::Basis, B2::Basis)
  return sum(element(term, B1, B2) for term in o.terms)
end

# matrix

function matrix(basisset::BasisSet)
  nₘₐₓ = length(basisset.basis)
  # S = [element(basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
  S = Array{Float64}(undef, nₘₐₓ, nₘₐₓ)
  for j in 1:nₘₐₓ
    for i in 1:j
      S[i,j] = element(basisset.basis[i], basisset.basis[j])
    end
  end
  return LinearAlgebra.Symmetric(S)
end

function matrix(hamiltonian::Hamiltonian, basisset::BasisSet)
  nₘₐₓ = length(basisset.basis)
  # H = [element(hamiltonian, basisset.basis[i], basisset.basis[j]) for i=1:nₘₐₓ, j=1:nₘₐₓ]
  H = Array{Float64}(undef, nₘₐₓ, nₘₐₓ)
  for j in 1:nₘₐₓ
    for i in 1:j
      H[i,j] = element(hamiltonian, basisset.basis[i], basisset.basis[j])
    end
  end
  return LinearAlgebra.Symmetric(H)
end

# docstring

@doc raw"""
`solve(hamiltonian::Hamiltonian, basisset::BasisSet)`

This function returns the eigenvalues ``E``  and eigenvectors ``\pmb{c}`` for
```math
\pmb{H} \pmb{c} = E \pmb{S} \pmb{c}.
```
The Hamiltonian matrix is defined as ``H_{ij} = \langle \phi_{i} | \hat{H} | \phi_{j} \rangle``. The overlap matrix is defined as ``S_{ij} = \langle \phi_{i} | \phi_{j} \rangle``.
""" solve(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4)

@doc raw"""
`function optimize(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)`

This function minimizes the energy by changing the exponents of the basis functions using Optim.jl.
```math
\frac{\partial E}{\partial a_i} = 0
```
""" optimize(hamiltonian::Hamiltonian, basisset::BasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)

@doc raw"""
`solve(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=4)`

This a solver for 1-basis calculations. This function returns `solve(hamiltonian, BasisSet(basis); perturbation=perturbation, info=info)`.
""" solve(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=4)

@doc raw"""
`optimize(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=4, optimizer=Optim.NelderMead())`

This a optimizer for 1-basis calculations. This function returns `optimize(hamiltonian, BasisSet(basis); perturbation=perturbation, info=info, progress=progress, optimizer=optimizer, options...)`.
""" optimize(hamiltonian::Hamiltonian, basis::Basis; perturbation=Hamiltonian(), info=1, progress=true, optimizer=Optim.NelderMead(), options...)

@doc raw"""
`solve(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4)`

This function is a wrapper for `solve(hamiltonian::Hamiltonian, basisset::BasisSet, ...)`.
""" solve(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4)

@doc raw"""
`optimize(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4, optimizer=Optim.NelderMead())`

This function minimizes the energy by optimizing $r_1$ and $r_n$ using Optim.jl.
```math
\frac{\partial E}{\partial r_1} = \frac{\partial E}{\partial r_n} = 0
```
""" optimize(hamiltonian::Hamiltonian, basisset::GeometricBasisSet; perturbation=Hamiltonian(), info=4, progress=true, optimizer=Optim.NelderMead(), options...)

@doc raw"""
`element(SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  S_{ij}
   = \langle \phi_{i} | \phi_{j} \rangle
  &= \int
     \phi_{i}^*(r)
     \phi_{j}(r)
     \mathrm{d} \pmb{r} \\
  &= \iiint
     \mathrm{e}^{-\alpha_i r^2}
     \mathrm{e}^{-\alpha_j r^2}
     ~r^2 \sin\theta ~
     \mathrm{d} r
     \mathrm{d} \theta
     \mathrm{d} \varphi \\
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
""" element(SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

@doc raw"""
`element(o::RestEnergy, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | mc^2 | \phi_{j} \rangle
  &= mc^2 \langle \phi_{i} | \phi_{j} \rangle \\
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
""" element(o::RestEnergy, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

@doc raw"""
`element(o::Laplacian, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | \nabla^2 | \phi_{j} \rangle
  = \underline{
      -6 \frac{\alpha_i \alpha_j \pi^{\frac{3}{2}}}{(\alpha_i + \alpha_j)^{\frac{5}{2}}}
    }
\end{aligned}
```
""" element(o::Laplacian, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

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
        \frac{\hbar^2}{2\mu}
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
""" element(o::NonRelativisticKinetic, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

@doc raw"""
`element(o::ConstantPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)`

```math
\begin{aligned}
  \langle \phi_{i} | c | \phi_{j} \rangle
  &= c \langle \phi_{i} | \phi_{j} \rangle \\
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
""" element(o::ConstantPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

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
""" element(o::LinearPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

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
""" element(o::CoulombPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

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
""" element(o::PowerLawPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

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
""" element(o::GaussianPotential, SGB1::SimpleGaussianBasis, SGB2::SimpleGaussianBasis)

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
""" element(o::Hamiltonian, B1::Basis, B2::Basis)

@doc raw"""
`matrix(basisset::BasisSet)`

This function returns the overlap matrix $\pmb{S}$. The element is written as ``S_{ij} = \langle \phi_{i} | \phi_{j} \rangle``.
""" matrix(basisset::BasisSet)

@doc raw"""
`matrix(hamiltonian::Hamiltonian, basisset::BasisSet)`

This function returns the Hamiltonian matrix $\pmb{H}$. The element is written as ``H_{ij} = \langle \phi_{i} | \hat{H} | \phi_{j} \rangle``.
""" matrix(hamiltonian::Hamiltonian, basisset::BasisSet)
