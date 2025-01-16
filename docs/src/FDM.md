```@meta
CurrentModule = TwoBody
```

# Finite Difference Method

```@example example
using TwoBody
```

```@example example
H = Hamiltonian(
  NonRelativisticKinetic(ℏ = 1 , m = 1),
  CoulombPotential(coefficient = -1),
)
nothing # hide
```

```@example example
FDM = FiniteDifferenceMethod(
  Δr = 0.1,
  rₘₐₓ = 50.0,
  l = 0,
  direction = :c,
  solver = :LinearAlgebra,
)
nothing # hide
```

```@repl example
solve(H, FDM)
```

```@example example
# solve
using TwoBody
H = Hamiltonian(NonRelativisticKinetic(1,1), CoulombPotential(-1))
FDM = TwoBody.FiniteDifferenceMethod()
res = solve(H, FDM, info=0, nₘₐₓ=4)

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
  X = res.method.R
  Y = 4π * X .^2 .* res.ψ[:,n] .^ 2
  scatter!(axis, X, Y, label="TwoBody.jl", markersize=6)
  lines!(axis, 0..50, r -> 4π * r^2 * abs(Antique.ψ(HA,r,0,0,n=n))^2, label="Antique.jl", color=:black)
  axislegend(axis, "n = $n", position=:rt, framevisible=false)
end
fig
```