using TwoBody
using Test
using QuadGK
using Printf
using Antique

@testset "Rayleigh–Ritz.jl" begin
  H = Hamiltonian(
    NonRelativisticKinetic(ℏ = 1 , m = 1),
    CoulombPotential(coefficient = -1),
  )
  BS = BasisSet(
    SimpleGaussianBasis(13.00773),
    SimpleGaussianBasis(1.962079),
    SimpleGaussianBasis(0.444529),
    SimpleGaussianBasis(0.1219492),
  )
  res = solve(H, BS)
  # 4π×∫|ψ(r)|²r²dr = 1
  for i in keys(res.E)
    @test isapprox(1, 4π *quadgk(r -> r^2 * abs(res.ψ[i](r))^2, 0, Inf, rtol=1e-12, maxevals=10^7)[1], atol=1e-6)
  end
  # <ψₙ|ψₙ> = cₙ' * S * cₙ = 1
  for i in keys(res.E)
    @test isapprox(1, res.C[:,i]' * res.S * res.C[:,i], atol=1e-6)
  end
  # |<ψₙ|H|ψₙ> - E| = |cₙ' * H * cₙ - E| = 0
  for i in keys(res.E)
    @test isapprox(0, abs((res.C[:,i]' * res.H * res.C[:,i]) - res.E[i]), atol=1e-6)
  end
  # Thijssen(2007)
  @test isapprox(-0.499278, res.E[1], atol=1e-6)
  # comparison with Antique.jl
  HA = HydrogenAtom(Z=1, mₑ=1.0, a₀=1.0, Eₕ=1.0, ℏ=1.0)
  for r in 0.2:0.1:2.0
    @show r
    @test isapprox(abs(Antique.ψ(HA, r, 0, 0)), abs(res.ψ[1](r)), rtol=1e-2)
  end
  # elements
  @testset "element()" begin
    for o in [
      ConstantPotential(),
      LinearPotential(),
      CoulombPotential(),
      PowerLawPotential(),
      PowerLawPotential(exponent=2),
      GaussianPotential(),
      # ExponentialPotential(),
      # YukawaPotential(),
      # DeltaPotential(),
    ]
      println(o)
      println("  i  j\tnumerical         \tanalytical        \t|error|")
      for j in keys(BS.basis)
      for i in keys(BS.basis)
        numerical  = 4*π*quadgk(r -> r^2 * TwoBody.V(o, r) * TwoBody.φ(BS.basis[i],r) * TwoBody.φ(BS.basis[j],r), 0, Inf, maxevals=10^3)[1]
        analytical = TwoBody.element(o, BS.basis[i], BS.basis[j])
        error = isinf(analytical) ? 0.0 : abs((numerical-analytical)/analytical)
        acceptance = error < 1e-5
        @printf("%3d%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, j, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
        @test acceptance
      end
      end
    end
  end
end
