using TwoBody
using Test
using QuadGK

@testset "TwoBody.jl" begin
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
  # Thijssen(2007)
  @test isapprox(-0.499278, res.E[1], atol=1e-6)
  # ψ
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
end
