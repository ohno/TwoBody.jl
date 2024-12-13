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
  @test isapprox(res.E[1], -0.499278, rtol=1e-6)
  # normalization
  @test quadgk(r -> 4π*r^2*res.ψ[1](r)^2, 0, Inf, rtol=1e-12, maxevals=10^7)[1] ≈ 1
end
