@testset "Rayleigh–Ritz.jl" begin

  # Testing Results

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
    ψ(r) = res.ψ[i](r)
    @test isapprox(1, 4π*quadgk(r -> r^2 * abs(ψ(r))^2, 0, Inf, rtol=1e-12, maxevals=10^7)[1], atol=1e-6)
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
  HA = Antique.HydrogenAtom(Z=1, mₑ=1.0, a₀=1.0, Eₕ=1.0, ℏ=1.0)
  for r in 0.2:0.1:2.0
    @test isapprox(abs(Antique.ψ(HA, r, 0, 0)), abs(res.ψ[1](r)), rtol=1e-2)
  end

  # Testing Matrix Elements

  @testset "element()" begin

    # kinetic terms
    o = NonRelativisticKinetic(ℏ=1, m=1)
    println(o)
    println("  i  j\tnumerical         \tanalytical        \t|error|")
    l = 0
    for j in keys(BS.basis)
    for i in keys(BS.basis)
      φᵢ(r) = TwoBody.φ(BS.basis[i], r)
      φⱼ(r) = TwoBody.φ(BS.basis[j], r)
      dφⱼ(r) = ForwardDiff.derivative(r -> φⱼ(r), r)
      d²φⱼ(r) = ForwardDiff.derivative(r -> dφⱼ(r), r)
      numerical  = 4*π*QuadGK.quadgk(r -> r^2 * φᵢ(r) * (-o.ℏ^2/2/o.m * (d²φⱼ(r) + 2/r*dφⱼ(r) - l*(l+1)/r^2*φⱼ(r))), 0, Inf, maxevals=10^3)[1]
      analytical = TwoBody.element(o, BS.basis[i], BS.basis[j])
      error = isinf(analytical) ? 0.0 : abs((numerical-analytical)/analytical)
      acceptance = error < 1e-5
      @printf("%3d%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, j, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
      @test acceptance
    end
    end

    # potential terms
    for o in [
      ConstantPotential(),
      LinearPotential(),
      CoulombPotential(),
      PowerLawPotential(),
      PowerLawPotential(exponent=2),
      GaussianPotential(),
      # ExponentialPotential(),
      # YukawaPotential(),
      # LogarithmicPotential(),
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
