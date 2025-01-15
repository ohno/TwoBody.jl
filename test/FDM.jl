@testset "FDM.jl" begin

  # Testing Results

  H = Hamiltonian(
    NonRelativisticKinetic(ℏ = 1 , m = 1),
    CoulombPotential(coefficient = -1),
  )
  FDM = FiniteDifferenceMethod(
    Δr = 0.1,
    rₘₐₓ = 50.0,
    l = 0,
    direction = :c,
    solver = :LinearAlgebra,
  )
  res = solve(H, FDM, info=4)
  
  println("<ψₙ|ψₙ> = cₙ' * cₙ = 1")
  println("  i\tnumerical         \tanalytical        \t|error|")
  for i in 1:res.nₘₐₓ
    numerical  = res.C[:,i]' * res.C[:,i]
    analytical = 1
    error = iszero(analytical) ? abs(numerical-analytical) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @printf("%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
    @test acceptance
  end
  
  println("|<ψₙ|H|ψₙ> - E| = |cₙ' * H * cₙ - E| = 0")
  println("  i\tnumerical         \tanalytical        \t|error|")
  for i in 1:res.nₘₐₓ
    numerical  = abs((res.C[:,i]' * res.H * res.C[:,i]) - res.E[i])
    analytical = 0
    error = iszero(analytical) ? abs(numerical-analytical) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-5
    @printf("%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
    @test acceptance
  end
  
  # comparison with Antique.jl
  HA = Antique.HydrogenAtom(Z=1, mₑ=1.0, a₀=1.0, Eₕ=1.0, ℏ=1.0)

  println("Energy")
  println("  i\tnumerical         \tanalytical        \t|error|")
  for i in 1:res.nₘₐₓ
    numerical  = res.E[i]
    analytical = Antique.E(HA, n=i)
    error = iszero(analytical) ? abs(numerical-analytical) : abs((numerical-analytical)/analytical)
    acceptance = error < 1e-2
    @printf("%3d\t%.16f\t%.16f\t%.16f%%\t%s\n", i, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
    @test acceptance
  end

  println("Wave Function")
  println("  i\t  r\tnumerical         \tanalytical        \t|error|")
  for n in 1:res.nₘₐₓ
    @show n
    for i in keys(res.method.R[begin:min(10,length(res.method.R))])
        r = res.method.R[i]
        numerical  = abs(res.ψ[i,n])
        analytical = abs(Antique.ψ(HA, r, 0, 0, n=n)) 
        error = iszero(analytical) ? abs(numerical-analytical) : abs((numerical-analytical)/analytical)
        acceptance = error < 5e-2
        @printf("%3d\t%.1f\t%.16f\t%.16f\t%.16f%%\t%s\n", i, r, numerical, analytical, error*100, acceptance ? "✔" :  "✗")
        @test acceptance
    end
  end

end
