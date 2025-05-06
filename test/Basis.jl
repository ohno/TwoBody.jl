@testset "Basis.jl" begin

  BS = BasisSet(
    SimpleGaussianBasis(13.00773),
    SimpleGaussianBasis(1.962079),
    SimpleGaussianBasis(0.444529),
    SimpleGaussianBasis(0.1219492),
  )

  println("φ(SGB,r) = exp(-a*r^2)")
  println("    a\t    r\tnumerical  \tanalytical")
  @testset "SimpleGaussianBasis" begin
    for b in BS.basis
      a = b.a
      for r in [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        numerical = TwoBody.φ(b,r)
        analytical = exp(-a*r^2)
        acceptance = abs(analytical)<1e-5 ? isapprox(analytical, numerical, atol=1e-2) : isapprox(analytical, numerical, rtol=1e-2)
        @printf("%5.2f\t%5.2f\t%.9f\t%.9f\t%s\n", a, r, numerical, analytical, acceptance ? "✔" :  "✗")
        @test acceptance
      end
    end
  end

  println("φ(r) = exp(-ar²)")
  println("φ(k) = exp(-k²/4a) / (2a)^(3/2)")
  println("φ(k) = 1/√2π³ ∫ φ(r) eⁱᵏʳ r²sin(θ)drdθdφ")
  println("    a\t    k\t θk\t φk\tnumerical  \tanalytical")
  for b in BS.basis
    a = b.a
    for k in [0.0, 3.0, 5.0]
    for θk in [0.0, 0.5]
    for φk in [0.0, 1.0]
      numerical = abs(
        quadgk(φ ->
        quadgk(θ ->
        quadgk(r ->
          (2π)^(-3/2) * TwoBody.φ(b,r) * TwoBody.expikr(k,θk,φk,r,θ,φ) * r^2 * sin(θ)
        , 0, Inf, maxevals=100)[1]
        , 0,   π, maxevals=20)[1]
        , 0,  2π, maxevals=40)[1]
      )
      analytical = exp(-k^2/4/a) / (2*a)^(3/2)
      acceptance = abs(analytical)<1e-5 ? isapprox(analytical, numerical, atol=1e-2) : isapprox(analytical, numerical, rtol=1e-2)
      @test acceptance
      @printf("%5.2f\t%5.2f\t%.1f\t%.1f\t%.9f\t%.9f\t%s\n", a, k, θk, φk, numerical, analytical, acceptance ? "✔" :  "✗")
    end
    end
    end
  end

  # @testset "specialbesselj(n,x)" begin
  #   for n in [0, 1, 2, 3]
  #     for x in [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
  #       numerical = SpecialFunctions.sphericalbesselj(n, x)
  #       analytical = [
  #         sin(x)/x
  #         sin(x)/x^2 - cos(x)/x
  #         (3/x^2-1)*sin(x)/x - 3*cos(x)/x^2
  #         (15/x^3-6/x)*sin(x)/x - (15/x^2-1)*cos(x)/x
  #       ][n+1]
  #       error = (analytical==0.0 || isnan(analytical) || isnan(numerical)) ? 0.0 : abs((numerical-analytical)/analytical)
  #       acceptance = error <1e-5
  #       @printf("%1d  %.1f\t%.9f\t%.9f\t%s\n", n, x, numerical, analytical, acceptance ? "✔" :  "✗")
  #       @test acceptance
  #     end
  #   end
  # end

  # @testset "expirk(r,θr,φr,k,θk,φk)" begin
  #   l_max = 10
  #   for  r in [0.0, 0.5, 1.0, 1.5]
  #   for θr in [0.0, 0.5]
  #   for φr in [0.0, 0.5]
  #   for  k in [0.0, 0.5, 1.0, 1.5]
  #   for θk in [0.0, 0.5]
  #   for φk in [0.0, 0.5]
  #     numerical = real(4*π*sum(sum(im^l * sphericalbesselj(l,r*k) * Y(θp,φp,l=l,m=m) * conj(Y(θr,φr,l=l,m=m)) for m in -l:l) for l in 0:l_max))
  #     analytical = real(expirp(r,θr,φr,k,θk,φk))
  #     error = abs((numerical-analytical)/analytical)
  #     acceptance = error <1e-5
  #     @test acceptance
  #     @printf("%.1f  %.1f  %.1f  %.1f  %.1f  %.1f\t%.9f\t%.9f\t%s\n", r, θr, φr, k, θk, φk, numerical, analytical, acceptance ? "✔" :  "✗")
  #   end
  #   end
  #   end
  #   end
  #   end
  #   end
  # end

  # @testset "φk(k,thetak,phip;ν=1.0,l=0,m=0)" begin
  #   println("  l\t ν\t  p\tθp\tφp\tnumerical        \tanalytical       ")
  #   for l in [0,1,2,3]
  #   for ν in [1.0, 7.0]
  #   for p in [0.0, 3.0, 5.0]
  #   for thetap in [0.0, 0.5]
  #   for phip in [0.0, 1.0]
  #     numerical = real(
  #       quadgk(phi ->
  #       quadgk(theta ->
  #       quadgk(r ->
  #         (2π)^(-3/2) * φ(r,theta,phi;ν=ν,l=l,m=0) * expirp(r,theta,phi,p,thetap,phip) * r^2 * sin(theta)
  #       , 0, Inf, maxevals=100)[1]
  #       , 0, π, maxevals=20)[1]
  #       , 0, 2π, maxevals=40)[1]
  #     )
  #     analytical = real(φp(p,thetap,phip,ν=ν,l=l,m=0))
  #     error = analytical==0.0 ? 0.0 : abs((numerical-analytical)/analytical)
  #     acceptance = error <1e-5
  #     @test acceptance
  #     @printf("%3d\t%.1f\t%.1f\t%.1f\t%.1f\t%.9f\t%.9f\t%s\n", l, ν, p, thetap, phip, numerical, analytical, acceptance ? "✔" :  "✗")
  #   end
  #   end
  #   end
  #   end
  #   end
  # end

end
