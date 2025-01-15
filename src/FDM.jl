export FiniteDifferenceMethod #, solve

import FiniteDifferenceMatrices
import SparseArrays
import ArnoldiMethod
import LinearAlgebra
import Subscripts
import Printf

struct FiniteDifferenceMethod
  Δr::Real          # Radial grid spacing
  rₘₐₓ::Real        # This value is not directly used in the calculation, but it is used to determine the `R`.
  R::StepRangeLen   # Radial grid
  l::Int            # Angular momentum quantum number
  direction::Symbol # :c for central, :f for forward, :b for backward
  solver::Symbol    # :ArnoldiMethod or :LinearAlgebra
  function FiniteDifferenceMethod(;Δr=0.1, rₘₐₓ=50.0, R=Δr:Δr:rₘₐₓ, l=0, direction=:c, solver=:LinearAlgebra)
    new(Δr, rₘₐₓ, R, l, direction, solver)
  end
end

Base.string(method::FiniteDifferenceMethod) = "FiniteDifferenceMethod(" * join(["$(symbol)=$(getproperty(method,symbol))" for symbol in fieldnames(typeof(method))], ", ") * ")"
Base.show(io::IO, method::FiniteDifferenceMethod) = print(io, Base.string(method))

function matrix(o::Hamiltonian, method::FiniteDifferenceMethod)
  return sum(matrix(term, method) for term in o.terms)
end

function matrix(o::RestEnergy, method::FiniteDifferenceMethod)
  return SparseArrays.spdiagm([o.m * o.c^2 for r in method.R])
end

function matrix(o::NonRelativisticKinetic, method::FiniteDifferenceMethod)
  D  = FiniteDifferenceMatrices.fdmatrix(Int64(length(method.R)), n=1, m=2, d=method.direction, h=method.Δr, t=typeof(method.Δr))
  D² = FiniteDifferenceMatrices.fdmatrix(Int64(length(method.R)), n=2, m=2, d=method.direction, h=method.Δr, t=typeof(method.Δr))
  return -o.ℏ^2/2/o.m * (D² + SparseArrays.spdiagm(2 ./ method.R) * D - method.l*(method.l+1) * SparseArrays.spdiagm(1 ./ method.R .^ 2))
end

function matrix(o::PotentialTerm, method::FiniteDifferenceMethod)
  return SparseArrays.spdiagm([V(o, r) for r in method.R])
end

function solve(hamiltonian::Hamiltonian, method::FiniteDifferenceMethod; perturbation=Hamiltonian(), info=4, nₘₐₓ=4)
  
  # Initialization
  nₘₐₓ = min(nₘₐₓ, length(method.R))

  # Hamiltonian
  H = matrix(hamiltonian, method)

  # Jacobian
  J = SparseArrays.spdiagm(method.R .^ 2)

  # Eigenvalues
  if method.solver == :LinearAlgebra
    E, C = LinearAlgebra.eigen(Matrix{Float64}(H))
  elseif method.solver == :ArnoldiMethod
    decomp, history = ArnoldiMethod.partialschur(H, nev=nₘₐₓ, tol=1e-9, which=:SR)
    E, C = ArnoldiMethod.partialeigen(decomp)
  else
    error("Unknown solver: $(method.solver)")
  end

  # Print
  if 0 < info
    println("\n# method\n")
    println(method)
    println("\n# eigenvalue\n")
    for n in 1:min(nₘₐₓ, info)
      Printf.@printf("E%s = %s\n", Subscripts.sub("$n"), "$(E[n])")
    end
    println("\n# others")
    println("\nn \tnorm, <ψₙ|ψₙ> = cₙ' * cₙ")
    for n in 1:min(nₘₐₓ, info)
      println("$n\t", C[:,n]' * C[:,n])
    end
    if !isempty(perturbation.terms)
      println("\nn \tperturbation")
      M = matrix(perturbation, method)
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
      M = matrix(term, method)
      for n in 1:min(nₘₐₓ, info)
        println("$n\t", C[:,n]' * M * C[:,n])
      end
    end
    println()
  end

  # Normalization
  N = LinearAlgebra.diagm([1 / sqrt(4 * π * method.Δr * C[:,n]' * J * C[:,n]) for n in 1:nₘₐₓ])
  ψ = C[:, 1:nₘₐₓ] * N

  # Return
  if 0 ≤ info
    return (
      hamiltonian = hamiltonian, 
      perturbation = perturbation,
      method = method,
      nₘₐₓ = nₘₐₓ,
      H = H,
      J = J,
      E = E[1:nₘₐₓ],
      C = C[:, 1:nₘₐₓ],
      ψ = ψ,
    )
  else
    return (
      E = E[1:nₘₐₓ],
    )
  end
end

function solve(hamiltonian::Hamiltonian, wavefunction::Function, method::FiniteDifferenceMethod, info=4, nₘₐₓ=4)

  # Initialization
  nₘₐₓ = min(nₘₐₓ, length(method.R))

  # Hamiltonian
  H = matrix(hamiltonian, method)

  # Jacobian
  J = SparseArrays.spdiagm(method.R .^ 2)

  # Wave Function
  ψ = wavefunction.(method.R)

  # Energy
  E = (ψ' * J * H * ψ) / (ψ' * J * ψ)

  # Return
  if 0 ≤ info
    return (
      hamiltonian = hamiltonian,
      method = method,
      nₘₐₓ = nₘₐₓ,
      H = H,
      J = J,
      E = E,
      ψ = ψ / sqrt(4 * π * method.Δr * ψ' * J * ψ),
    )
  else
    return (
      E = E,
    )
  end

end