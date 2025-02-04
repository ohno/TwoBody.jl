export BasisSet, Basis, GeometricBasisSet, PrimitiveBasis, ContractedBasis, SimpleGaussianBasis, GaussianBasis

# struct

@doc raw"""
`Basis` is an abstract type.
"""
abstract type Basis end

@doc raw"""
`PrimitiveBasis <: Basis` is an abstract type.
"""
abstract type PrimitiveBasis <: Basis end

@doc raw"""
`BasisSet(basis1, basis2, ...)`
```math
\{ \phi_1, \phi_2, \phi_3, \cdots  \}
```
The basis set is the input for Rayleigh–Ritz method. You can define the basis set like this:
```math
\begin{aligned}
  \phi_1(r) &= \exp(-13.00773 ~r^2), \\
  \phi_2(r) &= \exp(-1.962079 ~r^2), \\
  \phi_3(r) &= \exp(-0.444529 ~r^2), \\
  \phi_4(r) &= \exp(-0.1219492 ~r^2).
\end{aligned}
```
```@example
BS = BasisSet(
  SimpleGaussianBasis(13.00773),
  SimpleGaussianBasis(1.962079),
  SimpleGaussianBasis(0.444529),
  SimpleGaussianBasis(0.1219492),
)
```
"""
struct BasisSet
  basis::Vector{Basis}
  BasisSet(args...) = new([args...])
end
Base.getindex(BS::BasisSet, index) = BS.basis[index]
Base.length(BS::BasisSet) = length(BS.basis)

@doc raw"""
Exponents of Gaussian basis functions are given by geometric progression:
```math
\begin{aligned}
  & v_i = \frac{1}{r_i^2}, \\
  & r_i = r_1 a^{i-1}.
\end{aligned}
```

This function return array of $\nu_i$:
```math
(r_1, r_{n}, n, n_\mathrm{max}) \mapsto (\nu_1, \nu_2, \cdots, \nu_{n-1}, \nu_n, \nu_{n+1}, \cdots, \nu_{n_\mathrm{max}})
```

Usually $n = n_\mathrm{max}$. Set $n<n_\mathrm{max}$ if you want to extend the geometric progression.

Examples:
```jldoctest
julia> ν = TwoBody.geometric(0.1, 10.0, 5)
5-element Vector{Float64}:
 100.0
  10.0
   0.9999999999999997
   0.09999999999999996
   0.009999999999999995

julia> ν = TwoBody.geometric(0.1, 10.0, 5, nₘₐₓ = 10)
10-element Vector{Float64}:
 100.0
  10.0
   0.9999999999999997
   0.09999999999999996
   0.009999999999999995
   0.0009999999999999994
   9.999999999999994e-5
   9.999999999999992e-6
   9.999999999999991e-7
   9.999999999999988e-8
```
"""
function geometric(r₁, rₙ, n::Int; nₘₐₓ::Int=n, nₘᵢₙ::Int=1)
  rgp = try; abs(rₙ/r₁)^(1/(n-nₘᵢₙ)); catch; 0; end # ratio of geometric sequence
  ν = abs(r₁)^-2 .* rgp.^(-2.0*(nₘᵢₙ-1:nₘₐₓ-1))
  return ν
end

@doc raw"""
`GeometricBasisSet(basistype, r₁, rₙ, n; nₘᵢₙ=1, nₘₐₓ=n)`

This is a basis set with exponentials generated by `geometric()`.
"""
Base.@kwdef struct GeometricBasisSet
  basistype
  r₁::Real
  rₙ::Real
  n::Int
  nₘᵢₙ::Int
  nₘₐₓ::Int
  basis::Vector{Basis}
  function GeometricBasisSet(basistype, r₁, rₙ, n::Int; nₘᵢₙ::Int=1, nₘₐₓ::Int=n)
    new(basistype, r₁, rₙ, n, nₘᵢₙ, nₘₐₓ, [basistype(α) for α in geometric(r₁, rₙ, n, nₘᵢₙ=nₘᵢₙ, nₘₐₓ=nₘₐₓ)])
  end
  function GeometricBasisSet(basistype; r₁=0.1, rₙ=80.0, n::Int=20, nₘᵢₙ::Int=1, nₘₐₓ::Int=n)
    new(basistype, r₁, rₙ, n, nₘᵢₙ, nₘₐₓ, [basistype(α) for α in geometric(r₁, rₙ, n, nₘᵢₙ=nₘᵢₙ, nₘₐₓ=nₘₐₓ)])
  end
end
Base.getindex(GBS::GeometricBasisSet, index) = GBS.basis[index]
Base.length(GBS::GeometricBasisSet) = length(GBS.basis)

@doc raw"""
`SimpleGaussianBasis(a=1)`
```math
\phi_i(r) = \exp(-a_i r^2)
```
Note: This basis is not normalized. This is only for s-wave.
"""
Base.@kwdef struct SimpleGaussianBasis <: PrimitiveBasis
  a = 1
end

@doc raw"""
`GaussianBasis(a=1, l=0, m=0)`
```math
\phi_i(r, θ, φ) = N _{il} r^l \exp(-a_i r^2) Y_l^m(θ, φ)
```
"""
Base.@kwdef struct GaussianBasis <: PrimitiveBasis
  a = 1
  l = 0
  m = 0
end

@doc raw"""
`ContractedBasis([c1, c2, ...], [basis1, basis2, ...])`
```math
\phi' = \sum_i c_i \phi_i
```
"""
struct ContractedBasis <: Basis
  c::Vector
  φ::Vector
end

# display

Base.string(b::Basis) = "$(typeof(b))(" * join(["$(symbol)=$(getproperty(b,symbol))" for symbol in fieldnames(typeof(b))], ", ") * ")"
Base.show(io::IO, b::Basis) = print(io, Base.string(b))
Base.string(bs::BasisSet) = "BasisSet(" * join(["$(b)" for b in bs.basis], ", ") * ")"
# Base.string(bs::BasisSet) = "BasisSet(\n" * join(["  $(b)" for b in bs.basis], ",\n") * ",\n)"
Base.show(io::IO, bs::BasisSet) = print(io, Base.string(bs))

# function

φ(b::SimpleGaussianBasis, r) = exp(-b.a*r^2)
# φ(b::GaussianBasis, r, θ, φ) = N(b.l) * r^b.l * exp(-b.a*r^2) * Y(b.l, b.m, θ, φ)
φ(b::ContractedBasis, r, θ, φ) = sum(b.c[i] * φ(b.φ[i], r, θ, φ) for i in 1:length(b.c))
