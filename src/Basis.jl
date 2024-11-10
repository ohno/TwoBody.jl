export BasisSet, Basis, PrimitiveBasis, ContractedBasis, SimpleGaussianBasis, GaussianBasis

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
basisset = BasisSet(
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
Base.show(io::IO, bs::BasisSet) = print(io, Base.string(bs))

# function

φ(b::SimpleGaussianBasis, r) = exp(-b.a*r^2)
# φ(b::GaussianBasis, r, θ, φ) = N(b.l) * r^b.l * exp(-b.a*r^2) * Y(b.l, b.m, θ, φ)
φ(b::ContractedBasis, r, θ, φ) = sum(b.c[i] * φ(b.φ[i], r, θ, φ) for i in 1:length(b.c))
