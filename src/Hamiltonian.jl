export Operator, Hamiltonian, NonRelativisticKinetic, RestEnergy, RelativisticCorrection, RelativisticKinetic, ConstantPotential, LinearPotential, CoulombPotential, PowerLawPotential, GaussianPotential, ExponentialPotential, YukawaPotential, DeltaPotential, FunctionPotential, UniformGridPotential

# struct

@doc raw"""
`Operator` is an abstract type.
"""
abstract type Operator end

@doc raw"""
`PotentialTerm <: Operator` is an abstract type.
"""
abstract type PotentialTerm <: Operator end

@doc raw"""
`KineticTerm <: Operator` is an abstract type.
"""
abstract type KineticTerm <: Operator end

@doc raw"""
`Hamiltonian(operator1, operator2, ...)`
```math
\hat{H} = \sum_i \hat{o}_i
```

The Hamiltonian is the input for each solver. This is an example for the non-relativistic Hamiltonian of hydrogen atom in atomic units:

```math
\hat{H} = 
- \frac{1}{2} \nabla^2
- \frac{1}{r}
```

```@example
hamiltonian = Hamiltonian(
  NonRelativisticKinetic(ℏ =1 , m = 1),
  CoulombPotential(coefficient = -1),
)
```
"""
struct Hamiltonian
  terms::Array{Operator,1}
  Hamiltonian(args...) = new([args...])
end

@doc raw"""
`NonRelativisticKinetic(ℏ=1, m=1)`
```math
-\frac{\hbar^2}{2m} \nabla^2
```
"""
Base.@kwdef struct NonRelativisticKinetic <: KineticTerm
  ℏ = 1
  m = 1
end

@doc raw"""
`RestEnergy(c=1, m=1)`
```math
m c^2
```
Use `c = 137.035999177` (from [2022 CODATA](https://physics.nist.gov/cgi-bin/cuu/Value?alphinv)) in the atomic units.
"""
Base.@kwdef struct RestEnergy <: KineticTerm
  c = 1
  m = 1
end

@doc raw"""
`RelativisticCorrection(c=1, m=1, n=2)`
The p^{2n} term of the Taylor expansion:
```math
\begin{aligned}
  \sqrt{p^2 c^2 + m^2 c^4}
  =& m \times c^2 \\
  &+ 1 / 2   / m         \times p^2 (n=1) \\
  &- 1 / 8   / m^3 / c^2 \times p^4 (n=2) \\
  &+ 1 / 16  / m^5 / c^4 \times p^6 (n=3) \\
  &- 5 / 128 / m^7 / c^6 \times p^8 (n=4) \\
  &+ \cdots
\end{aligned}
```
Use `c = 137.035999177` (from [2022 CODATA](https://physics.nist.gov/cgi-bin/cuu/Value?alphinv)) in the atomic units.
"""
Base.@kwdef struct RelativisticCorrection <: KineticTerm
  c = 1
  m = 1
  n = 2
end

@doc raw"""
`RelativisticKinetic(c=1, m=1)`
```math
\sqrt{p^2 c^2 + m^2 c^4} - m c^2
```
Use `c = 137.035999177` (from [2022 CODATA](https://physics.nist.gov/cgi-bin/cuu/Value?alphinv)) in the atomic units.
"""
Base.@kwdef struct RelativisticKinetic <: KineticTerm
  c = 1
  m = 1
end

@doc raw"""
`ConstantPotential(constant=1)`
```math
+ c
```
| Arguments | Description |
| :-- | :-- |
| `constant` | ``c`` |
"""
Base.@kwdef struct ConstantPotential <: PotentialTerm
  constant = 1
end

@doc raw"""
`LinearPotential(coefficient=1)`
```math
+ ar
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
"""
Base.@kwdef struct LinearPotential <: PotentialTerm
  coefficient = 1
end

@doc raw"""
`CoulombPotential(coefficient=1)`
```math
+ \frac{a}{r}
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
"""
Base.@kwdef struct CoulombPotential <: PotentialTerm
  coefficient = 1
end

@doc raw"""
`PowerLawPotential(coefficient=1, exponent=1)`
```math
+ ar^n
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
| `exponent` | ``n`` |
"""
Base.@kwdef struct PowerLawPotential <: PotentialTerm
  coefficient = 1
  exponent = 1
end

@doc raw"""
`GaussianPotential(coefficient=1, exponent=1)`
```math
+ a \exp(- b r^2)
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
| `exponent`    | ``b`` |
"""
Base.@kwdef struct GaussianPotential <: PotentialTerm
  coefficient = 1
  exponent = 1
end

@doc raw"""
`ExponentialPotential(coefficient=1, exponent=1)`
```math
+ a \exp(- b r)
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
| `exponent`    | ``b`` |
"""
Base.@kwdef struct ExponentialPotential <: PotentialTerm
  coefficient = 1
  exponent = 1
end

@doc raw"""
`YukawaPotential(coefficient=1, exponent=1)`
```math
+ \frac{a}{r} \exp(- b r)
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
| `exponent`    | ``b`` |
"""
Base.@kwdef struct YukawaPotential <: PotentialTerm
  coefficient = 1
  exponent = 1
end

@doc raw"""
`DeltaPotential(coefficient=1)`
```math
+ a δ(r)
```
| Arguments | Description |
| :-- | :-- |
| `coefficient` | ``a`` |
"""
Base.@kwdef struct DeltaPotential <: PotentialTerm
  coefficient = 1
end

@doc raw"""
`FunctionPotential(f)`
```math
+ f(r)
```
"""
Base.@kwdef struct FunctionPotential <: PotentialTerm
  f::Function
end

@doc raw"""
`UniformGridPotential(R, V)`
"""
Base.@kwdef struct UniformGridPotential <: PotentialTerm
  R::StepRangeLen
  V::Array{Number,1}
end

# display

Base.string(t::Operator) = "$(typeof(t))(" * join(["$(symbol)=$(getproperty(t,symbol))" for symbol in fieldnames(typeof(t))], ", ") * ")"
Base.show(io::IO, t::Operator) = print(io, Base.string(t))
Base.string(H::Hamiltonian) = "Hamiltonian(" * join(["$(term)" for term in H.terms], ", ") * ")"
Base.show(io::IO, H::Hamiltonian) = print(io, Base.string(H))

# function

V(p::ConstantPotential   , r) = p.constant
V(p::LinearPotential     , r) = p.coefficient * r
V(p::CoulombPotential    , r) = p.coefficient / r
V(p::PowerLawPotential   , r) = p.coefficient * r ^ p.exponent
V(p::GaussianPotential   , r) = p.coefficient * exp(- p.exponent * r ^ 2)
V(p::ExponentialPotential, r) = p.coefficient * exp(- p.exponent * r)
V(p::YukawaPotential     , r) = p.coefficient * exp(- p.exponent * r) / r
# V(p::DeltaPotential      , r) = 
V(p::FunctionPotential   , r) = p.f(r)
V(p::UniformGridPotential, r) = p.V[findfirst(p.R, r)]
