module NaturalUnits
# Natural units (https://en.wikipedia.org/wiki/Natural_units) signify
# different systems of units, depending on the application. Here is
# implemented a way to easily define a system of units by providing
# expressions for how the basic units of the system are constructed
# from natural constants. Also provided are implementations for four
# common systems, Hartree atomic units (ℏ = e = mₑ = 1/(4πε₀) = kB =
# 1), Rydberg atomic units (ℏ = e/√2 = 2mₑ = 1/(4πε₀) = kB = 1), both
# used in atomic physics, "natural units" (ℏ = c = kB = ε₀ = μ₀ = eV =
# 1) used in particle physics and cosmology, and geometrized units (c
# = G = kB = ε₀ = 1), used in general relativity.
#
# Overview of the algorithm to transfer from one system of units, to
# another:
#
# SI units (Système International d'Unités) consist of the following
# base units: m, kg, s, A, K, mol, and cd. The two last units don't
# really have a use in natural units, so they will be excluded from
# the following, and we are left with the set {uᵢ} = {m, kg, s, A, K}.
# Let y be the unit of any quantity. We can then express y as
#
# y = ∏ uᵢ^{yᵢ}, yᵢ ∈ ℚ,
#     ᵢ
#
# e.g., W/m² = kg·s⁻³. We define the vector 𝐲 = [yᵢ]. For the case of
# W/cm², we would have 𝐲 = [0 1 -3 0 0]ᵀ.
#
# For the case of Hartree atomic units, we have the following
# constants normalized: ℏ = mₑ = e = 1/(4πε₀) = kᴮ = 1. If we label
# the units of these quantities as {𝐚ⱼ}, we get
#
# 𝐚₁ = Js  = [ 2  1 -1  0  0]ᵀ
# 𝐚₂ = kg  = [ 0  1  0  0  0]ᵀ
# 𝐚₃ = C   = [ 0  0  1  1  0]ᵀ
# 𝐚₄ = m/F = [ 3  1 -4 -2  0]ᵀ
# 𝐚₅ = J/K = [ 2  1 -2  0 -1]ᵀ
#
# 𝐚ⱼ ∈ ℚ⁵
#
# We now want to express the unit y in the system spanned by {𝐚ⱼ}:
#
# y = ∏ aⱼ^{pⱼ}, pᵢ ∈ ℚ.
#     ⱼ
#
# We find the vector 𝐩 = [pⱼ] by solving the equation system
#
# 𝖠𝐩 = 𝐲,
#
# where 𝖠 = [𝐚ⱼ]. Since 𝖠 ∈ ℚ⁵ˣ⁵, we can't use inv(A) to solve the
# system, since that implicitly converts to Float64. However,
# factorize(A) returns a rational LU factorization, so we may
# calculate the inverse as Ainv = factorize(A)\eye(A).
#
# When we have solved for 𝐩, we know the unit of y in the target unit
# system. To calculate the conversion factor, we perform the dot
# product c = 𝐚ᵀ𝐩, and quantity y in the target unit system is given
# by y/c. We retain the composite SI unit for possible back-conversion
# later.
using Physical
using Quantities

NaturalUnitSystems = Dict{Symbol,Type}()

import Quantities.*, Quantities./

# {m, kg, s, A, K}
si_units = map(u -> first(keys(u.unit.d)),
               [Meter, Kilogram, Second, Ampere, Kelvin])
gram_unit = first(keys(Gram.unit.d))

# This function takes a quantity, and returns its value in the base
# units of SI, {m, kg, s, A, K} along with the powers of these base
# units. The powers are expected to be rational.
function si_powers(q)
    si = asbase(q)
    val = si.value
    unit_powers =
    map(function (u)
        k = keys(si.unit.d)
        convert(Rational{Int},
                if u in k
                si.unit.d[u]
                elseif u == si_units[2] && gram_unit in k
                val /= 1000
                si.unit.d[gram_unit]
                else
                0
                end)
    end, si_units)
    val,unit_powers
end

# Given a set of unit, construct the system matrix 𝖠 = [𝐚ⱼ]
function unit_system_matrix(units)
    units_si = map(si_powers, units)
    A = hcat(map(c -> c[2], units_si)...)
    A,factorize(A)\eye(A)
end

macro define_natural_unit_system(name, basic_unit, units)
    quantity_sym = Symbol("$(name)UnitsQuantity")
    units_sym = Symbol("$(name)Units")
    matrix_sym = Symbol("$(name)Matrix")
    inv_matrix_sym = Symbol("$(name)InvMatrix")
    units_v = eval(units)
    name_str = string(name)
    conv_fun_sym = Symbol("conv_fac_$(name)")

    @eval begin
        $units_sym = $units_v
        $matrix_sym,$inv_matrix_sym = unit_system_matrix($units_v)

        type $quantity_sym{T<:Quantities.QValue}
            value::T # Value in natural unit system
            unit::Quantities.Quantity{T} # Conversion factor, in SI units
        end

        import Base.*
        # Multiplication by scalar
        function *{T<:Quantities.QValue,U<:Quantities.QValue}(a::T, b::$quantity_sym{U})
            V = promote_type(T,U)
            $quantity_sym{V}(a*b.value, convert(Quantities.Quantity{V}, b.unit))
        end
        function *{T<:Quantities.QValue,U<:Quantities.QValue}(a::$quantity_sym{T}, b::U)
            V = promote_type(T,U)
            $quantity_sym{V}(a.value*b, convert(Quantities.Quantity{V}, a.unit))
        end
        # Multiplication by other quantity
        function *{T<:Quantities.QValue,U<:Quantities.QValue}(a::$quantity_sym{T}, b::$quantity_sym{U})
            V = promote_type(T,U)
            $quantity_sym{V}(a.value*b.value, convert(Quantities.Quantity{V}, a.unit*b.unit))
        end

        import Base./
        # Division by scalar
        function /{T<:Quantities.QValue,U<:Quantities.QValue}(a::T, b::$quantity_sym{U})
            V = promote_type(T,U)
            $quantity_sym{V}(a/b.value, convert(Quantities.Quantity{V}, b.unit^(-1)))
        end
        function /{T<:Quantities.QValue,U<:Quantities.QValue}(a::$quantity_sym{T}, b::U)
            V = promote_type(T,U)
            $quantity_sym{V}(a.value/b, convert(Quantities.Quantity{V}, a.unit))
        end
        # Division by other quantity
        function /{T<:Quantities.QValue,U<:Quantities.QValue}(a::$quantity_sym{T}, b::$quantity_sym{U})
            V = promote_type(T,U)
            $quantity_sym{V}(a.value/b.value, convert(Quantities.Quantity{V}, a.unit/b.unit))
        end

        import Base.+
        import Base.-

        function +{T<:Quantities.QValue,U<:Quantities.QValue}(a::$quantity_sym{T}, b::$quantity_sym{U})
            V = promote_type(T,U)
            assert(a.unit == b.unit)
            $quantity_sym{V}(a.value+b.value, a.unit)
        end

        import Base.√
        function √{T<:Quantities.QValue}(a::$quantity_sym{T})
            $quantity_sym{T}(√(a.value), √(a.unit))
        end

        function $conv_fun_sym{T}(si::Quantities.Quantity{T})
            units_si = si_powers(si)
            powers = $inv_matrix_sym*units_si[2] # 𝐩 = 𝖠⁻¹𝐲
            reduce(*, map(a -> a[1]^a[2], zip($units_v,powers)))
        end

        import Base.convert
        function convert{T<:Quantities.QValue,U<:Quantities.QValue}(to::Type{$quantity_sym{U}}, from::Quantities.Quantity{T})
            conv_fac = $conv_fun_sym(from)
            $quantity_sym{U}(from/conv_fac,conv_fac)
        end

        # For idempotence
        function convert{T<:Quantities.QValue,U<:Quantities.QValue}(to::Type{$quantity_sym{U}}, from::$quantity_sym{T})
            $quantity_sym{U}(from.value,from.unit)
        end

        function si{T<:Quantities.QValue}(from::$quantity_sym{T})
            from.value*from.unit
        end

        export $quantity_sym, $units_sym, *, √, as
        NaturalUnitSystems[Symbol($name_str)] = $quantity_sym
    end
end

function get_natural_unity(name, T = Float64)
    NaturalUnitSystems[name]{T}
end

@define_natural_unit_system(Hartree, "au",
                            [Hbar_planck,
                             E_electron,
                             M_electron,
                             K_coulomb,
                             K_boltzmann])
@define_natural_unit_system(Rydberg, "au",
                            [Hbar_planck,
                             E_electron/√2,
                             2M_electron,
                             K_coulomb,
                             K_boltzmann])

@define_natural_unit_system(Planck, "planck",
                            [G_gravitation,
                             Hbar_planck,
                             C_light,
                             K_coulomb,
                             K_boltzmann])

# # Overdetermined
# @define_natural_unit_system(Natural, ElectronVolt,
#                             [Hbar_planck,
#                              C_light,
#                              K_boltzmann,
#                              E0_electric,
#                              U0_magnetic,
#                              ElectronVolt])

# # Underdetermined
# @define_natural_unit_system(Geometrized, Meter,
#                             [C_light,
#                              G_gravitation,
#                              K_boltzmann,
#                              E0_electric])

# # Convenience alias
# const AtomicUnitsOne = get_natural_unity(:Hartree)
# export AtomicUnitsOne

export @define_natural_unit_system, get_natural_unity
end
