"""Isotopic properties, nuclear mass excess data, and composition conservation rules."""

module Species

using Printf

export Isotope, NuclearNetwork, ISOTOPE_DATABASE
export mass_fractions, electron_fraction, mean_molecular_weights, validate_conservation

"""
    Isotope

Nuclear isotope specification with atomic number, mass number, and binding energetics.

# Fields
- `name::String`: Full isotope designation (e.g., "Carbon-12").
- `symbol::String`: Canonical chemical symbol (e.g., "12C").
- `Z::Int`: Atomic number (number of protons).
- `A::Int`: Mass number (total nucleons).
- `N::Int`: Neutron number (A - Z).
- `mass_excess_mev::Float64`: Nuclear mass excess Δ in MeV (AME2020 evaluation).
- `spin::Float64`: Ground state nuclear spin J.
"""
struct Isotope
    name::String
    symbol::String
    Z::Int
    A::Int
    N::Int
    mass_excess_mev::Float64
    spin::Float64
end

function Isotope(name::String, symbol::String, Z::Int, A::Int, mass_excess::Float64, spin::Float64=0.0)
    return Isotope(name, symbol, Z, A, A - Z, mass_excess, spin)
end

# Accurate nuclear mass excess values from Atomic Mass Evaluation (AME2020) in MeV
const ISOTOPE_DATABASE = Dict{String, Isotope}(
    "n"    => Isotope("Neutron",       "n",    0,  1,  8.071317, 0.5),
    "p"    => Isotope("Proton",        "p",    1,  1,  7.288971, 0.5),
    "d"    => Isotope("Deuterium",     "2H",   1,  2, 13.135722, 1.0),
    "t"    => Isotope("Tritium",       "3H",   1,  3, 14.949810, 0.5),
    "he3"  => Isotope("Helium-3",      "3He",  2,  3, 14.931218, 0.5),
    "he4"  => Isotope("Helium-4",      "4He",  2,  4,  2.424916, 0.0),
    "li7"  => Isotope("Lithium-7",     "7Li",  3,  7, 14.907720, 1.5),
    "be7"  => Isotope("Beryllium-7",   "7Be",  4,  7, 15.769000, 1.5),
    "b8"   => Isotope("Boron-8",       "8B",   5,  8, 22.921600, 2.0),
    "c12"  => Isotope("Carbon-12",     "12C",  6, 12,  0.000000, 0.0), # Defined as 0 MeV reference
    "c13"  => Isotope("Carbon-13",     "13C",  6, 13,  3.125009, 0.5),
    "n13"  => Isotope("Nitrogen-13",   "13N",  7, 13,  5.345480, 0.5),
    "n14"  => Isotope("Nitrogen-14",   "14N",  7, 14,  2.863417, 1.0),
    "n15"  => Isotope("Nitrogen-15",   "15N",  7, 15,  0.101439, 0.5),
    "o15"  => Isotope("Oxygen-15",     "15O",  8, 15,  2.855600, 0.5),
    "o16"  => Isotope("Oxygen-16",     "16O",  8, 16, -4.737001, 0.0),
    "o17"  => Isotope("Oxygen-17",     "17O",  8, 17, -0.808764, 2.5),
    "f17"  => Isotope("Fluorine-17",   "17F",  9, 17,  1.951710, 2.5),
    "ne20" => Isotope("Neon-20",       "20Ne", 10, 20, -7.041930, 0.0),
    "mg24" => Isotope("Magnesium-24",  "24Mg", 12, 24,-13.933570, 0.0),
    "si28" => Isotope("Silicon-28",    "28Si", 14, 28,-21.492790, 0.0),
    "ni56" => Isotope("Nickel-56",     "56Ni", 28, 56,-39.407000, 0.0),
)

"""
    NuclearNetwork

Ordered set of nuclear species included in a reaction network.
"""
struct NuclearNetwork
    species::Vector{Isotope}
    symbol_to_index::Dict{String, Int}
    A_vector::Vector{Float64}
    Z_vector::Vector{Float64}
    mass_excess::Vector{Float64}

    function NuclearNetwork(symbols::Vector{String})
        species = [ISOTOPE_DATABASE[s] for s in symbols]
        sym_map = Dict(symbols[i] => i for i in 1:length(symbols))
        A_vec = [Float64(iso.A) for iso in species]
        Z_vec = [Float64(iso.Z) for iso in species]
        m_vec = [iso.mass_excess_mev for iso in species]
        new(species, sym_map, A_vec, Z_vec, m_vec)
    end
end

Base.length(net::NuclearNetwork) = length(net.species)

"""
    mass_fractions(net::NuclearNetwork, Y::AbstractVector{<:Real})

Compute mass fractions X_i = A_i * Y_i from molar abundances Y_i.
"""
function mass_fractions(net::NuclearNetwork, Y::AbstractVector{<:Real})
    return net.A_vector .* Y
end

"""
    electron_fraction(net::NuclearNetwork, Y::AbstractVector{<:Real})

Calculate total electron-to-nucleon ratio Y_e = sum(Z_i * Y_i).
"""
function electron_fraction(net::NuclearNetwork, Y::AbstractVector{<:Real})
    return sum(net.Z_vector .* Y)
end

"""
    mean_molecular_weights(net::NuclearNetwork, Y::AbstractVector{<:Real})

Compute ion mean molecular weight μ_ion and electron mean molecular weight μ_e.
"""
function mean_molecular_weights(net::NuclearNetwork, Y::AbstractVector{<:Real})
    sum_Y = sum(Y)
    y_e = electron_fraction(net, Y)
    mu_ion = sum_Y > 0 ? 1.0 / sum_Y : 1.0
    mu_e = y_e > 0 ? 1.0 / y_e : 1.0
    mu_total = 1.0 / (sum_Y + y_e)
    return mu_ion, mu_e, mu_total
end

"""
    validate_conservation(net::NuclearNetwork, Y::AbstractVector{<:Real}; tol=1e-5)

Ensure total mass fraction sum(A_i * Y_i) equals 1.0 and all abundances are non-negative.
"""
function validate_conservation(net::NuclearNetwork, Y::AbstractVector{<:Real}; tol=1e-5)
    any(Y .< -1e-12) && return false
    total_X = sum(mass_fractions(net, Y))
    return abs(total_X - 1.0) <= tol
end

end # module Species
