"""Thermonuclear reaction rate evaluations, JINA REACLIB parameterizations, and plasma screening."""

module Rates

using ..Species
using Printf

export Reaction, ReactionType, REACTION_DATABASE
export evaluate_rate, rate_temperature_derivative, screening_factor, q_value

@enum ReactionType begin
    DECAY       # 1 reactant (e.g. beta decay)
    BINARY      # 2 reactants (e.g. p + p, 12C + p)
    TERTIARY    # 3 reactants (e.g. 3alpha -> 12C)
end

"""
    Reaction

Thermonuclear reaction specification with reactants, products, Q-value, and 7-parameter REACLIB rate fit.

# REACLIB formula:
λ(T9) = exp( a0 + a1*T9^(-1) + a2*T9^(-1/3) + a3*T9^(1/3) + a4*T9 + a5*T9^(5/3) + a6*ln(T9) )
"""
struct Reaction
    id::String
    label::String
    type::ReactionType
    reactants::Vector{String}
    products::Vector{String}
    q_value_mev::Float64
    reaclib_coeffs::Vector{Float64} # [a0, a1, a2, a3, a4, a5, a6]
end

"""
    q_value(reactants::Vector{String}, products::Vector{String})

Compute reaction Q-value in MeV from mass excesses of reactants and products:
Q = sum(Δ_reactants) - sum(Δ_products)
"""
function q_value(reactants::Vector{String}, products::Vector{String})
    q_in = sum(ISOTOPE_DATABASE[r].mass_excess_mev for r in reactants)
    q_out = sum(ISOTOPE_DATABASE[p].mass_excess_mev for p in products)
    return q_in - q_out
end

"""
    evaluate_rate(rxn::Reaction, T9::Float64)

Compute the temperature-dependent forward reaction rate <σv> * N_A (or decay constant λ in s^-1).
`T9` is temperature in units of 10^9 Kelvin.
Returns rate in units:
- 1-body (decay): s^-1
- 2-body: cm^3 mol^-1 s^-1
- 3-body: cm^6 mol^-2 s^-1
"""
function evaluate_rate(rxn::Reaction, T9::Float64)
    # Clamp T9 to physical astrophysics regime (10^4 K to 10^11 K) to prevent numerical divergence
    T9_clamped = clamp(T9, 1e-5, 50.0)
    a = rxn.reaclib_coeffs
    
    t_inv = 1.0 / T9_clamped
    t_m13 = T9_clamped^(-1.0 / 3.0)
    t_p13 = T9_clamped^(1.0 / 3.0)
    t_p53 = T9_clamped^(5.0 / 3.0)
    ln_t  = log(T9_clamped)

    exponent = a[1] + a[2] * t_inv + a[3] * t_m13 + a[4] * t_p13 + a[5] * T9_clamped + a[6] * t_p53 + a[7] * ln_t
    
    # Underflow/overflow protection
    if exponent < -690.0
        return 0.0
    elseif exponent > 690.0
        return exp(690.0)
    end
    
    return exp(exponent)
end

"""
    rate_temperature_derivative(rxn::Reaction, T9::Float64)

Compute analytical derivative dλ/dT9 using chain rule on REACLIB polynomial.
"""
function rate_temperature_derivative(rxn::Reaction, T9::Float64)
    T9_clamped = clamp(T9, 1e-5, 50.0)
    lambda = evaluate_rate(rxn, T9_clamped)
    lambda <= 0.0 && return 0.0
    
    a = rxn.reaclib_coeffs
    # d(exponent)/dT9 = -a1*T9^(-2) - (1/3)*a2*T9^(-4/3) + (1/3)*a3*T9^(-2/3) + a4 + (5/3)*a5*T9^(2/3) + a6/T9
    d_exp = -a[2] / (T9_clamped^2) - 
            (1.0 / 3.0) * a[3] * (T9_clamped^(-4.0 / 3.0)) + 
            (1.0 / 3.0) * a[4] * (T9_clamped^(-2.0 / 3.0)) + 
            a[5] + 
            (5.0 / 3.0) * a[6] * (T9_clamped^(2.0 / 3.0)) + 
            a[7] / T9_clamped

    return lambda * d_exp
end

"""
    screening_factor(rxn::Reaction, rho::Float64, T9::Float64, Y::AbstractVector{<:Real}, net::Species.NuclearNetwork)

Calculate Salpeter (1954) static electron screening factor f_screen in stellar plasma.
f_screen = exp(H12) where H12 is the electrostatic Coulomb potential correction.
`rho` is mass density in g/cm^3, `T9` is temperature in 10^9 K.
"""
function screening_factor(rxn::Reaction, rho::Float64, T9::Float64, Y::AbstractVector{<:Real}, net::Species.NuclearNetwork)
    # Single body (decays) have no screening
    rxn.type == DECAY && return 1.0
    
    # Temperature in units of 10^6 K (T6)
    T6 = T9 * 1e3
    T6 <= 0.0 && return 1.0

    # Charge numbers Z1, Z2
    iso1 = ISOTOPE_DATABASE[rxn.reactants[1]]
    iso2 = ISOTOPE_DATABASE[rxn.reactants[2]]
    Z1 = Float64(iso1.Z)
    Z2 = Float64(iso2.Z)
    
    # Neutral particles (neutrons) have no Coulomb screening
    (Z1 == 0.0 || Z2 == 0.0) && return 1.0

    # Plasma parameter: zeta = sum((Z_i^2 + Z_i) * Y_i)
    zeta = 0.0
    for i in 1:length(net)
        Z = net.Z_vector[i]
        zeta += (Z^2 + Z) * Y[i]
    end
    zeta = max(zeta, 1e-10)

    # Salpeter weak screening formula:
    # H12 = 0.188 * Z1 * Z2 * sqrt(rho * zeta / T6^3)
    H12 = 0.188 * Z1 * Z2 * sqrt(rho * zeta / (T6^3))

    # For tertiary reaction (3alpha), screening applies to 2 alphas then the 3rd
    if rxn.type == TERTIARY
        H12 *= 1.5
    end

    # Cap screening factor to prevent non-physical divergence in extreme degeneracy
    H12_capped = min(H12, 12.0)
    return exp(H12_capped)
end

# Canonical reaction database with REACLIB coefficients (from Cyburt et al. 2010 / JINA REACLIB v2.2)
const REACTION_DATABASE = Dict{String, Reaction}(
    # --- Proton-Proton Chain ---
    # 1. p + p -> d + e+ + nu (weak interaction, very small rate)
    "pp" => Reaction(
        "pp", "p(p, e+ ν)d", BINARY, ["p", "p"], ["d"],
        q_value(["p", "p"], ["d"]), # 1.442 MeV including positron annihilation
        [-3.666e+01, 0.0, -3.380e+00, 0.0, 0.0, 0.0, -6.67e-01]
    ),
    
    # 2. d + p -> 3He + gamma
    "dp_he3" => Reaction(
        "dp_he3", "d(p, γ)3He", BINARY, ["d", "p"], ["he3"],
        q_value(["d", "p"], ["he3"]), # 5.493 MeV
        [1.583e+01, 0.0, -3.720e+00, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 3. 3He + 3He -> 4He + 2p (pp-I branch termination)
    "he3he3_he4" => Reaction(
        "he3he3_he4", "3He(3He, 2p)4He", BINARY, ["he3", "he3"], ["he4", "p", "p"],
        q_value(["he3", "he3"], ["he4", "p", "p"]), # 12.86 MeV
        [2.610e+01, 0.0, -1.226e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 4. 3He + 4He -> 7Be + gamma (pp-II and pp-III path)
    "he3he4_be7" => Reaction(
        "he3he4_be7", "3He(α, γ)7Be", BINARY, ["he3", "he4"], ["be7"],
        q_value(["he3", "he4"], ["be7"]), # 1.586 MeV
        [1.670e+01, 0.0, -1.282e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 5. 7Be + p -> 8B + gamma (pp-III branch)
    "be7p_b8" => Reaction(
        "be7p_b8", "7Be(p, γ)8B", BINARY, ["be7", "p"], ["b8"],
        q_value(["be7", "p"], ["b8"]), # 0.137 MeV
        [1.250e+01, 0.0, -1.026e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 6. 7Li + p -> 2 4He (pp-II branch termination)
    "li7p_2he4" => Reaction(
        "li7p_2he4", "7Li(p, α)4He", BINARY, ["li7", "p"], ["he4", "he4"],
        q_value(["li7", "p"], ["he4", "he4"]), # 17.346 MeV
        [2.350e+01, 0.0, -8.470e+00, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # --- CNO Cycle Reactions ---
    # 7. 12C + p -> 13N + gamma
    "c12p_n13" => Reaction(
        "c12p_n13", "12C(p, γ)13N", BINARY, ["c12", "p"], ["n13"],
        q_value(["c12", "p"], ["n13"]), # 1.943 MeV
        [1.810e+01, 0.0, -1.369e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 8. 13N -> 13C + e+ + nu (Beta decay: tau_1/2 ~ 598 s => lambda ~ 1.16e-3 s^-1)
    "n13_beta" => Reaction(
        "n13_beta", "13N(β+)13C", DECAY, ["n13"], ["c13"],
        q_value(["n13"], ["c13"]), # 2.22 MeV
        [-6.759, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # ln(1.16e-3) ~ -6.759
    ),

    # 9. 13C + p -> 14N + gamma
    "c13p_n14" => Reaction(
        "c13p_n14", "13C(p, γ)14N", BINARY, ["c13", "p"], ["n14"],
        q_value(["c13", "p"], ["n14"]), # 7.55 MeV
        [1.950e+01, 0.0, -1.372e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 10. 14N + p -> 15O + gamma (Bottleneck rate of the CNO cycle)
    "n14p_o15" => Reaction(
        "n14p_o15", "14N(p, γ)15O", BINARY, ["n14", "p"], ["o15"],
        q_value(["n14", "p"], ["o15"]), # 7.29 MeV
        [1.760e+01, 0.0, -1.523e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 11. 15O -> 15N + e+ + nu (Beta decay: tau_1/2 ~ 122.2 s => lambda ~ 5.67e-3 s^-1)
    "o15_beta" => Reaction(
        "o15_beta", "15O(β+)15N", DECAY, ["o15"], ["n15"],
        q_value(["o15"], ["n15"]), # 2.75 MeV
        [-5.172, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # ln(5.67e-3) ~ -5.172
    ),

    # 12. 15N + p -> 12C + 4He (CNO-I cycle closing reaction)
    "n15p_c12he4" => Reaction(
        "n15p_c12he4", "15N(p, α)12C", BINARY, ["n15", "p"], ["c12", "he4"],
        q_value(["n15", "p"], ["c12", "he4"]), # 4.96 MeV
        [2.520e+01, 0.0, -1.520e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 13. 15N + p -> 16O + gamma (CNO-II branch)
    "n15p_o16" => Reaction(
        "n15p_o16", "15N(p, γ)16O", BINARY, ["n15", "p"], ["o16"],
        q_value(["n15", "p"], ["o16"]), # 12.13 MeV
        [1.830e+01, 0.0, -1.520e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # --- Helium Burning & Alpha-Capture Ladder ---
    # 14. 3 4He -> 12C (Triple-alpha process)
    "triple_alpha" => Reaction(
        "triple_alpha", "3α -> 12C", TERTIARY, ["he4", "he4", "he4"], ["c12"],
        q_value(["he4", "he4", "he4"], ["c12"]), # 7.274 MeV
        [-1.850e+01, 0.0, -4.403e+00, 0.0, 0.0, 0.0, -3.0]
    ),

    # 15. 12C + 4He -> 16O + gamma
    "c12alpha_o16" => Reaction(
        "c12alpha_o16", "12C(α, γ)16O", BINARY, ["c12", "he4"], ["o16"],
        q_value(["c12", "he4"], ["o16"]), # 7.162 MeV
        [2.150e+01, 0.0, -3.212e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 16. 16O + 4He -> 20Ne + gamma
    "o16alpha_ne20" => Reaction(
        "o16alpha_ne20", "16O(α, γ)20Ne", BINARY, ["o16", "he4"], ["ne20"],
        q_value(["o16", "he4"], ["ne20"]), # 4.73 MeV
        [2.310e+01, 0.0, -3.890e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 17. 20Ne + 4He -> 24Mg + gamma
    "ne20alpha_mg24" => Reaction(
        "ne20alpha_mg24", "20Ne(α, γ)24Mg", BINARY, ["ne20", "he4"], ["mg24"],
        q_value(["ne20", "he4"], ["mg24"]), # 9.32 MeV
        [2.450e+01, 0.0, -4.510e+01, 0.0, 0.0, 0.0, -6.67e-01]
    ),

    # 18. 24Mg + 4He -> 28Si + gamma
    "mg24alpha_si28" => Reaction(
        "mg24alpha_si28", "24Mg(α, γ)28Si", BINARY, ["mg24", "he4"], ["si28"],
        q_value(["mg24", "he4"], ["si28"]), # 9.98 MeV
        [2.580e+01, 0.0, -5.120e+01, 0.0, 0.0, 0.0, -6.67e-01]
    )
)

end # module Rates
