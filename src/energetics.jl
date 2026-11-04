"""Nuclear energetics, specific energy generation rate ε_nuc (erg/g/s), and neutrino energy loss calculations."""

module Energetics

using ..Species
using ..Rates
using ..Networks

export compute_energy_generation, compute_neutrino_losses, net_thermal_energy_rate
export ERG_PER_MEV, AVOGADRO_NA

# Physical Constants (CODATA 2018)
const ERG_PER_MEV = 1.602176634e-6           # erg / MeV
const AVOGADRO_NA = 6.02214076e23            # mol^-1
const MEV_MOL_TO_ERG_G = AVOGADRO_NA * ERG_PER_MEV # ~ 9.6485e17 erg/(mol*MeV)

# Average neutrino energy loss per reaction in MeV (escape without depositing heat in standard stellar plasma)
const NEUTRINO_LOSS_MEV = Dict{String, Float64}(
    "pp"       => 0.265,   # Average pp-neutrino energy
    "n13_beta" => 0.707,   # 13N beta-decay neutrino
    "o15_beta" => 0.997,   # 15O beta-decay neutrino
    "be7p_b8"  => 6.711    # 8B high-energy solar neutrino
)

"""
    compute_energy_generation(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Compute the total specific nuclear energy release rate ε_nuc (in erg g^-1 s^-1).
Calculated via reaction fluxes: ε_nuc = N_A * c_unit * sum(r_k * Q_k).
"""
function compute_energy_generation(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    fluxes = Networks.compute_fluxes(net, Y, rho, T9)
    total_q_dot = 0.0 # MeV * mol / (g * s)

    for (k, rxn) in enumerate(net.reactions)
        total_q_dot += fluxes[k] * rxn.q_value_mev
    end

    eps_nuc = total_q_dot * MEV_MOL_TO_ERG_G
    return eps_nuc
end

"""
    compute_neutrino_losses(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Compute energy lost to escaping neutrinos ε_nu (in erg g^-1 s^-1).
"""
function compute_neutrino_losses(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    fluxes = Networks.compute_fluxes(net, Y, rho, T9)
    total_loss_dot = 0.0

    for (k, rxn) in enumerate(net.reactions)
        if haskey(NEUTRINO_LOSS_MEV, rxn.id)
            e_nu = NEUTRINO_LOSS_MEV[rxn.id]
            total_loss_dot += fluxes[k] * e_nu
        end
    end

    return total_loss_dot * MEV_MOL_TO_ERG_G
end

"""
    net_thermal_energy_rate(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Compute net thermal power deposited in stellar gas: ε_thermal = ε_nuc - ε_nu (in erg g^-1 s^-1).
"""
function net_thermal_energy_rate(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    eps_nuc = compute_energy_generation(net, Y, rho, T9)
    eps_nu = compute_neutrino_losses(net, Y, rho, T9)
    return max(0.0, eps_nuc - eps_nu)
end

end # module Energetics
