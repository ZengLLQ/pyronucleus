"""Network topologies, stoichiometric relations, reaction fluxes, and analytical Jacobian construction."""

module Networks

using ..Species
using ..Rates
using LinearAlgebra
using SparseArrays

export ReactionNetwork, compute_rhs!, compute_rhs, compute_jacobian!, compute_jacobian
export pp_chain_network, cno_cycle_network, triple_alpha_network, full_stellar_network

"""
    ReactionNetwork

A complete coupled thermonuclear reaction network containing species and reactions.

# Fields
- `species_net::NuclearNetwork`: Species definitions and index mappings.
- `reactions::Vector{Reaction}`: Included reaction channels.
- `stoichiometry::Matrix{Float64}`: Matrix of net stoichiometric changes c_{i,k} for species i in reaction k.
"""
struct ReactionNetwork
    species_net::NuclearNetwork
    reactions::Vector{Reaction}
    stoichiometry::Matrix{Float64} # (num_species, num_reactions)
    reactant_indices::Vector{Vector{Int}}
    product_indices::Vector{Vector{Int}}
    symbol_to_index::Dict{String, Int}

    function ReactionNetwork(species_symbols::Vector{String}, rxn_keys::Vector{String})
        sp_net = NuclearNetwork(species_symbols)
        rxns = [REACTION_DATABASE[k] for k in rxn_keys]
        n_sp = length(sp_net)
        n_rx = length(rxns)
        
        stoich = zeros(Float64, n_sp, n_rx)
        r_indices = Vector{Vector{Int}}(undef, n_rx)
        p_indices = Vector{Vector{Int}}(undef, n_rx)

        for (k, rxn) in enumerate(rxns)
            r_idx = Int[]
            for r in rxn.reactants
                if haskey(sp_net.symbol_to_index, r)
                    idx = sp_net.symbol_to_index[r]
                    push!(r_idx, idx)
                    stoich[idx, k] -= 1.0
                end
            end
            r_indices[k] = r_idx

            p_idx = Int[]
            for p in rxn.products
                if haskey(sp_net.symbol_to_index, p)
                    idx = sp_net.symbol_to_index[p]
                    push!(p_idx, idx)
                    stoich[idx, k] += 1.0
                end
            end
            p_indices[k] = p_idx
        end

        new(sp_net, rxns, stoich, r_indices, p_indices, sp_net.symbol_to_index)
    end
end

Base.length(net::ReactionNetwork) = length(net.species_net)

"""
    compute_fluxes(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Calculate reaction rates/fluxes r_k (in units of mol g^-1 s^-1) for all reactions.
Properly accounts for identical reactant symmetry factors (1/2! for 2 identical, 1/3! for 3 identical).
"""
function compute_fluxes(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    n_rx = length(net.reactions)
    fluxes = zeros(Float64, n_rx)
    
    for k in 1:n_rx
        rxn = net.reactions[k]
        raw_rate = evaluate_rate(rxn, T9)
        f_scr = screening_factor(rxn, rho, T9, Y, net.species_net)
        eff_rate = raw_rate * f_scr

        if rxn.type == Rates.DECAY
            idx1 = net.reactant_indices[k][1]
            fluxes[k] = eff_rate * max(0.0, Y[idx1])
        elseif rxn.type == Rates.BINARY
            idx1 = net.reactant_indices[k][1]
            idx2 = net.reactant_indices[k][2]
            y1 = max(0.0, Y[idx1])
            y2 = max(0.0, Y[idx2])
            
            if idx1 == idx2
                # Identical reactants: 1/2 * rho * rate * Y^2
                fluxes[k] = 0.5 * rho * eff_rate * (y1^2)
            else
                fluxes[k] = rho * eff_rate * y1 * y2
            end
        elseif rxn.type == Rates.TERTIARY
            # E.g. Triple-alpha: 1/6 * rho^2 * rate * Y^3
            idx1 = net.reactant_indices[k][1]
            y1 = max(0.0, Y[idx1])
            fluxes[k] = (1.0 / 6.0) * (rho^2) * eff_rate * (y1^3)
        end
    end
    
    return fluxes
end

"""
    compute_rhs!(dY::AbstractVector{Float64}, net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

In-place evaluation of the right-hand side dY/dt = S * r for abundance ODE integration.
"""
function compute_rhs!(dY::AbstractVector{Float64}, net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    fluxes = compute_fluxes(net, Y, rho, T9)
    mul!(dY, net.stoichiometry, fluxes)
    return dY
end

"""
    compute_rhs(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Out-of-place evaluation of the abundance derivatives dY/dt.
"""
function compute_rhs(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    dY = zeros(Float64, length(net))
    return compute_rhs!(dY, net, Y, rho, T9)
end

"""
    compute_jacobian!(J::AbstractMatrix{Float64}, net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Compute the exact analytical Jacobian matrix J_{i,j} = ∂(dY_i/dt) / ∂Y_j.
"""
function compute_jacobian!(J::AbstractMatrix{Float64}, net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    fill!(J, 0.0)
    n_sp = length(net)
    n_rx = length(net.reactions)

    for k in 1:n_rx
        rxn = net.reactions[k]
        raw_rate = evaluate_rate(rxn, T9)
        f_scr = screening_factor(rxn, rho, T9, Y, net.species_net)
        eff_rate = raw_rate * f_scr

        if rxn.type == Rates.DECAY
            j = net.reactant_indices[k][1]
            df_dyj = eff_rate
            for i in 1:n_sp
                c = net.stoichiometry[i, k]
                c != 0.0 && (J[i, j] += c * df_dyj)
            end
        elseif rxn.type == Rates.BINARY
            j1 = net.reactant_indices[k][1]
            j2 = net.reactant_indices[k][2]
            y1 = max(0.0, Y[j1])
            y2 = max(0.0, Y[j2])

            if j1 == j2
                # ∂(0.5 * rho * rate * Y1^2)/∂Y1 = rho * rate * Y1
                df_dy1 = rho * eff_rate * y1
                for i in 1:n_sp
                    c = net.stoichiometry[i, k]
                    c != 0.0 && (J[i, j1] += c * df_dy1)
                end
            else
                # ∂(rho * rate * Y1 * Y2)/∂Y1 = rho * rate * Y2
                df_dy1 = rho * eff_rate * y2
                df_dy2 = rho * eff_rate * y1
                for i in 1:n_sp
                    c = net.stoichiometry[i, k]
                    if c != 0.0
                        J[i, j1] += c * df_dy1
                        J[i, j2] += c * df_dy2
                    end
                end
            end
        elseif rxn.type == Rates.TERTIARY
            j1 = net.reactant_indices[k][1]
            y1 = max(0.0, Y[j1])
            # ∂((1/6) * rho^2 * rate * Y1^3)/∂Y1 = 0.5 * rho^2 * rate * Y1^2
            df_dy1 = 0.5 * (rho^2) * eff_rate * (y1^2)
            for i in 1:n_sp
                c = net.stoichiometry[i, k]
                c != 0.0 && (J[i, j1] += c * df_dy1)
            end
        end

        # Plasma electron screening abundance derivative: ∂(f_scr)/∂Y_j = f_scr * (H12 / (2*zeta)) * (Z_j^2 + Z_j)
        if rxn.type != Rates.DECAY
            iso1 = ISOTOPE_DATABASE[rxn.reactants[1]]
            iso2 = ISOTOPE_DATABASE[rxn.reactants[2]]
            Z1 = Float64(iso1.Z)
            Z2 = Float64(iso2.Z)
            if Z1 > 0.0 && Z2 > 0.0
                T6 = T9 * 1e3
                if T6 > 0.0
                    zeta = 0.0
                    for s in 1:n_sp
                        Zs = net.species_net.Z_vector[s]
                        zeta += (Zs^2 + Zs) * max(0.0, Y[s])
                    end
                    if zeta > 1e-10
                        H12 = 0.188 * Z1 * Z2 * sqrt(rho * zeta / (T6^3))
                        rxn.type == Rates.TERTIARY && (H12 *= 1.5)
                        if H12 < 12.0 # within uncapped screening regime
                            # Compute flux for this reaction
                            flux_k = 0.0
                            if rxn.type == Rates.BINARY
                                idx1 = net.reactant_indices[k][1]
                                idx2 = net.reactant_indices[k][2]
                                y1 = max(0.0, Y[idx1])
                                y2 = max(0.0, Y[idx2])
                                flux_k = idx1 == idx2 ? 0.5 * rho * eff_rate * (y1^2) : rho * eff_rate * y1 * y2
                            elseif rxn.type == Rates.TERTIARY
                                idx1 = net.reactant_indices[k][1]
                                y1 = max(0.0, Y[idx1])
                                flux_k = (1.0 / 6.0) * (rho^2) * eff_rate * (y1^3)
                            end

                            prefac = flux_k * (H12 / (2.0 * zeta))
                            for j in 1:n_sp
                                Zj = net.species_net.Z_vector[j]
                                if Zj > 0.0
                                    d_flux_scr = prefac * (Zj^2 + Zj)
                                    for i in 1:n_sp
                                        c = net.stoichiometry[i, k]
                                        c != 0.0 && (J[i, j] += c * d_flux_scr)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return J
end

"""
    compute_jacobian(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)

Out-of-place analytical Jacobian calculation.
"""
function compute_jacobian(net::ReactionNetwork, Y::AbstractVector{<:Real}, rho::Float64, T9::Float64)
    n = length(net)
    J = zeros(Float64, n, n)
    return compute_jacobian!(J, net, Y, rho, T9)
end

# --- Pre-configured Stellar Reaction Network Presets ---

"""
    pp_chain_network()

Stellar hydrogen burning via the Proton-Proton chains (PP-I, PP-II, PP-III).
Species: p, d, 3He, 4He, 7Li, 7Be, 8B
"""
function pp_chain_network()
    species = ["p", "d", "he3", "he4", "li7", "be7", "b8"]
    rxns = ["pp", "dp_he3", "he3he3_he4", "he3he4_be7", "be7p_b8", "li7p_2he4"]
    return ReactionNetwork(species, rxns)
end

"""
    cno_cycle_network()

Carbon-Nitrogen-Oxygen catalytic bi-cycle (CNO-I and CNO-II).
Species: p, 4He, 12C, 13C, 13N, 14N, 15N, 15O, 16O
"""
function cno_cycle_network()
    species = ["p", "he4", "c12", "c13", "n13", "n14", "n15", "o15", "o16"]
    rxns = ["c12p_n13", "n13_beta", "c13p_n14", "n14p_o15", "o15_beta", "n15p_c12he4", "n15p_o16"]
    return ReactionNetwork(species, rxns)
end

"""
    triple_alpha_network()

Helium burning via 3α and sequential alpha-capture ladder up to 28Si.
Species: 4He, 12C, 16O, 20Ne, 24Mg, 28Si
"""
function triple_alpha_network()
    species = ["he4", "c12", "o16", "ne20", "mg24", "si28"]
    rxns = ["triple_alpha", "c12alpha_o16", "o16alpha_ne20", "ne20alpha_mg24", "mg24alpha_si28"]
    return ReactionNetwork(species, rxns)
end

"""
    full_stellar_network()

Combined 22-isotope stellar reaction network spanning PP, CNO, and Alpha captures.
"""
function full_stellar_network()
    species = ["p", "d", "he3", "he4", "li7", "be7", "b8", 
               "c12", "c13", "n13", "n14", "n15", "o15", "o16", "ne20", "mg24", "si28"]
    rxns = ["pp", "dp_he3", "he3he3_he4", "he3he4_be7", "be7p_b8", "li7p_2he4",
            "c12p_n13", "n13_beta", "c13p_n14", "n14p_o15", "o15_beta", "n15p_c12he4", "n15p_o16",
            "triple_alpha", "c12alpha_o16", "o16alpha_ne20", "ne20alpha_mg24", "mg24alpha_si28"]
    return ReactionNetwork(species, rxns)
end

end # module Networks
