"""
Pyronucleus: Thermonuclear reaction network and stellar nucleosynthesis simulation framework in Julia.
"""
module Pyronucleus

include("species.jl")
include("rates.jl")
include("networks.jl")
include("energetics.jl")
include("integrator.jl")
include("nse.jl")
include("cli.jl")

using .Species
using .Rates
using .Networks
using .Energetics
using .Integrator
using .NSE
using .CLI

# Re-export core APIs
export Isotope, NuclearNetwork, ISOTOPE_DATABASE
export mass_fractions, electron_fraction, mean_molecular_weights, validate_conservation

export Reaction, ReactionType, REACTION_DATABASE
export evaluate_rate, rate_temperature_derivative, screening_factor, q_value

export ReactionNetwork, compute_rhs!, compute_rhs, compute_jacobian!, compute_jacobian
export pp_chain_network, cno_cycle_network, triple_alpha_network, full_stellar_network

export compute_energy_generation, compute_neutrino_losses, net_thermal_energy_rate

export NetworkIntegrator, IntegrationResult, solve_burn, step_linearly_implicit!

export solve_nse, NSESolution, binding_energy_mev

export run_cli

end # module Pyronucleus
