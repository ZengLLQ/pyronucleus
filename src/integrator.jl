"""Stiff ordinary differential equation integrator for nuclear reaction networks."""

module Integrator

using ..Species
using ..Rates
using ..Networks
using ..Energetics
using LinearAlgebra

export NetworkIntegrator, IntegrationResult, solve_burn, step_linearly_implicit!

"""
    IntegrationResult

Trajectory solution of the nuclear reaction network integration.

# Fields
- `times::Vector{Float64}`: Time grid in seconds.
- `abundances::Matrix{Float64}`: Abundances Y_i at each time point (num_species × num_times).
- `energy_rates::Vector{Float64}`: Specific nuclear energy release rate ε_nuc (erg/g/s) over time.
- `temperatures::Vector{Float64}`: Temperature T9 trajectory.
- `densities::Vector{Float64}`: Density ρ trajectory in g/cm^3.
- `accepted_steps::Int`: Count of accepted solver steps.
- `rejected_steps::Int`: Count of rejected steps due to local error violation.
"""
struct IntegrationResult
    times::Vector{Float64}
    abundances::Matrix{Float64}
    energy_rates::Vector{Float64}
    temperatures::Vector{Float64}
    densities::Vector{Float64}
    accepted_steps::Int
    rejected_steps::Int
end

"""
    NetworkIntegrator

Configuration for adaptive stiff ODE integration.
"""
struct NetworkIntegrator
    rtol::Float64
    atol::Float64
    h_min::Float64
    h_max::Float64
    max_steps::Int
    project_mass::Bool

    function NetworkIntegrator(;
        rtol::Float64=1e-5,
        atol::Float64=1e-12,
        h_min::Float64=1e-25,
        h_max::Float64=3.15e16, # ~1 Gyr
        max_steps::Int=100000,
        project_mass::Bool=true
    )
        new(rtol, atol, h_min, h_max, max_steps, project_mass)
    end
end

"""
    step_linearly_implicit!(Y_next::Vector{Float64}, net::ReactionNetwork, Y::Vector{Float64}, 
                           rho::Float64, T9::Float64, h::Float64)

Perform a single linearly implicit Euler step:
(I - h*J) * ΔY = h * f(Y)
Y_{n+1} = Y_n + ΔY
"""
function step_linearly_implicit!(
    Y_next::Vector{Float64},
    net::ReactionNetwork,
    Y::Vector{Float64},
    rho::Float64,
    T9::Float64,
    h::Float64
)
    n = length(net)
    f = compute_rhs(net, Y, rho, T9)
    J = compute_jacobian(net, Y, rho, T9)

    # Construct system matrix: A = I - h * J
    A = Matrix{Float64}(I, n, n) - h * J
    
    # Solve linear system A * ΔY = h * f
    rhs = h * f
    delta_Y = A \ rhs

    for i in 1:n
        Y_next[i] = max(0.0, Y[i] + delta_Y[i])
    end

    return Y_next
end

"""
    project_mass_conservation!(net::ReactionNetwork, Y::Vector{Float64})

Renormalize abundances to strictly enforce sum(A_i * Y_i) == 1.0.
"""
function project_mass_conservation!(net::ReactionNetwork, Y::Vector{Float64})
    total_mass = sum(net.species_net.A_vector .* Y)
    if total_mass > 0.0 && abs(total_mass - 1.0) > 1e-15
        Y ./= total_mass
    end
    return Y
end

"""
    solve_burn(net::ReactionNetwork, Y0::Vector{Float64}, t_end::Float64; 
               rho::Float64, T9::Float64, integrator::NetworkIntegrator=NetworkIntegrator(),
               dt_init::Float64=1e-8)

Integrate nuclear burning from t=0 to t=t_end under conditions (rho, T9).
Uses adaptive step size semi-implicit Euler with Richardson extrapolation and mass projection.
"""
function solve_burn(
    net::ReactionNetwork,
    Y0::Vector{Float64},
    t_end::Float64;
    rho::Float64,
    T9::Float64,
    integrator::NetworkIntegrator=NetworkIntegrator(),
    dt_init::Float64=1e-8
)
    n = length(net)
    Y_curr = copy(Y0)
    integrator.project_mass && project_mass_conservation!(net, Y_curr)

    t = 0.0
    h = min(dt_init, t_end)

    # Output storage
    times = Float64[0.0]
    abundances_list = Vector{Float64}[copy(Y_curr)]
    energy_rates = Float64[compute_energy_generation(net, Y_curr, rho, T9)]
    temperatures = Float64[T9]
    densities = Float64[rho]

    accepted = 0
    rejected = 0
    step = 0

    # Temp buffers
    Y_full = zeros(Float64, n)
    Y_half1 = zeros(Float64, n)
    Y_half2 = zeros(Float64, n)

    while t < t_end && step < integrator.max_steps
        step += 1

        # Prevent stepping past t_end
        if t + h > t_end
            h = t_end - t
        end

        # 1. Full step of size h
        step_linearly_implicit!(Y_full, net, Y_curr, rho, T9, h)

        # 2. Two half steps of size h/2
        step_linearly_implicit!(Y_half1, net, Y_curr, rho, T9, 0.5 * h)
        step_linearly_implicit!(Y_half2, net, Y_half1, rho, T9, 0.5 * h)

        # 3. Local error estimation (compare full step to double-half step)
        err = 0.0
        for i in 1:n
            scale = integrator.atol + integrator.rtol * max(abs(Y_curr[i]), abs(Y_half2[i]))
            diff = abs(Y_half2[i] - Y_full[i])
            err = max(err, diff / scale)
        end

        if err <= 1.0 || h <= integrator.h_min
            # Accept step
            accepted += 1
            t += h

            # Richardson extrapolation for 2nd order accuracy: Y = 2*Y_half2 - Y_full
            for i in 1:n
                Y_curr[i] = max(0.0, 2.0 * Y_half2[i] - Y_full[i])
            end
            integrator.project_mass && project_mass_conservation!(net, Y_curr)

            push!(times, t)
            push!(abundances_list, copy(Y_curr))
            push!(energy_rates, compute_energy_generation(net, Y_curr, rho, T9))
            push!(temperatures, T9)
            push!(densities, rho)

            # Propose larger step
            factor = min(2.0, max(0.2, 0.9 * (max(err, 1e-4)^(-0.5))))
            h = clamp(h * factor, integrator.h_min, integrator.h_max)
        else
            # Reject step
            rejected += 1
            factor = max(0.1, 0.9 * (err^(-0.5)))
            h = max(integrator.h_min, h * factor)
        end
    end

    # Format abundance matrix (species × times)
    n_pts = length(times)
    abundances_matrix = zeros(Float64, n, n_pts)
    for j in 1:n_pts
        abundances_matrix[:, j] .= abundances_list[j]
    end

    return IntegrationResult(times, abundances_matrix, energy_rates, temperatures, densities, accepted, rejected)
end

end # module Integrator
