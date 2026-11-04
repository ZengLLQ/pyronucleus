"""Nuclear Statistical Equilibrium (NSE) solver via 2D Newton-Raphson on chemical potentials."""

module NSE

using ..Species
using LinearAlgebra
using Printf

export solve_nse, NSESolution, binding_energy_mev

"""
    NSESolution

Result of Nuclear Statistical Equilibrium calculation at given (rho, T9, Ye).

# Fields
- `Y::Vector{Float64}`: Equilibrium molar abundances Y_i.
- `X::Vector{Float64}`: Mass fractions X_i = A_i * Y_i.
- `Y_p::Float64`: Free proton abundance.
- `Y_n::Float64`: Free neutron abundance.
- `iterations::Int`: Number of Newton-Raphson iterations to reach convergence.
- `converged::Bool`: Convergence status flag.
"""
struct NSESolution
    Y::Vector{Float64}
    X::Vector{Float64}
    Y_p::Float64
    Y_n::Float64
    iterations::Int
    converged::Bool
end

# Physical constants for Saha-NSE equation
# k_B in MeV / K: k_B = 8.617333262e-11 MeV/K => k_B * 10^9 K = 0.08617333 MeV
const KB_T9_MEV = 0.08617333262 # MeV for T9 = 1.0

# Prefactor C_0: N_A * (2 * pi * hbar^2 / (m_u * k_B * 10^9 K))^(3/2) ~ 1.01333e-10 cm^3 / g
const SAHA_PREFACTOR = 1.01332933932e-10

"""
    binding_energy_mev(iso::Isotope)

Compute nuclear binding energy B = Z * Δ_p + N * Δ_n - Δ_i in MeV.
"""
function binding_energy_mev(iso::Isotope)
    delta_p = ISOTOPE_DATABASE["p"].mass_excess_mev
    delta_n = ISOTOPE_DATABASE["n"].mass_excess_mev
    return iso.Z * delta_p + iso.N * delta_n - iso.mass_excess_mev
end

"""
    compute_saha_coefficients(net::NuclearNetwork, rho::Float64, T9::Float64)

Compute pre-exponential Saha factors C_i in logarithmic space for Y_i = C_i * Y_p^Z * Y_n^N.
"""
function compute_saha_coefficients(net::NuclearNetwork, rho::Float64, T9::Float64)
    n = length(net)
    log_C = zeros(Float64, n)
    k_T = KB_T9_MEV * T9

    # Quantum concentration scale theta ~ rho * SAHA_PREFACTOR * T9^(-3/2)
    log_theta = log(max(rho, 1e-10) * SAHA_PREFACTOR) - 1.5 * log(T9)

    for i in 1:n
        iso = net.species[i]
        A = iso.A
        Z = iso.Z
        N = iso.N

        # Protons and neutrons are the reference species (log C = 0)
        if (Z == 1 && A == 1) || (Z == 0 && A == 1)
            log_C[i] = 0.0
            continue
        end

        spin = iso.spin
        G_i = 2.0 * spin + 1.0 # Nuclear spin degeneracy
        B_i = binding_energy_mev(iso)

        # log C_i = ln(G_i) + 1.5*ln(A) - A*ln(2) + (A-1)*log_theta + B_i / (k_B * T)
        log_val = log(G_i) + 1.5 * log(Float64(A)) - Float64(A) * log(2.0) + (Float64(A) - 1.0) * log_theta + (B_i / k_T)
        log_C[i] = log_val
    end

    return log_C
end

"""
    solve_nse(net::NuclearNetwork, rho::Float64, T9::Float64, Ye::Float64; 
              max_iters::Int=50, tol::Float64=1e-10)

Compute Nuclear Statistical Equilibrium composition for density `rho` (g/cm^3), 
temperature `T9` (10^9 K), and electron fraction `Ye`.
Solves nonlinear constraints sum(A_i * Y_i) = 1.0 and sum(Z_i * Y_i) = Ye via 2D Newton-Raphson.
"""
function solve_nse(
    net::NuclearNetwork,
    rho::Float64,
    T9::Float64,
    Ye::Float64;
    max_iters::Int=50,
    tol::Float64=1e-10
)
    log_C = compute_saha_coefficients(net, rho, T9)
    n = length(net)

    # Initial trial guess derived from dominant binding nucleus
    vals = [log_C[i] / Float64(net.species[i].A) for i in 1:n]
    dom_idx = argmax(vals)
    dom_iso = net.species[dom_idx]
    
    xi_0 = (-log(Float64(dom_iso.A)) - log_C[dom_idx]) / Float64(dom_iso.A)
    xi_p = xi_0 + log(clamp(2.0 * Ye, 0.01, 1.99))
    xi_n = xi_0 + log(clamp(2.0 * (1.0 - Ye), 0.01, 1.99))

    converged = false
    iter = 0

    Y = zeros(Float64, n)

    while iter < max_iters
        iter += 1

        # 1. Evaluate abundances Y_i from log abundances
        for i in 1:n
            iso = net.species[i]
            log_Yi = log_C[i] + Float64(iso.Z) * xi_p + Float64(iso.N) * xi_n
            Y[i] = log_Yi < -600.0 ? 0.0 : (log_Yi > 10.0 ? exp(10.0) : exp(log_Yi))
        end

        # 2. Residual vector: [sum(A_i * Y_i) - 1.0, sum(Z_i * Y_i) - Ye]
        f1 = sum(net.A_vector .* Y) - 1.0
        f2 = sum(net.Z_vector .* Y) - Ye

        if abs(f1) < tol && abs(f2) < tol
            converged = true
            break
        end

        # 3. 2x2 Jacobian: J_11 = ∂f1/∂xi_p = sum(A_i * Z_i * Y_i), etc.
        J11 = sum(net.A_vector .* net.Z_vector .* Y)
        J12 = sum(net.A_vector .* [Float64(s.N) for s in net.species] .* Y)
        J21 = sum(net.Z_vector .* net.Z_vector .* Y)
        J22 = sum(net.Z_vector .* [Float64(s.N) for s in net.species] .* Y)

        det_J = J11 * J22 - J12 * J21
        if abs(det_J) < 1e-30
            break
        end

        # Newton step: delta_xi = - J^-1 * f
        d_xi_p = -(J22 * f1 - J12 * f2) / det_J
        d_xi_n = -(-J21 * f1 + J11 * f2) / det_J

        step_norm = max(abs(d_xi_p), abs(d_xi_n))
        scale = step_norm > 2.0 ? 2.0 / step_norm : 1.0
        xi_p += scale * d_xi_p
        xi_n += scale * d_xi_n
    end

    # Enforce strict mass conservation on final result
    X = net.A_vector .* Y
    total_X = sum(X)
    if total_X > 0.0
        X ./= total_X
        Y .= X ./ net.A_vector
    end

    p_idx = get(net.symbol_to_index, "p", 0)
    n_idx = get(net.symbol_to_index, "n", 0)
    Y_p = p_idx > 0 ? Y[p_idx] : 0.0
    Y_n = n_idx > 0 ? Y[n_idx] : 0.0

    return NSESolution(Y, X, Y_p, Y_n, iter, converged)
end

end # module NSE
