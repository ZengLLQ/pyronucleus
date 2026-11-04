"""Command-line interface for Pyronucleus stellar reaction simulations."""

module CLI

using ..Species
using ..Rates
using ..Networks
using ..Energetics
using ..Integrator
using ..NSE
using Printf

export run_cli

function print_header()
    println("================================================================================")
    println("  PYRONUCLEUS: Thermonuclear Reaction Networks & Stellar Nucleosynthesis")
    println("  Stiff ODE Solver | JINA REACLIB | Salpeter Screening | Saha NSE Equilibrium")
    println("================================================================================")
end

function cmd_solar_core(args::Vector{String})
    print_header()
    println("Mode: Solar Core Hydrogen Burning (PP-Chain Dominant)")
    
    # Standard Solar Model (SSM) core conditions
    rho = 150.0       # g / cm^3
    T9  = 0.0157      # 1.57e7 K
    t_end = 3.1536e13 # ~ 1 million years (in seconds)

    net = pp_chain_network()
    println(@sprintf("Network: %d species, %d reaction channels", length(net), length(net.reactions)))
    println(@sprintf("Conditions: Density ρ = %.1f g/cm³, Temperature T = %.3f × 10⁶ K (T9 = %.4f)", rho, T9 * 1e3, T9))

    # Initial Solar composition: X(1H)=0.73, X(4He)=0.25, traces
    Y0 = zeros(Float64, length(net))
    Y0[net.symbol_to_index["p"]]   = 0.73 / 1.0
    Y0[net.symbol_to_index["he4"]] = 0.25 / 4.0
    Y0[net.symbol_to_index["he3"]] = 1e-5 / 3.0
    Y0[net.symbol_to_index["d"]]   = 1e-10 / 2.0

    integ = NetworkIntegrator(rtol=1e-5, atol=1e-12)
    println("\nIntegrating stiff thermonuclear network...")
    res = solve_burn(net, Y0, t_end; rho=rho, T9=T9, integrator=integ, dt_init=1.0)

    println(@sprintf("Completed in %d accepted steps (%d rejected steps)", res.accepted_steps, res.rejected_steps))
    println("\nFinal Composition after 1 Million Years:")
    println("--------------------------------------------------------------------------------")
    println(" Isotope | Mass Number A | Molar Abundance Y_i | Mass Fraction X_i")
    println("--------------------------------------------------------------------------------")
    final_Y = res.abundances[:, end]
    for (i, iso) in enumerate(net.species_net.species)
        X_i = iso.A * final_Y[i]
        println(@sprintf(" %-7s | %13d | %19.6e | %17.6e", iso.symbol, iso.A, final_Y[i], X_i))
    end
    println("--------------------------------------------------------------------------------")
    
    eps_nuc = compute_energy_generation(net, final_Y, rho, T9)
    eps_nu  = compute_neutrino_losses(net, final_Y, rho, T9)
    eps_th  = net_thermal_energy_rate(net, final_Y, rho, T9)

    println(@sprintf("Specific Nuclear Power (ε_nuc):     %12.4e erg g⁻¹ s⁻¹", eps_nuc))
    println(@sprintf("Solar Neutrino Loss Rate (ε_ν):     %12.4e erg g⁻¹ s⁻¹", eps_nu))
    println(@sprintf("Net Thermal Gas Deposition (ε_th):  %12.4e erg g⁻¹ s⁻¹", eps_th))
    println(@sprintf("Solar Core Luminosity density:       %12.4e erg cm⁻³ s⁻¹", eps_th * rho))
end

function cmd_cno_burn(args::Vector{String})
    print_header()
    println("Mode: Massive Star Core Hydrogen Burning (CNO Bi-Cycle)")

    # Massive main sequence star conditions (M > 1.5 M_sun)
    rho = 50.0        # g / cm^3
    T9  = 0.025       # 25 million K
    t_end = 3.1536e11 # ~ 10,000 years

    net = cno_cycle_network()
    println(@sprintf("Network: %d species, %d reaction channels", length(net), length(net.reactions)))
    println(@sprintf("Conditions: Density ρ = %.1f g/cm³, Temperature T = %.1f × 10⁶ K (T9 = %.4f)", rho, T9 * 1e3, T9))

    # Initial composition: X(H) = 0.70, X(He) = 0.28, X(12C) = 0.01, X(16O) = 0.01
    Y0 = zeros(Float64, length(net))
    Y0[net.symbol_to_index["p"]]   = 0.70 / 1.0
    Y0[net.symbol_to_index["he4"]] = 0.28 / 4.0
    Y0[net.symbol_to_index["c12"]] = 0.01 / 12.0
    Y0[net.symbol_to_index["o16"]] = 0.01 / 16.0

    integ = NetworkIntegrator(rtol=1e-5, atol=1e-12)
    println("\nIntegrating CNO catalytic cycle...")
    res = solve_burn(net, Y0, t_end; rho=rho, T9=T9, integrator=integ, dt_init=1e-2)

    println(@sprintf("Completed in %d accepted steps (%d rejected steps)", res.accepted_steps, res.rejected_steps))
    println("\nEquilibrated CNO Abundances at 10,000 Years:")
    println("--------------------------------------------------------------------------------")
    println(" Isotope | Mass Number A | Molar Abundance Y_i | Mass Fraction X_i")
    println("--------------------------------------------------------------------------------")
    final_Y = res.abundances[:, end]
    for (i, iso) in enumerate(net.species_net.species)
        X_i = iso.A * final_Y[i]
        println(@sprintf(" %-7s | %13d | %19.6e | %17.6e", iso.symbol, iso.A, final_Y[i], X_i))
    end
    println("--------------------------------------------------------------------------------")

    eps_nuc = compute_energy_generation(net, final_Y, rho, T9)
    eps_nu  = compute_neutrino_losses(net, final_Y, rho, T9)
    eps_th  = net_thermal_energy_rate(net, final_Y, rho, T9)

    println(@sprintf("Specific Nuclear Power (ε_nuc):     %12.4e erg g⁻¹ s⁻¹", eps_nuc))
    println(@sprintf("CNO Neutrino Loss Rate (ε_ν):       %12.4e erg g⁻¹ s⁻¹", eps_nu))
    println(@sprintf("Net Thermal Gas Deposition (ε_th):  %12.4e erg g⁻¹ s⁻¹", eps_th))
end

function cmd_triple_alpha(args::Vector{String})
    print_header()
    println("Mode: Red Giant Core Helium Burning & Alpha-Capture Ladder")

    # Helium flash / core helium burning conditions
    rho = 1.0e4       # g / cm^3
    T9  = 0.20        # 200 million K
    t_end = 3.1536e12 # ~ 100,000 years

    net = triple_alpha_network()
    println(@sprintf("Network: %d species, %d reaction channels", length(net), length(net.reactions)))
    println(@sprintf("Conditions: Density ρ = %.1e g/cm³, Temperature T = %.1f × 10⁸ K (T9 = %.2f)", rho, T9 * 10, T9))

    # Initial composition: Pure Helium-4
    Y0 = zeros(Float64, length(net))
    Y0[net.symbol_to_index["he4"]] = 1.0 / 4.0

    integ = NetworkIntegrator(rtol=1e-5, atol=1e-12)
    println("\nIntegrating 3α and sequential α-capture ladder...")
    res = solve_burn(net, Y0, t_end; rho=rho, T9=T9, integrator=integ, dt_init=1.0)

    println(@sprintf("Completed in %d accepted steps (%d rejected steps)", res.accepted_steps, res.rejected_steps))
    println("\nSynthesized Nucleosynthesis Products:")
    println("--------------------------------------------------------------------------------")
    println(" Isotope | Mass Number A | Molar Abundance Y_i | Mass Fraction X_i")
    println("--------------------------------------------------------------------------------")
    final_Y = res.abundances[:, end]
    for (i, iso) in enumerate(net.species_net.species)
        X_i = iso.A * final_Y[i]
        println(@sprintf(" %-7s | %13d | %19.6e | %17.6e", iso.symbol, iso.A, final_Y[i], X_i))
    end
    println("--------------------------------------------------------------------------------")

    eps_nuc = compute_energy_generation(net, final_Y, rho, T9)
    println(@sprintf("Helium Burning Power (ε_nuc):       %12.4e erg g⁻¹ s⁻¹", eps_nuc))
end

function cmd_nse(args::Vector{String})
    print_header()
    println("Mode: Nuclear Statistical Equilibrium (NSE) Supernova Core")

    # Core-collapse supernova / Silicon burning conditions
    rho = 1.0e7       # g / cm^3
    T9  = 5.0         # 5 Billion Kelvin
    Ye  = 0.50        # Symmetric nuclear matter

    net = NuclearNetwork(["n", "p", "he4", "c12", "o16", "ne20", "mg24", "si28", "ni56"])
    println(@sprintf("Species: %d isotopes evaluated in thermodynamic equilibrium", length(net)))
    println(@sprintf("State: Density ρ = %.1e g/cm³, T9 = %.2f (5 GK), Electron Fraction Ye = %.3f", rho, T9, Ye))

    println("\nSolving 2D Newton-Raphson Saha-NSE system...")
    sol = solve_nse(net, rho, T9, Ye)

    println(@sprintf("Convergence: %s after %d iterations", sol.converged ? "SUCCESS" : "FAILED", sol.iterations))
    println(@sprintf("Free Proton Abundance (Y_p):  %.6e", sol.Y_p))
    println(@sprintf("Free Neutron Abundance (Y_n): %.6e", sol.Y_n))
    println("\nEquilibrium Abundances (Iron-Peak / Nickel-56 Dominance):")
    println("--------------------------------------------------------------------------------")
    println(" Isotope | Z | N | Mass A | Binding (MeV) | Mass Fraction X_i")
    println("--------------------------------------------------------------------------------")
    for (i, iso) in enumerate(net.species)
        b_mev = binding_energy_mev(iso)
        println(@sprintf(" %-7s | %2d| %2d| %6d | %13.3f | %17.6e", iso.symbol, iso.Z, iso.N, iso.A, b_mev, sol.X[i]))
    end
    println("--------------------------------------------------------------------------------")
    println(@sprintf("Total Mass Fraction Conservation:   %.10f", sum(sol.X)))
end

function cmd_benchmark(args::Vector{String})
    print_header()
    println("Mode: Stiff Solver Performance & Throughput Benchmark")

    net = full_stellar_network()
    n = length(net)
    println(@sprintf("Full stellar network: %d species, %d reaction channels", n, length(net.reactions)))

    Y = fill(1.0 / n, n)
    rho = 100.0
    T9 = 0.05

    dY = zeros(Float64, n)
    J = zeros(Float64, n, n)
    Y_next = zeros(Float64, n)

    # 1. Warm-up
    compute_rhs!(dY, net, Y, rho, T9)
    compute_jacobian!(J, net, Y, rho, T9)
    step_linearly_implicit!(Y_next, net, Y, rho, T9, 1.0)

    # 2. Benchmark RHS evaluations
    n_rhs = 100_000
    t0 = time_ns()
    for _ in 1:n_rhs
        compute_rhs!(dY, net, Y, rho, T9)
    end
    t_rhs = (time_ns() - t0) / 1e9
    rate_rhs = n_rhs / t_rhs

    # 3. Benchmark Jacobian evaluations
    n_jac = 50_000
    t0 = time_ns()
    for _ in 1:n_jac
        compute_jacobian!(J, net, Y, rho, T9)
    end
    t_jac = (time_ns() - t0) / 1e9
    rate_jac = n_jac / t_jac

    # 4. Benchmark Full Linearly Implicit Steps
    n_steps = 25_000
    t0 = time_ns()
    for _ in 1:n_steps
        step_linearly_implicit!(Y_next, net, Y, rho, T9, 1.0)
    end
    t_steps = (time_ns() - t0) / 1e9
    rate_steps = n_steps / t_steps

    println("\nBenchmark Results:")
    println(@sprintf("  RHS Evaluations:        %10d in %6.3f s (%.2e evals/sec)", n_rhs, t_rhs, rate_rhs))
    println(@sprintf("  Jacobian Evaluations:   %10d in %6.3f s (%.2e evals/sec)", n_jac, t_jac, rate_jac))
    println(@sprintf("  Linearly-Implicit Steps: %9d in %6.3f s (%.2e steps/sec)", n_steps, t_steps, rate_steps))
    println("All kernels executed with optimal analytical linear algebra.")
end

function run_cli(args::Vector{String}=ARGS)
    if isempty(args) || args[1] in ["-h", "--help", "help"]
        print_header()
        println("Usage: julia --project=. bin/pyronucleus <command> [options]\n")
        println("Commands:")
        println("  solar-core     Simulate solar core hydrogen burning via the PP-chain")
        println("  cno-burn       Simulate massive star hydrogen burning via the CNO cycle")
        println("  triple-alpha   Simulate red giant helium burning and the alpha ladder")
        println("  nse            Calculate Nuclear Statistical Equilibrium (NSE) composition")
        println("  benchmark      Run stiff solver and Jacobian performance benchmarks")
        return
    end

    cmd = args[1]
    cmd_args = args[2:end]

    if cmd == "solar-core"
        cmd_solar_core(cmd_args)
    elseif cmd == "cno-burn"
        cmd_cno_burn(cmd_args)
    elseif cmd == "triple-alpha"
        cmd_triple_alpha(cmd_args)
    elseif cmd == "nse"
        cmd_nse(cmd_args)
    elseif cmd == "benchmark"
        cmd_benchmark(cmd_args)
    else
        println("Unknown command: $cmd")
        println("Run with --help to see available commands.")
    end
end

end # module CLI
