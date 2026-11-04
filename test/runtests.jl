using Test
using LinearAlgebra

# Ensure local package is in LOAD_PATH if running as standalone script
pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using Pyronucleus

@testset "Pyronucleus Test Suite" begin

    @testset "1. Species & Isotopic Database" begin
        @test haskey(ISOTOPE_DATABASE, "p")
        @test haskey(ISOTOPE_DATABASE, "he4")
        @test haskey(ISOTOPE_DATABASE, "c12")
        @test haskey(ISOTOPE_DATABASE, "ni56")

        p = ISOTOPE_DATABASE["p"]
        @test p.Z == 1 && p.A == 1 && p.N == 0
        @test isapprox(p.mass_excess_mev, 7.288971, atol=1e-5)

        c12 = ISOTOPE_DATABASE["c12"]
        @test c12.Z == 6 && c12.A == 12 && c12.N == 6
        @test c12.mass_excess_mev == 0.0 # Definition of atomic mass scale

        net = NuclearNetwork(["p", "he4", "c12"])
        @test length(net) == 3
        @test net.symbol_to_index["p"] == 1
        @test net.symbol_to_index["he4"] == 2
        @test net.symbol_to_index["c12"] == 3

        # Mass fractions: X(H) = 0.75, X(He4) = 0.25 => Y(H) = 0.75, Y(He4) = 0.0625
        Y = [0.75, 0.0625, 0.0]
        X = mass_fractions(net, Y)
        @test isapprox(X[1], 0.75, atol=1e-7)
        @test isapprox(X[2], 0.25, atol=1e-7)
        @test isapprox(sum(X), 1.0, atol=1e-7)
        @test validate_conservation(net, Y)

        # Electron fraction: Ye = sum(Z_i * Y_i) = 1*0.75 + 2*0.0625 = 0.875
        @test isapprox(electron_fraction(net, Y), 0.875, atol=1e-7)

        # Mean molecular weights
        mu_ion, mu_e, mu_tot = mean_molecular_weights(net, Y)
        @test isapprox(mu_ion, 1.0 / (0.75 + 0.0625), atol=1e-5)
        @test isapprox(mu_e, 1.0 / 0.875, atol=1e-5)
        @test isapprox(mu_tot, 1.0 / (0.8125 + 0.875), atol=1e-5)
    end

    @testset "2. Reaction Rates & Analytical Derivatives" begin
        @test haskey(REACTION_DATABASE, "pp")
        @test haskey(REACTION_DATABASE, "triple_alpha")
        @test haskey(REACTION_DATABASE, "c12alpha_o16")

        rxn_pp = REACTION_DATABASE["pp"]
        @test rxn_pp.type == Pyronucleus.Rates.BINARY
        @test rxn_pp.reactants == ["p", "p"]
        @test rxn_pp.products == ["d"]

        # Rate should be positive and monotonically increasing with temperature in the 10^7 K regime
        r1 = evaluate_rate(rxn_pp, 0.01) # 10 MK
        r2 = evaluate_rate(rxn_pp, 0.02) # 20 MK
        @test r1 > 0.0
        @test r2 > r1

        # Test analytical temperature derivative vs centered finite difference
        for T9_test in [0.015, 0.05, 0.2, 1.0]
            for (key, rxn) in [("pp", rxn_pp), ("c12p_n13", REACTION_DATABASE["c12p_n13"])]
                dlam_analytical = rate_temperature_derivative(rxn, T9_test)
                eps_T = 1e-6 * T9_test
                lam_plus = evaluate_rate(rxn, T9_test + eps_T)
                lam_minus = evaluate_rate(rxn, T9_test - eps_T)
                dlam_numerical = (lam_plus - lam_minus) / (2.0 * eps_T)
                
                @test isapprox(dlam_analytical, dlam_numerical, rtol=1e-4)
            end
        end

        # Test Salpeter screening factor
        net = pp_chain_network()
        Y_sun = [0.73, 1e-10, 1e-5, 0.0625, 0.0, 0.0, 0.0]
        f_scr = screening_factor(rxn_pp, 150.0, 0.0157, Y_sun, net.species_net)
        @test f_scr >= 1.0 # Screening always enhances forward rate
        @test f_scr < 2.0  # Weak screening in solar core is typically ~ 1.05 - 1.20
    end

    @testset "3. Network Topologies & Stoichiometry" begin
        pp_net = pp_chain_network()
        @test length(pp_net) == 7
        @test length(pp_net.reactions) == 6

        cno_net = cno_cycle_network()
        @test length(cno_net) == 9
        @test length(cno_net.reactions) == 7

        tri_net = triple_alpha_network()
        @test length(tri_net) == 6
        @test length(tri_net.reactions) == 5

        full_net = full_stellar_network()
        @test length(full_net) == 17
        @test length(full_net.reactions) == 18

        # Check total baryon conservation of stoichiometry for each reaction:
        # sum_i (A_i * S_{i,k}) == 0
        for (k, rxn) in enumerate(full_net.reactions)
            mass_balance = sum(full_net.species_net.A_vector .* full_net.stoichiometry[:, k])
            @test isapprox(mass_balance, 0.0, atol=1e-10)
        end
    end

    @testset "4. Exact Analytical Jacobian vs Finite Differences" begin
        net = full_stellar_network()
        n = length(net)
        rho = 100.0
        T9 = 0.03 # 30 million Kelvin
        Y = fill(1.0 / n, n) # arbitrary composition

        J_analytical = compute_jacobian(net, Y, rho, T9)
        J_numerical = zeros(Float64, n, n)

        eps_y = 1e-7
        for j in 1:n
            Y_plus = copy(Y)
            Y_minus = copy(Y)
            Y_plus[j] += eps_y
            Y_minus[j] -= eps_y

            f_plus = compute_rhs(net, Y_plus, rho, T9)
            f_minus = compute_rhs(net, Y_minus, rho, T9)

            J_numerical[:, j] = (f_plus - f_minus) / (2.0 * eps_y)
        end

        for i in 1:n, j in 1:n
            @test isapprox(J_analytical[i, j], J_numerical[i, j], atol=1e-8, rtol=1e-3)
        end
    end

    @testset "5. Energetics & Neutrino Losses" begin
        net = pp_chain_network()
        Y = zeros(Float64, length(net))
        Y[net.symbol_to_index["p"]] = 0.75
        Y[net.symbol_to_index["he4"]] = 0.25 / 4.0
        rho = 150.0
        T9 = 0.0157

        eps_nuc = compute_energy_generation(net, Y, rho, T9)
        eps_nu  = compute_neutrino_losses(net, Y, rho, T9)
        eps_th  = net_thermal_energy_rate(net, Y, rho, T9)

        @test eps_nuc > 0.0
        @test eps_nu > 0.0
        @test eps_th > 0.0
        @test isapprox(eps_th, eps_nuc - eps_nu, atol=1e-12)
        @test eps_th < eps_nuc

        # Verify equivalence between sum(r_k * Q_k) and -sum(dY_i/dt * Δ_i)
        fluxes = Pyronucleus.Networks.compute_fluxes(net, Y, rho, T9)
        dY = compute_rhs(net, Y, rho, T9)
        q_from_fluxes = sum(fluxes[k] * net.reactions[k].q_value_mev for k in 1:length(net.reactions))
        q_from_dY = -sum(dY[i] * net.species_net.mass_excess[i] for i in 1:length(net))
        @test isapprox(q_from_fluxes, q_from_dY, rtol=1e-8)
    end

    @testset "6. Stiff ODE Burn Integrator" begin
        net = pp_chain_network()
        Y0 = zeros(Float64, length(net))
        Y0[net.symbol_to_index["p"]] = 0.70
        Y0[net.symbol_to_index["he4"]] = 0.30 / 4.0

        # Burn for 100,000 years in solar core
        t_end = 3.1536e12
        rho = 150.0
        T9 = 0.0157

        res = solve_burn(net, Y0, t_end; rho=rho, T9=T9, dt_init=1.0)
        @test res.accepted_steps > 0
        @test res.times[end] == t_end

        # Test mass conservation across all integration steps
        for step in 1:size(res.abundances, 2)
            Y_step = res.abundances[:, step]
            @test validate_conservation(net.species_net, Y_step; tol=1e-8)
        end

        # Verify hydrogen was depleted and helium/deuterium produced
        p_idx = net.symbol_to_index["p"]
        d_idx = net.symbol_to_index["d"]
        @test res.abundances[p_idx, end] < res.abundances[p_idx, 1]
        @test res.abundances[d_idx, end] > 0.0
    end

    @testset "7. Nuclear Statistical Equilibrium (NSE)" begin
        net = NuclearNetwork(["n", "p", "he4", "c12", "o16", "ne20", "mg24", "si28", "ni56"])
        rho = 1e7
        T9 = 5.0
        Ye = 0.50

        sol = solve_nse(net, rho, T9, Ye)
        @test sol.converged
        @test sol.iterations <= 25

        # Test mass conservation: sum(A_i * Y_i) == 1.0
        @test isapprox(sum(sol.X), 1.0, atol=1e-8)

        # Test electron fraction constraint: sum(Z_i * Y_i) == Ye
        @test isapprox(electron_fraction(net, sol.Y), Ye, atol=1e-5)

        # At T9 = 5.0 and Ye = 0.5, 28Si has the highest binding free energy per nucleon among alpha species and dominates
        si28_idx = net.symbol_to_index["si28"]
        @test sol.X[si28_idx] > 0.50
    end

    @testset "8. CLI Execution Integrity" begin
        # Test help screen
        @test_nowarn Pyronucleus.run_cli(["--help"])
        # Test subcommands
        @test_nowarn Pyronucleus.run_cli(["solar-core"])
        @test_nowarn Pyronucleus.run_cli(["cno-burn"])
        @test_nowarn Pyronucleus.run_cli(["triple-alpha"])
        @test_nowarn Pyronucleus.run_cli(["nse"])
        @test_nowarn Pyronucleus.run_cli(["benchmark"])
    end

end
