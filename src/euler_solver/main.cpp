#include "io.hpp"

// Main solver
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
    try {
        Config c;
        parse_args(c, argc, argv);
        validate_config(c);

        const double dx = (c.xhi - c.xlo) / c.nx;
        const double dy = (c.yhi - c.ylo) / c.ny;
        const int N = c.nx * c.ny;

        const Prim Winf = make_freestream(c);
        const Cons Uinf = prim_to_cons(Winf, c.gamma);
        const Cons Uscale = conservative_reference_scale(c);

        const auto geom = build_geometry(c, dx, dy);
        const auto panels = build_panels(c);

        FaceData base;
        base.build(c, geom, panels, dx, dy);

        std::vector<Cons> U(N, Uinf), U1(N, Uinf), U2(N, Uinf), RHS(N);
        std::vector<Prim> W(N);
        std::vector<int> tlag(N, 0);

        std::vector<bool> tx_stage((c.nx + 1) * c.ny);
        std::vector<bool> ty_stage(c.nx * (c.ny + 1));
        std::vector<bool> tx_all((c.nx + 1) * c.ny);
        std::vector<bool> ty_all(c.nx * (c.ny + 1));

        std::vector<double> r2_hist, ri_hist, nr_hist, nre_hist, du_hist, due_hist;
        std::vector<ForceHistoryEntry> fhist;
        std::vector<DiagnosticsHistoryEntry> dhist;

        r2_hist.reserve(c.max_iters);
        ri_hist.reserve(c.max_iters);
        nr_hist.reserve(c.max_iters);
        nre_hist.reserve(c.max_iters);
        du_hist.reserve(c.max_iters);
        due_hist.reserve(c.max_iters);

        double initial_r2 = -1.0;
        double nre = -1.0;
        double due = -1.0;
        const double ema = 0.08;

        const double beta_ack = std::sqrt(c.M_inf * c.M_inf - 1.0);
        const double theta_ack = c.t_over_c;
        const double alpha = aoa_rad(c);
        const double Cd_ack = 4.0 * (theta_ack * theta_ack + alpha * alpha) / beta_ack;
        const double Cl_ack = 4.0 * alpha / beta_ack;
        std::cout << "Supersonic diamond airfoil Euler solver\n";
        std::cout << "Grid: " << c.nx << " x " << c.ny << "\n";
        std::cout << "Mach: " << c.M_inf << "  AoA: " << c.aoa_deg << " deg\n";
        std::cout << "CFL: " << c.cfl << "  Max iterations: " << c.max_iters << "\n";
        std::cout << "Output: " << c.output_dir << "\n";
        std::cout << "Ackeret: Cd=" << Cd_ack << " Cl=" << Cl_ack << "\n";
        std::cout << "Prandtl-Meyer angle: "
                  << prandtl_meyer_angle(c.M_inf, c.gamma) * 180.0 / PI << " deg\n\n";

        for (int iter = 0; iter < c.max_iters; ++iter) {
            FaceData fd = base;
            fd.apply_trouble_lag(tlag);

            std::fill(tx_all.begin(), tx_all.end(), false);
            std::fill(ty_all.begin(), ty_all.end(), false);

            update_primitive_cache(W, U, geom, c, Winf);
            const double dt = compute_dt(W, geom, c, dx, dy);
            Diagnostics diag;

            const bool do_diag = ((iter + 1) % c.diagnostics_interval == 0 || iter == 0);

            const ResidualInfo res0 = compute_rhs(
                W, geom, Winf, fd, RHS, tx_stage, ty_stage, c, dx, dy, do_diag ? &diag : nullptr
            );

            for (std::size_t q = 0; q < tx_all.size(); ++q) tx_all[q] = tx_all[q] || tx_stage[q];
            for (std::size_t q = 0; q < ty_all.size(); ++q) ty_all[q] = ty_all[q] || ty_stage[q];

            for (int k = 0; k < N; ++k) {
                if (geom[k].solid) {
                    U1[k] = U[k];
                    continue;
                }
                U1[k] = U[k] + RHS[k] * dt;
                if (!physical(cons_to_prim(U1[k], c.gamma), c)) {
                    U1[k] = U[k];
                }
            }

            update_primitive_cache(W, U1, geom, c, Winf);
            compute_rhs(W, geom, Winf, fd, RHS, tx_stage, ty_stage, c, dx, dy, nullptr);

            for (std::size_t q = 0; q < tx_all.size(); ++q) tx_all[q] = tx_all[q] || tx_stage[q];
            for (std::size_t q = 0; q < ty_all.size(); ++q) ty_all[q] = ty_all[q] || ty_stage[q];

            for (int k = 0; k < N; ++k) {
                if (geom[k].solid) {
                    U2[k] = U[k];
                    continue;
                }
                U2[k] = U[k] * 0.75 + (U1[k] + RHS[k] * dt) * 0.25;
                if (!physical(cons_to_prim(U2[k], c.gamma), c)) {
                    U2[k] = U[k];
                }
            }

            update_primitive_cache(W, U2, geom, c, Winf);
            compute_rhs(W, geom, Winf, fd, RHS, tx_stage, ty_stage, c, dx, dy, nullptr);

            for (std::size_t q = 0; q < tx_all.size(); ++q) tx_all[q] = tx_all[q] || tx_stage[q];
            for (std::size_t q = 0; q < ty_all.size(); ++q) ty_all[q] = ty_all[q] || ty_stage[q];

            double du = 0.0;

            for (int k = 0; k < N; ++k) {
                if (geom[k].solid) {
                    continue;
                }

                Cons Un = U[k] * (1.0 / 3.0) + (U2[k] + RHS[k] * dt) * (2.0 / 3.0);

                if (!physical(cons_to_prim(Un, c.gamma), c)) {
                    Un = U[k];
                }

                du = std::max({
                    du,
                    std::abs(Un.r - U[k].r) / Uscale.r,
                    std::abs(Un.ru - U[k].ru) / Uscale.ru,
                    std::abs(Un.rv - U[k].rv) / Uscale.rv,
                    std::abs(Un.rE - U[k].rE) / Uscale.rE
                });

                U[k] = Un;
            }

            apply_sponge(U, geom, Uinf, c, dx, dy, dt);

            for (auto& v : tlag) {
                if (v > 0) --v;
            }

            for (int j = 0; j < c.ny; ++j) {
                for (int iface = 1; iface < c.nx; ++iface) {
                    const int fi = iface * c.ny + j;
                    if (tx_all[fi]) {
                        tlag[idx(iface - 1, j, c.nx)] = c.trouble_lag;
                        tlag[idx(iface, j, c.nx)] = c.trouble_lag;
                    }
                }
            }

            for (int i = 0; i < c.nx; ++i) {
                for (int jf = 1; jf < c.ny; ++jf) {
                    const int fi = i * (c.ny + 1) + jf;
                    if (ty_all[fi]) {
                        tlag[idx(i, jf - 1, c.nx)] = c.trouble_lag;
                        tlag[idx(i, jf, c.nx)] = c.trouble_lag;
                    }
                }
            }

            int ntr = 0;
            for (auto v : tlag) {
                if (v > 0) ++ntr;
            }

            if (do_diag) {
                diag.n_troubled_cells = ntr;
                dhist.push_back({ iter + 1, dt, diag });
            }

            if (initial_r2 < 0.0) {
                initial_r2 = std::max(res0.l2, 1.0e-30);
            }

            const double nr = res0.l2 / initial_r2;
            nre = (nre < 0.0) ? nr : ema * nr + (1.0 - ema) * nre;
            due = (due < 0.0) ? du : ema * du + (1.0 - ema) * due;

            r2_hist.push_back(res0.l2);
            ri_hist.push_back(res0.linf);
            nr_hist.push_back(nr);
            nre_hist.push_back(nre);
            du_hist.push_back(du);
            due_hist.push_back(due);

            if ((iter + 1) % c.force_interval == 0 || iter == 0) {
                const Forces F = compute_forces(U, geom, panels, c, dx, dy);
                fhist.push_back({ iter + 1, res0.l2, res0.linf, nr, nre, du, due, F.Cd, F.Cl, F.Cm_LE, F.Cm_c4, F.L_over_D });
            }

            write_snapshot_if_needed(iter, U, geom, c);

            if ((iter + 1) % 100 == 0 || iter == 0) {
                std::cout << "Iter " << std::setw(7) << (iter + 1)
                          << " nres_l2=" << std::scientific << std::setprecision(3) << nr
                          << " nres_ema=" << nre
                          << " rhs_inf=" << res0.linf
                          << " du=" << du
                          << " du_ema=" << due << "\n";
                std::cout.unsetf(std::ios::floatfield);
            }

            if ((iter + 1) >= c.min_iters && force_history_stable(fhist, c) && due < c.solution_change_tol) {
                std::cout << "\nStopping criteria met at iteration " << (iter + 1) << "\n";
                break;
            }
        }

        std::vector<PanelForce> pfs_panel;
        std::vector<PanelForce> pfs_wall;

        const Forces F_panel = compute_forces_panel_sampling(U, geom, panels, c, dx, dy, &pfs_panel);
        const Forces F_wall = compute_forces_wall_pressure(U, geom, panels, c, dx, dy, &pfs_wall);
        const Forces F = c.use_panel_force_for_primary_coefficients ? F_panel : F_wall;
        const std::vector<PanelForce>& pfs_primary =
            c.use_panel_force_for_primary_coefficients ? pfs_panel : pfs_wall;

        auto pct_diff = [](double ref, double val) {
            return 100.0 * (val - ref) / std::max(std::abs(ref), 1.0e-30);
        };

        std::cout << "\nFinal forces\n"
                  << std::fixed << std::setprecision(6)
                  << " Method = " << (c.use_panel_force_for_primary_coefficients ? "PANEL_SAMPLING" : "WALL_PRESSURE") << "\n"
                  << " Fx = " << F.Fx << " N/m\n"
                  << " Fy = " << F.Fy << " N/m\n"
                  << " Drag = " << F.Drag << " N/m\n"
                  << " Lift = " << F.Lift << " N/m\n"
                  << " Cd = " << F.Cd << "\n"
                  << " Cl = " << F.Cl << "\n"
                  << " Cm_LE = " << F.Cm_LE << "\n"
                  << " Cm_c4 = " << F.Cm_c4 << "\n"
                  << " L/D = " << F.L_over_D << "\n"
                  << " Ackeret: Cd=" << Cd_ack << " Cl=" << Cl_ack << "\n"
                  << " Error: Cd=" << std::showpos << 100.0 * (F.Cd - Cd_ack) / Cd_ack
                  << "% Cl=" << 100.0 * (F.Cl - Cl_ack) / Cl_ack << "%\n"
                  << std::noshowpos;

        std::cout << "\nForce method comparison\n"
                  << " Panel Cd = " << F_panel.Cd << "  Cl = " << F_panel.Cl << "\n"
                  << " Wall   Cd = " << F_wall.Cd << "  Cl = " << F_wall.Cl << "\n"
                  << " Difference wall vs panel: Cd = " << pct_diff(F_panel.Cd, F_wall.Cd)
                  << "%, Cl = " << pct_diff(F_panel.Cl, F_wall.Cl) << "%\n";

        if (std::abs(pct_diff(F_panel.Cd, F_wall.Cd)) > 2.0
            || std::abs(pct_diff(F_panel.Cl, F_wall.Cl)) > 2.0) {
            std::cout << "\nWarning: force methods differ by more than 2%.\n"
                      << "Check the pressure sampling before using the wall-pressure force result.\n";
        }

        write_solution_csv(make_output_path(c, c.solution_csv), U, geom, c);
        write_residual_csv(make_output_path(c, c.residual_csv), r2_hist, ri_hist, nr_hist, nre_hist, du_hist, due_hist);
        write_force_csv(make_output_path(c, c.force_csv), F);
        write_force_csv(make_output_path(c, c.force_panel_csv), F_panel);
        write_force_csv(make_output_path(c, c.force_wall_pressure_csv), F_wall);
        write_force_comparison_csv(make_output_path(c, c.force_comparison_csv), F_panel, F_wall);
        write_force_history_csv(make_output_path(c, c.force_history_csv), fhist);
        write_panel_forces_csv(make_output_path(c, c.panel_forces_csv), pfs_primary);
        write_surface_cp_csv(make_output_path(c, c.surface_cp_csv), U, geom, panels, c, dx, dy);
        write_wall_distance_csv(make_output_path(c, c.wall_distance_csv), geom, panels, c);
        write_diagnostics_csv(make_output_path(c, c.diagnostics_csv), dhist);

        std::cout << "\nWrote outputs to: " << c.output_dir << "\n";

        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
