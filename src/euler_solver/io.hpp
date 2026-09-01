#pragma once

#include "forces.hpp"

// Output files
// -----------------------------------------------------------------------------

static void write_solution_csv(const std::string& path, const std::vector<Cons>& U, const std::vector<CellGeom>& geom, const Config& c) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    const Prim Winf = make_freestream(c);
    const double q = 0.5 * Winf.rho * velocity_magnitude(Winf) * velocity_magnitude(Winf);

    out << "x,y,is_solid,rho,u,v,p,Cp,T,a,Mach,entropy_like\n";

    for (int j = 0; j < c.ny; ++j) {
        for (int i = 0; i < c.nx; ++i) {
            const int k = idx(i, j, c.nx);

            if (geom[k].solid) {
                out << geom[k].x << "," << geom[k].y << ",1,nan,nan,nan,nan,nan,nan,nan,nan,nan\n";
                continue;
            }

            const Prim W = cons_to_prim(U[k], c.gamma);
            const double T = W.p / (W.rho * c.gas_R);
            const double a = sound_speed(W, c.gamma);
            const double M = velocity_magnitude(W) / a;
            const double Cp = (W.p - c.p_inf) / q;
            const double S = W.p / std::pow(W.rho, c.gamma);

            out << std::setprecision(12)
                << geom[k].x << "," << geom[k].y << ",0,"
                << W.rho << "," << W.u << "," << W.v << "," << W.p << ","
                << Cp << "," << T << "," << a << "," << M << "," << S << "\n";
        }
    }
}

static void write_residual_csv(
    const std::string& path,
    const std::vector<double>& r2,
    const std::vector<double>& ri,
    const std::vector<double>& nr,
    const std::vector<double>& nre,
    const std::vector<double>& du,
    const std::vector<double>& due
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "iteration,residual_l2,residual_linf,normalized_residual_l2,normalized_residual_l2_ema,solution_change,solution_change_ema\n";

    for (std::size_t i = 0; i < r2.size(); ++i) {
        out << (i + 1) << "," << std::setprecision(12)
            << r2[i] << "," << ri[i] << "," << nr[i] << "," << nre[i] << "," << du[i] << "," << due[i] << "\n";
    }
}

static void write_force_csv(const std::string& path, const Forces& F) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "Fx,Fy,Drag,Lift,Cd,Cl,Mz_LE,Mz_c4,Cm_LE,Cm_c4,L_over_D\n"
        << std::setprecision(12)
        << F.Fx << "," << F.Fy << "," << F.Drag << "," << F.Lift << ","
        << F.Cd << "," << F.Cl << "," << F.Mz_LE << "," << F.Mz_c4 << ","
        << F.Cm_LE << "," << F.Cm_c4 << "," << F.L_over_D << "\n";
}

static void write_force_comparison_csv(
    const std::string& path,
    const Forces& panel,
    const Forces& wall
) {
    std::ofstream out(path);

    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    auto pct = [](double ref, double val) {
        return 100.0 * (val - ref) / std::max(std::abs(ref), 1.0e-30);
    };

    out << "quantity,panel_sampling,wall_pressure,percent_difference_wall_vs_panel\n";
    out << std::setprecision(12)
        << "Cd," << panel.Cd << "," << wall.Cd << "," << pct(panel.Cd, wall.Cd) << "\n"
        << "Cl," << panel.Cl << "," << wall.Cl << "," << pct(panel.Cl, wall.Cl) << "\n"
        << "Cm_LE," << panel.Cm_LE << "," << wall.Cm_LE << "," << pct(panel.Cm_LE, wall.Cm_LE) << "\n"
        << "Cm_c4," << panel.Cm_c4 << "," << wall.Cm_c4 << "," << pct(panel.Cm_c4, wall.Cm_c4) << "\n"
        << "L_over_D," << panel.L_over_D << "," << wall.L_over_D << "," << pct(panel.L_over_D, wall.L_over_D) << "\n";
}

static void write_force_history_csv(const std::string& path, const std::vector<ForceHistoryEntry>& h) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "iteration,residual_l2,residual_linf,normalized_residual_l2,normalized_residual_l2_ema,solution_change,solution_change_ema,Cd,Cl,Cm_LE,Cm_c4,L_over_D\n";

    for (const auto& e : h) {
        out << e.iteration << "," << std::setprecision(12)
            << e.residual_l2 << "," << e.residual_linf << ","
            << e.normalized_residual_l2 << "," << e.normalized_residual_l2_ema << ","
            << e.solution_change << "," << e.solution_change_ema << ","
            << e.Cd << "," << e.Cl << "," << e.Cm_LE << "," << e.Cm_c4 << "," << e.L_over_D << "\n";
    }
}

static void write_panel_forces_csv(const std::string& path, const std::vector<PanelForce>& pf) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "panel_id,surface,x_mid,y_mid,nx,ny,length,pressure_avg,Cp_avg,Fx,Fy,Drag,Lift,Cd,Cl,Mz_LE,Mz_c4,Cm_LE,Cm_c4,avg_sample_distance,avg_valid_samples,fallback_count\n";

    for (const auto& p : pf) {
        out << p.panel_id << "," << p.surface << "," << std::setprecision(12)
            << p.x_mid << "," << p.y_mid << "," << p.nx << "," << p.ny << "," << p.length << ","
            << p.pressure_avg << "," << p.Cp_avg << "," << p.Fx << "," << p.Fy << ","
            << p.Drag << "," << p.Lift << "," << p.Cd << "," << p.Cl << ","
            << p.Mz_LE << "," << p.Mz_c4 << "," << p.Cm_LE << "," << p.Cm_c4 << ","
            << p.avg_sample_distance << "," << p.avg_valid_samples << "," << p.fallback_count << "\n";
    }
}

static void write_surface_cp_csv(
    const std::string& path,
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const std::vector<Panel>& panels,
    const Config& c,
    double dx,
    double dy
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "panel_id,surface,s_over_panel,x,y,nx,ny,pressure,Cp,sample_distance,valid_samples,used_fallback\n";

    const int ns = std::max(2, c.surface_samples_per_panel);

    for (std::size_t pid = 0; pid < panels.size(); ++pid) {
        const Panel& pan = panels[pid];
        const Vec2 ab = sub(pan.b, pan.a);

        for (int s = 0; s < ns; ++s) {
            const double eta = static_cast<double>(s) / static_cast<double>(ns - 1);
            const Vec2 pnt = { pan.a.x + eta * ab.x, pan.a.y + eta * ab.y };

            const PressureSample ps = sample_wall_pressure(pnt, pan.n, U, geom, c, dx, dy);

            out << pid << "," << panel_surface_name(static_cast<int>(pid)) << "," << std::setprecision(12)
                << eta << "," << pnt.x << "," << pnt.y << "," << pan.n.x << "," << pan.n.y << ","
                << ps.pressure << "," << ps.Cp << "," << ps.sample_distance << ","
                << ps.valid_samples << "," << ps.used_fallback << "\n";
        }
    }
}

static void write_wall_distance_csv(const std::string& path, const std::vector<CellGeom>& geom, const std::vector<Panel>& panels, const Config& c) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "x,y,is_solid,wall_distance,wall_normal_x,wall_normal_y\n";

    for (int j = 0; j < c.ny; ++j) {
        for (int i = 0; i < c.nx; ++i) {
            const int k = idx(i, j, c.nx);
            const Vec2 p{ geom[k].x, geom[k].y };
            const Vec2 n = wall_normal(p, panels);
            out << std::setprecision(12)
                << p.x << "," << p.y << "," << (geom[k].solid ? 1 : 0) << ","
                << wall_distance(p, panels) << "," << n.x << "," << n.y << "\n";
        }
    }
}

static void write_diagnostics_csv(const std::string& path, const std::vector<DiagnosticsHistoryEntry>& h) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path);
    }

    out << "iteration,dt,n_faces_x,n_faces_y,n_hllc_x,n_hllc_y,n_llf_x,n_llf_y,n_weno_x,n_weno_y,n_first_order_x,n_first_order_y,n_wall_x,n_wall_y,n_troubled_x,n_troubled_y,n_strong_shock_x,n_strong_shock_y,n_bad_reconstruction_x,n_bad_reconstruction_y,n_troubled_cells\n";

    for (const auto& e : h) {
        const auto d = e.d;
        out << e.iteration << "," << std::setprecision(12) << e.dt << ","
            << d.n_faces_x << "," << d.n_faces_y << ","
            << d.n_hllc_x << "," << d.n_hllc_y << ","
            << d.n_llf_x << "," << d.n_llf_y << ","
            << d.n_weno_x << "," << d.n_weno_y << ","
            << d.n_first_order_x << "," << d.n_first_order_y << ","
            << d.n_wall_x << "," << d.n_wall_y << ","
            << d.n_troubled_x << "," << d.n_troubled_y << ","
            << d.n_strong_shock_x << "," << d.n_strong_shock_y << ","
            << d.n_bad_reconstruction_x << "," << d.n_bad_reconstruction_y << ","
            << d.n_troubled_cells << "\n";
    }
}

// -----------------------------------------------------------------------------
// Convergence check
// -----------------------------------------------------------------------------

static bool force_history_stable(const std::vector<ForceHistoryEntry>& h, const Config& c) {
    const int n = static_cast<int>(h.size());
    if (n < c.force_window) {
        return false;
    }

    const int s = n - c.force_window;

    double cdmin = HUGE_VAL_NUM;
    double cdmax = -HUGE_VAL_NUM;
    double clmin = HUGE_VAL_NUM;
    double clmax = -HUGE_VAL_NUM;
    double cds = 0.0;
    double cls = 0.0;

    for (int i = s; i < n; ++i) {
        cdmin = std::min(cdmin, h[i].Cd);
        cdmax = std::max(cdmax, h[i].Cd);
        clmin = std::min(clmin, h[i].Cl);
        clmax = std::max(clmax, h[i].Cl);
        cds += h[i].Cd;
        cls += h[i].Cl;
    }

    const double cdm = cds / c.force_window;
    const double clm = cls / c.force_window;

    const double cd_band = (cdmax - cdmin) / std::max(std::abs(cdm), 1.0e-12);
    const double cl_band = (clmax - clmin) / std::max(std::abs(clm), 1.0e-12);

    return cd_band < c.force_tol && cl_band < c.force_tol;
}

static void write_snapshot_if_needed(int iter, const std::vector<Cons>& U, const std::vector<CellGeom>& geom, const Config& c) {
    if (!c.write_snapshots || (iter + 1) % c.snapshot_interval != 0) {
        return;
    }

    std::ostringstream os;
    os << c.snapshot_prefix << std::setw(6) << std::setfill('0') << (iter + 1) << ".csv";
    write_solution_csv(make_output_path(c, os.str()), U, geom, c);
}

// -----------------------------------------------------------------------------
// Command line
// -----------------------------------------------------------------------------

static void print_usage(const char* exe) {
    std::cout
        << "Usage:\n"
        << "  " << exe << " [options]\n\n"
        << "Main options:\n"
        << "  --mach <value>              Freestream Mach number, must be > 1\n"
        << "  --aoa <deg>                 Angle of attack in degrees\n"
        << "  --nx <cells> --ny <cells>   Grid resolution\n"
        << "  --max-iters <n>             Maximum pseudo-time iterations\n"
        << "  --min-iters <n>             Minimum iterations before convergence check\n"
        << "  --output <folder>           Output directory\n"
        << "  --primary-panel             Use panel sampling for forces\n"
        << "  --primary-wall              Use wall-pressure integration for forces\n"
        << "  --snapshots                 Enable solution snapshots\n"
        << "  --no-snapshots              Disable solution snapshots\n"
        << "  --help                      Print this message\n";
}

static void parse_args(Config& c, int argc, char** argv) {
    for (int a = 1; a < argc; ++a) {
        const std::string key = argv[a];

        auto val = [&](const std::string& name) -> std::string {
            if (a + 1 >= argc) {
                throw std::runtime_error("Missing value for argument " + name);
            }
            return std::string(argv[++a]);
        };

        if (key == "--help" || key == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        }
        else if (key == "--nx") c.nx = std::stoi(val(key));
        else if (key == "--ny") c.ny = std::stoi(val(key));
        else if (key == "--mach") c.M_inf = std::stod(val(key));
        else if (key == "--aoa") c.aoa_deg = std::stod(val(key));
        else if (key == "--max-iters") c.max_iters = std::stoi(val(key));
        else if (key == "--min-iters") c.min_iters = std::stoi(val(key));
        else if (key == "--cfl") c.cfl = std::stod(val(key));
        else if (key == "--output") c.output_dir = val(key);
        else if (key == "--snapshots") c.write_snapshots = true;
        else if (key == "--no-snapshots") c.write_snapshots = false;
        else if (key == "--trouble-ratio") c.trouble_ratio = std::stod(val(key));
        else if (key == "--trouble-lag") c.trouble_lag = std::stoi(val(key));
        else if (key == "--body-exclusion") c.body_exclusion_chord = std::stod(val(key));
        else if (key == "--wall-first-offset") c.wall_pressure_first_offset_cells = std::stod(val(key));
        else if (key == "--wall-step") c.wall_pressure_step_cells = std::stod(val(key));
        else if (key == "--wall-max-offset") c.wall_pressure_max_offset_cells = std::stod(val(key));
        else if (key == "--wall-rays") c.wall_pressure_normal_rays = std::stoi(val(key));
        else if (key == "--wall-min-samples") c.wall_pressure_min_valid_samples = std::stoi(val(key));
        else if (key == "--wall-max-samples") c.wall_pressure_max_valid_samples = std::stoi(val(key));
        else if (key == "--wall-tangent-radius") c.wall_pressure_tangent_radius = std::stoi(val(key));
        else if (key == "--panel-skip") c.panel_force_skip_cells = std::stoi(val(key));
        else if (key == "--panel-collect") c.panel_force_collect_cells = std::stoi(val(key));
        else if (key == "--primary-wall") c.use_panel_force_for_primary_coefficients = false;
        else if (key == "--primary-panel") c.use_panel_force_for_primary_coefficients = true;
        else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }
}

// -----------------------------------------------------------------------------
// Input checks
// -----------------------------------------------------------------------------

static void validate_config(const Config& c) {
    if (c.nx < 40 || c.ny < 40) {
        throw std::runtime_error("Grid is too small for this immersed-boundary Euler solver.");
    }

    if (!(c.M_inf > c.minimum_allowed_mach)) {
        throw std::runtime_error("This solver configuration is supersonic only. Use --mach greater than 1.");
    }

    if (!(c.cfl > 0.0 && c.cfl < 1.0)) {
        throw std::runtime_error("CFL must be positive and below 1 for this explicit RK solver.");
    }

    if (c.force_window < 2 || c.force_interval < 1) {
        throw std::runtime_error("Invalid force-history convergence settings.");
    }

    if (c.max_iters < 1 || c.min_iters < 0 || c.min_iters > c.max_iters) {
        throw std::runtime_error("Invalid iteration settings.");
    }

    if (!(c.gamma > 1.0) || !(c.gas_R > 0.0) || !(c.p_inf > 0.0) || !(c.T_inf > 0.0)) {
        throw std::runtime_error("Invalid gas or freestream settings.");
    }

    if (!(c.chord > 0.0) || !(c.t_over_c > 0.0)) {
        throw std::runtime_error("Invalid airfoil geometry.");
    }

    if (c.panel_force_collect_cells < 1 || c.surface_samples_per_panel < 2 || c.force_samples_per_panel < 1) {
        throw std::runtime_error("Invalid force-sampling settings.");
    }
}

// -----------------------------------------------------------------------------
