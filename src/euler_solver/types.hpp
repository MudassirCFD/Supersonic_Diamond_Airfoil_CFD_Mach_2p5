#pragma once

// Supersonic Diamond Airfoil
// Two-dimensional compressible Euler solver

// -----------------------------------------------------------------------------
// Includes
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Constants and data types
// -----------------------------------------------------------------------------

static constexpr double PI = 3.14159265358979323846;
static constexpr double HUGE_VAL_NUM = 1.0e30;

struct Prim {
    double rho{}, u{}, v{}, p{};
};

struct Cons {
    double r{}, ru{}, rv{}, rE{};

    Cons operator+(const Cons& o) const {
        return { r + o.r, ru + o.ru, rv + o.rv, rE + o.rE };
    }

    Cons operator-(const Cons& o) const {
        return { r - o.r, ru - o.ru, rv - o.rv, rE - o.rE };
    }

    Cons operator*(double s) const {
        return { r * s, ru * s, rv * s, rE * s };
    }

    Cons& operator+=(const Cons& o) {
        r += o.r;
        ru += o.ru;
        rv += o.rv;
        rE += o.rE;
        return *this;
    }
};

struct Vec2 {
    double x{}, y{};
};

struct CellGeom {
    double x{}, y{};
    bool solid{ false };
};

struct Panel {
    Vec2 a{}, b{}, n{};
    double length{};
};

struct Forces {
    double Fx{}, Fy{}, Drag{}, Lift{}, Cd{}, Cl{};
    double Mz_LE{}, Mz_c4{}, Cm_LE{}, Cm_c4{}, L_over_D{};
};

struct PanelForce {
    int panel_id{};
    std::string surface{};
    double x_mid{}, y_mid{}, nx{}, ny{}, length{}, pressure_avg{}, Cp_avg{};
    double Fx{}, Fy{}, Drag{}, Lift{}, Cd{}, Cl{}, Mz_LE{}, Mz_c4{}, Cm_LE{}, Cm_c4{};
    double avg_sample_distance{}, avg_valid_samples{};
    int fallback_count{};
};

struct PressureSample {
    double pressure = 101325.0;
    double Cp = 0.0;
    double sample_distance = 0.0;
    int valid_samples = 0;
    int used_fallback = 1;
};

struct ForceHistoryEntry {
    int iteration{};
    double residual_l2{}, residual_linf{}, normalized_residual_l2{}, normalized_residual_l2_ema{};
    double solution_change{}, solution_change_ema{};
    double Cd{}, Cl{}, Cm_LE{}, Cm_c4{}, L_over_D{};
};

struct ResidualInfo {
    double l2{}, linf{};
};

struct Diagnostics {
    int n_faces_x{}, n_faces_y{};
    int n_hllc_x{}, n_hllc_y{};
    int n_llf_x{}, n_llf_y{};
    int n_weno_x{}, n_weno_y{};
    int n_first_order_x{}, n_first_order_y{};
    int n_wall_x{}, n_wall_y{};
    int n_troubled_x{}, n_troubled_y{};
    int n_strong_shock_x{}, n_strong_shock_y{};
    int n_bad_reconstruction_x{}, n_bad_reconstruction_y{};
    int n_troubled_cells{};
};

struct DiagnosticsHistoryEntry {
    int iteration{};
    double dt{};
    Diagnostics d{};
};

// -----------------------------------------------------------------------------
// Solver configuration
// -----------------------------------------------------------------------------

struct Config {
    int nx = 720;
    int ny = 360;

    double xlo = -1.0, xhi = 3.0, ylo = -1.5, yhi = 1.5;

    double gamma = 1.4, gas_R = 287.0;
    double M_inf = 2.5, p_inf = 101325.0, T_inf = 288.15, aoa_deg = 5.0;


    double minimum_allowed_mach = 1.0001;
    double chord = 1.0, t_over_c = 0.10, x_le = 0.0, y_centre = 0.0;

    int max_iters = 40000;
    int min_iters = 8000;
    double cfl = 0.20;

    double solution_change_tol = 5.0e-6;
    double force_tol = 1.0e-4;
    int force_interval = 100;
    int force_window = 12;
    int diagnostics_interval = 100;

    double weno_eps = 1.0e-6;
    double trouble_ratio = 1.5e4;
    int trouble_lag = 5;

    double body_exclusion_chord = 0.12;
    int solid_check_radius = 1;

    bool use_sponge = true;
    double sponge_fraction = 0.10;
    double sponge_sigma = 8.0;

    double rho_floor = 1.0e-10;
    double p_floor = 1.0e-8;

    int force_samples_per_panel = 120;
    int surface_samples_per_panel = 240;


    int panel_force_skip_cells = 2;
    int panel_force_collect_cells = 4;
    bool use_panel_force_for_primary_coefficients = true;


    int wall_pressure_min_valid_samples = 3;
    int wall_pressure_max_valid_samples = 6;
    int wall_pressure_normal_rays = 12;
    int wall_pressure_tangent_radius = 1;

    double wall_pressure_first_offset_cells = 0.75;
    double wall_pressure_step_cells = 0.35;
    double wall_pressure_max_offset_cells = 8.0;
    double wall_pressure_idw_power = 2.0;


    double freestream_pressure_reject_rel = 2.0e-5;

    bool write_snapshots = false;
    int snapshot_interval = 500;
    std::string snapshot_prefix = "weno5_snap_";

    std::string output_dir = ".";
    std::string solution_csv = "weno5_solution.csv";
    std::string residual_csv = "weno5_residual.csv";
    std::string force_csv = "weno5_forces.csv";
    std::string force_panel_csv = "weno5_forces_panel.csv";
    std::string force_wall_pressure_csv = "weno5_forces_wall_pressure.csv";
    std::string force_comparison_csv = "force_comparison.csv";
    std::string force_history_csv = "weno5_force_history.csv";
    std::string panel_forces_csv = "panel_forces.csv";
    std::string surface_cp_csv = "surface_cp.csv";
    std::string wall_distance_csv = "wall_distance.csv";
    std::string diagnostics_csv = "weno5_diagnostics.csv";
};
