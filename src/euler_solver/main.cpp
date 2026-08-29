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

    int max_iters = 25000;
    int min_iters = 6000;
    double cfl = 0.20;

    double solution_change_tol = 5.0e-6;
    double force_tol = 5.0e-5;
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
    int snapshot_interval = 2000;
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

// -----------------------------------------------------------------------------
// File paths
// -----------------------------------------------------------------------------

static std::string make_output_path(const Config& c, const std::string& filename) {
    if (c.output_dir.empty() || c.output_dir == ".") {
        return filename;
    }
    std::filesystem::create_directories(c.output_dir);
    return (std::filesystem::path(c.output_dir) / filename).string();
}

// -----------------------------------------------------------------------------
// Basic utilities
// -----------------------------------------------------------------------------

static inline int idx(int i, int j, int nx) { return j * nx + i; }
static inline double deg_to_rad(double deg) { return deg * PI / 180.0; }
static inline double aoa_rad(const Config& c) { return deg_to_rad(c.aoa_deg); }
static inline double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
static inline double norm(Vec2 a) { return std::sqrt(dot(a, a)); }
static inline Vec2 add(Vec2 a, Vec2 b) { return { a.x + b.x, a.y + b.y }; }
static inline Vec2 sub(Vec2 a, Vec2 b) { return { a.x - b.x, a.y - b.y }; }
static inline Vec2 mul(Vec2 a, double s) { return { a.x * s, a.y * s }; }

static inline Vec2 unit(Vec2 a) {
    const double n = norm(a);
    return n < 1.0e-14 ? Vec2{ 0.0, 0.0 } : Vec2{ a.x / n, a.y / n };
}

// -----------------------------------------------------------------------------
// Flow state and thermodynamics
// -----------------------------------------------------------------------------

static Prim make_freestream(const Config& c) {
    const double alpha = aoa_rad(c);
    const double a = std::sqrt(c.gamma * c.gas_R * c.T_inf);
    const double V = c.M_inf * a;
    return { c.p_inf / (c.gas_R * c.T_inf), V * std::cos(alpha), V * std::sin(alpha), c.p_inf };
}

static inline double sound_speed(const Prim& W, double gamma) {
    return std::sqrt(std::max(gamma * W.p / std::max(W.rho, 1.0e-30), 1.0e-12));
}

static inline double velocity_magnitude(const Prim& W) {
    return std::sqrt(W.u * W.u + W.v * W.v);
}

static Cons prim_to_cons(const Prim& W, double gamma) {
    const double ke = 0.5 * (W.u * W.u + W.v * W.v);
    const double E = W.p / ((gamma - 1.0) * W.rho) + ke;
    return { W.rho, W.rho * W.u, W.rho * W.v, W.rho * E };
}

static Prim cons_to_prim(const Cons& U, double gamma) {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    if (!std::isfinite(U.r) || U.r <= 0.0) {
        return { nan, nan, nan, nan };
    }

    Prim W;
    W.rho = U.r;
    W.u = U.ru / W.rho;
    W.v = U.rv / W.rho;

    const double kinetic = 0.5 * (W.u * W.u + W.v * W.v);
    W.p = (gamma - 1.0) * (U.rE - W.rho * kinetic);

    return W;
}

static bool physical(const Prim& W, const Config& c) {
    return std::isfinite(W.rho)
        && std::isfinite(W.u)
        && std::isfinite(W.v)
        && std::isfinite(W.p)
        && W.rho > c.rho_floor
        && W.p > c.p_floor;
}

// -----------------------------------------------------------------------------
// Residual scaling
// -----------------------------------------------------------------------------

static Cons conservative_reference_scale(const Config& c) {
    const Prim W = make_freestream(c);
    const Cons U = prim_to_cons(W, c.gamma);
    const double rho_ref = std::max(std::abs(U.r), 1.0e-30);
    const double mom_ref = std::max(W.rho * velocity_magnitude(W), 1.0e-30);
    const double ene_ref = std::max(std::abs(U.rE), 1.0e-30);
    return { rho_ref, mom_ref, mom_ref, ene_ref };
}

// -----------------------------------------------------------------------------
// Analytical reference
// -----------------------------------------------------------------------------

static double prandtl_meyer_angle(double M, double gamma) {
    if (M <= 1.0) {
        return 0.0;
    }

    const double gm1 = gamma - 1.0;
    const double gp1 = gamma + 1.0;

    return std::sqrt(gp1 / gm1) * std::atan(std::sqrt(gm1 / gp1 * (M * M - 1.0)))
        - std::atan(std::sqrt(M * M - 1.0));
}

// -----------------------------------------------------------------------------
// Airfoil geometry
// -----------------------------------------------------------------------------

static std::array<Vec2, 4> diamond_vertices(const Config& c) {
    const double chord = c.chord;
    const double t = c.t_over_c * chord;

    return {
        Vec2{ c.x_le, c.y_centre },
        Vec2{ c.x_le + 0.5 * chord, c.y_centre + 0.5 * t },
        Vec2{ c.x_le + chord, c.y_centre },
        Vec2{ c.x_le + 0.5 * chord, c.y_centre - 0.5 * t }
    };
}

static bool inside_diamond(double x, double y, const Config& c) {
    const double xl = x - c.x_le;
    const double yl = y - c.y_centre;
    const double chord = c.chord;
    const double half_t = 0.5 * c.t_over_c * chord;

    if (xl < 0.0 || xl > chord) {
        return false;
    }

    const double yu = (xl <= 0.5 * chord)
        ? (half_t / (0.5 * chord)) * xl
        : (half_t / (0.5 * chord)) * (chord - xl);

    return std::abs(yl) <= yu;
}

static std::vector<CellGeom> build_geometry(const Config& c, double dx, double dy) {
    std::vector<CellGeom> geom(c.nx * c.ny);

    for (int j = 0; j < c.ny; ++j) {
        for (int i = 0; i < c.nx; ++i) {
            const int k = idx(i, j, c.nx);
            geom[k].x = c.xlo + (i + 0.5) * dx;
            geom[k].y = c.ylo + (j + 0.5) * dy;
            geom[k].solid = inside_diamond(geom[k].x, geom[k].y, c);
        }
    }

    return geom;
}

static std::vector<Panel> build_panels(const Config& c) {
    const auto V = diamond_vertices(c);
    Vec2 p[5] = { V[0], V[1], V[2], V[3], V[0] };
    const Vec2 centre{ c.x_le + 0.5 * c.chord, c.y_centre };

    std::vector<Panel> panels;
    panels.reserve(4);

    for (int k = 0; k < 4; ++k) {
        const Vec2 ab = sub(p[k + 1], p[k]);
        const double L = norm(ab);
        const Vec2 tangent = unit(ab);
        Vec2 n{ -tangent.y, tangent.x };

        const Vec2 mid{ 0.5 * (p[k].x + p[k + 1].x), 0.5 * (p[k].y + p[k + 1].y) };

        if (dot(n, sub(mid, centre)) < 0.0) {
            n.x = -n.x;
            n.y = -n.y;
        }

        panels.push_back({ p[k], p[k + 1], n, L });
    }

    return panels;
}

static std::string panel_surface_name(int id) {
    switch (id) {
    case 0: return "upper_forward";
    case 1: return "upper_rear";
    case 2: return "lower_rear";
    case 3: return "lower_forward";
    default: return "unknown";
    }
}

static double point_segment_distance(Vec2 p, Vec2 a, Vec2 b) {
    const Vec2 ab = sub(b, a);
    const Vec2 ap = sub(p, a);
    const double ab2 = dot(ab, ab);
    const double s = ab2 > 0.0 ? std::max(0.0, std::min(1.0, dot(ap, ab) / ab2)) : 0.0;
    return norm(sub(p, add(a, mul(ab, s))));
}

static double wall_distance(Vec2 p, const std::vector<Panel>& panels) {
    double best = HUGE_VAL_NUM;
    for (const auto& pan : panels) {
        best = std::min(best, point_segment_distance(p, pan.a, pan.b));
    }
    return best;
}

static Vec2 wall_normal(Vec2 p, const std::vector<Panel>& panels) {
    double best = HUGE_VAL_NUM;
    Vec2 best_n{ 1.0, 0.0 };

    for (const auto& pan : panels) {
        const Vec2 ab = sub(pan.b, pan.a);
        const Vec2 ap = sub(p, pan.a);
        const double ab2 = dot(ab, ab);
        const double s = ab2 > 0.0 ? std::max(0.0, std::min(1.0, dot(ap, ab) / ab2)) : 0.0;
        const Vec2 d = sub(p, add(pan.a, mul(ab, s)));
        const double d2 = dot(d, d);

        if (d2 < best) {
            best = d2;
            best_n = pan.n;
        }
    }

    return best_n;
}

static Prim reflected_state(const Prim& W, Vec2 n) {
    const double vn = W.u * n.x + W.v * n.y;
    return { W.rho, W.u - 2.0 * vn * n.x, W.v - 2.0 * vn * n.y, W.p };
}

// -----------------------------------------------------------------------------
// Euler fluxes
// -----------------------------------------------------------------------------

static Cons physical_flux(const Prim& W, double gamma, double nx, double ny) {
    const Cons U = prim_to_cons(W, gamma);
    const double un = W.u * nx + W.v * ny;

    return {
        U.r * un,
        U.ru * un + W.p * nx,
        U.rv * un + W.p * ny,
        (U.rE + W.p) * un
    };
}

static Cons llf_flux(const Prim& WL, const Prim& WR, double gamma, double nx, double ny) {
    const Cons UL = prim_to_cons(WL, gamma);
    const Cons UR = prim_to_cons(WR, gamma);
    const Cons FL = physical_flux(WL, gamma, nx, ny);
    const Cons FR = physical_flux(WR, gamma, nx, ny);

    const double aL = sound_speed(WL, gamma);
    const double aR = sound_speed(WR, gamma);
    const double unL = WL.u * nx + WL.v * ny;
    const double unR = WR.u * nx + WR.v * ny;
    const double smax = std::max(std::abs(unL) + aL, std::abs(unR) + aR);

    return (FL + FR) * 0.5 - (UR - UL) * (0.5 * smax);
}

static Cons hllc_flux(const Prim& WL, const Prim& WR, double gamma, double nx, double ny) {
    const Cons UL = prim_to_cons(WL, gamma);
    const Cons UR = prim_to_cons(WR, gamma);
    const Cons FL = physical_flux(WL, gamma, nx, ny);
    const Cons FR = physical_flux(WR, gamma, nx, ny);

    const double aL = sound_speed(WL, gamma);
    const double aR = sound_speed(WR, gamma);
    const double unL = WL.u * nx + WL.v * ny;
    const double unR = WR.u * nx + WR.v * ny;

    const double sL = std::min(unL - aL, unR - aR);
    const double sR = std::max(unL + aL, unR + aR);

    if (sL >= 0.0) {
        return FL;
    }

    if (sR <= 0.0) {
        return FR;
    }

    const double num = (WR.p - WL.p) + WL.rho * unL * (sL - unL) - WR.rho * unR * (sR - unR);
    const double den = WL.rho * (sL - unL) - WR.rho * (sR - unR);

    if (std::abs(den) < 1.0e-30) {
        return llf_flux(WL, WR, gamma, nx, ny);
    }

    const double sS = num / den;

    auto star = [&](const Prim& WK, const Cons& UK, const Cons& FK, double unK, double sK) -> Cons {
        const double wave = sK - unK;
        const double denom = sK - sS;

        if (std::abs(wave) < 1.0e-14 || std::abs(denom) < 1.0e-14) {
            return llf_flux(WL, WR, gamma, nx, ny);
        }

        const double factor = WK.rho * wave / denom;

        if (!(factor > 0.0) || !std::isfinite(factor)) {
            return llf_flux(WL, WR, gamma, nx, ny);
        }

        const double EK = UK.rE / WK.rho;
        const double pterm = WK.p / (WK.rho * wave);

        Cons Us{
            factor,
            factor * (WK.u + (sS - unK) * nx),
            factor * (WK.v + (sS - unK) * ny),
            factor * (EK + (sS - unK) * (sS + pterm))
        };

        if (!std::isfinite(Us.r) || !std::isfinite(Us.ru) || !std::isfinite(Us.rv) || !std::isfinite(Us.rE)) {
            return llf_flux(WL, WR, gamma, nx, ny);
        }

        return FK + (Us - UK) * sK;
    };

    return sS >= 0.0 ? star(WL, UL, FL, unL, sL) : star(WR, UR, FR, unR, sR);
}

// -----------------------------------------------------------------------------
// WENO5 reconstruction and shock detection
// -----------------------------------------------------------------------------

static bool strong_shock_face(const Prim& WL, const Prim& WR, const Config& c) {
    const double pmin = std::max(std::min(WL.p, WR.p), c.p_floor);
    const double rmin = std::max(std::min(WL.rho, WR.rho), c.rho_floor);
    const double pressure_jump = std::abs(WR.p - WL.p) / pmin;
    const double density_jump = std::abs(WR.rho - WL.rho) / rmin;
    const double ML = velocity_magnitude(WL) / std::max(sound_speed(WL, c.gamma), 1.0e-12);
    const double MR = velocity_magnitude(WR) / std::max(sound_speed(WR, c.gamma), 1.0e-12);
    const double mach_jump = std::abs(MR - ML);
    return pressure_jump > 2.50 || density_jump > 1.50 || mach_jump > 1.25;
}

static inline double sqr(double x) { return x * x; }

static inline double weno5_left(double vm2, double vm1, double vi, double vp1, double vp2, double eps) {
    const double q0 = (1.0 / 3.0) * vm2 - (7.0 / 6.0) * vm1 + (11.0 / 6.0) * vi;
    const double q1 = -(1.0 / 6.0) * vm1 + (5.0 / 6.0) * vi + (1.0 / 3.0) * vp1;
    const double q2 = (1.0 / 3.0) * vi + (5.0 / 6.0) * vp1 - (1.0 / 6.0) * vp2;

    const double d00 = vm2 - 2.0 * vm1 + vi;
    const double d01 = vm2 - 4.0 * vm1 + 3.0 * vi;
    const double d10 = vm1 - 2.0 * vi + vp1;
    const double d11 = vm1 - vp1;
    const double d20 = vi - 2.0 * vp1 + vp2;
    const double d21 = 3.0 * vi - 4.0 * vp1 + vp2;

    const double b0 = (13.0 / 12.0) * d00 * d00 + 0.25 * d01 * d01;
    const double b1 = (13.0 / 12.0) * d10 * d10 + 0.25 * d11 * d11;
    const double b2 = (13.0 / 12.0) * d20 * d20 + 0.25 * d21 * d21;

    const double a0 = 0.1 / sqr(eps + b0);
    const double a1 = 0.6 / sqr(eps + b1);
    const double a2 = 0.3 / sqr(eps + b2);
    const double s = a0 + a1 + a2;

    return (a0 * q0 + a1 * q1 + a2 * q2) / s;
}

static inline double weno5_right(double vm1, double vi, double vp1, double vp2, double vp3, double eps) {
    return weno5_left(vp3, vp2, vp1, vi, vm1, eps);
}

static inline double beta_ratio(double vm2, double vm1, double vi, double vp1, double vp2, double eps) {
    const double d00 = vm2 - 2.0 * vm1 + vi;
    const double d01 = vm2 - 4.0 * vm1 + 3.0 * vi;
    const double d10 = vm1 - 2.0 * vi + vp1;
    const double d11 = vm1 - vp1;
    const double d20 = vi - 2.0 * vp1 + vp2;
    const double d21 = 3.0 * vi - 4.0 * vp1 + vp2;

    const double b0 = (13.0 / 12.0) * d00 * d00 + 0.25 * d01 * d01;
    const double b1 = (13.0 / 12.0) * d10 * d10 + 0.25 * d11 * d11;
    const double b2 = (13.0 / 12.0) * d20 * d20 + 0.25 * d21 * d21;

    return std::max({ b0, b1, b2 }) / (std::min({ b0, b1, b2 }) + eps);
}

// -----------------------------------------------------------------------------
// Immersed boundary
// -----------------------------------------------------------------------------

struct FaceData {
    int nx{}, ny{};
    std::vector<bool> ok_x, ok_y, wall_x, wall_y, fL_x, fR_x, fB_y, fT_y;
    std::vector<Vec2> nrm_x, nrm_y;

    void build(const Config& c, const std::vector<CellGeom>& geom, const std::vector<Panel>& panels, double dx, double dy) {
        nx = c.nx;
        ny = c.ny;

        const int nxf = (nx + 1) * ny;
        const int nyf = nx * (ny + 1);

        ok_x.assign(nxf, false);
        ok_y.assign(nyf, false);
        wall_x.assign(nxf, false);
        wall_y.assign(nyf, false);
        fL_x.assign(nxf, false);
        fR_x.assign(nxf, false);
        fB_y.assign(nyf, false);
        fT_y.assign(nyf, false);
        nrm_x.assign(nxf, { 1.0, 0.0 });
        nrm_y.assign(nyf, { 0.0, 1.0 });

        const int R = c.solid_check_radius;

        for (int j = 0; j < ny; ++j) {
            for (int iface = 0; iface <= nx; ++iface) {
                const int fi = iface * ny + j;
                const bool fL = (iface > 0) && !geom[idx(iface - 1, j, nx)].solid;
                const bool fR = (iface < nx) && !geom[idx(iface, j, nx)].solid;
                fL_x[fi] = fL;
                fR_x[fi] = fR;

                if (!fL && !fR) {
                    continue;
                }

                const Vec2 fc{ c.xlo + iface * dx, c.ylo + (j + 0.5) * dy };

                if (fL && fR) {
                    if (wall_distance(fc, panels) < c.body_exclusion_chord * c.chord) {
                        continue;
                    }

                    bool clear = true;
                    for (int dj = -R; dj <= R && clear; ++dj) {
                        for (int di = -R; di <= R && clear; ++di) {
                            const int ii = iface + di;
                            const int jj = j + dj;
                            if (ii < 0 || ii >= nx || jj < 0 || jj >= ny) {
                                continue;
                            }
                            if (geom[idx(ii, jj, nx)].solid) {
                                clear = false;
                            }
                        }
                    }
                    ok_x[fi] = clear;
                }
                else {
                    wall_x[fi] = true;
                    nrm_x[fi] = wall_normal(fc, panels);
                }
            }
        }

        for (int i = 0; i < nx; ++i) {
            for (int jf = 0; jf <= ny; ++jf) {
                const int fi = i * (ny + 1) + jf;
                const bool fB = (jf > 0) && !geom[idx(i, jf - 1, nx)].solid;
                const bool fT = (jf < ny) && !geom[idx(i, jf, nx)].solid;
                fB_y[fi] = fB;
                fT_y[fi] = fT;

                if (!fB && !fT) {
                    continue;
                }

                const Vec2 fc{ c.xlo + (i + 0.5) * dx, c.ylo + jf * dy };

                if (fB && fT) {
                    if (wall_distance(fc, panels) < c.body_exclusion_chord * c.chord) {
                        continue;
                    }

                    bool clear = true;
                    for (int dj = -R; dj <= R && clear; ++dj) {
                        for (int di = -R; di <= R && clear; ++di) {
                            const int ii = i + di;
                            const int jj = jf + dj;
                            if (ii < 0 || ii >= nx || jj < 0 || jj >= ny) {
                                continue;
                            }
                            if (geom[idx(ii, jj, nx)].solid) {
                                clear = false;
                            }
                        }
                    }
                    ok_y[fi] = clear;
                }
                else {
                    wall_y[fi] = true;
                    nrm_y[fi] = wall_normal(fc, panels);
                }
            }
        }
    }

    void apply_trouble_lag(const std::vector<int>& lag) {
        for (int j = 0; j < ny; ++j) {
            for (int iface = 1; iface < nx; ++iface) {
                const int fi = iface * ny + j;
                if (!ok_x[fi]) {
                    continue;
                }
                if (lag[idx(iface - 1, j, nx)] > 0 || lag[idx(iface, j, nx)] > 0) {
                    ok_x[fi] = false;
                }
            }
        }

        for (int i = 0; i < nx; ++i) {
            for (int jf = 1; jf < ny; ++jf) {
                const int fi = i * (ny + 1) + jf;
                if (!ok_y[fi]) {
                    continue;
                }
                if (lag[idx(i, jf - 1, nx)] > 0 || lag[idx(i, jf, nx)] > 0) {
                    ok_y[fi] = false;
                }
            }
        }
    }
};

// -----------------------------------------------------------------------------
// Primitive state cache
// -----------------------------------------------------------------------------

static void update_primitive_cache(
    std::vector<Prim>& W,
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const Config& c,
    const Prim& Winf
) {
    const int N = c.nx * c.ny;
    for (int k = 0; k < N; ++k) {
        W[k] = geom[k].solid ? Winf : cons_to_prim(U[k], c.gamma);
    }
}

// -----------------------------------------------------------------------------
// Spatial residual
// -----------------------------------------------------------------------------

static ResidualInfo compute_rhs(
    const std::vector<Prim>& W,
    const std::vector<CellGeom>& geom,
    const Prim& Winf,
    const FaceData& fd,
    std::vector<Cons>& RHS,
    std::vector<bool>& troubled_x,
    std::vector<bool>& troubled_y,
    const Config& c,
    double dx,
    double dy,
    Diagnostics* diag = nullptr
) {
    const int nx = c.nx;
    const int ny = c.ny;

    std::fill(RHS.begin(), RHS.end(), Cons{ 0.0, 0.0, 0.0, 0.0 });
    std::fill(troubled_x.begin(), troubled_x.end(), false);
    std::fill(troubled_y.begin(), troubled_y.end(), false);

    if (diag) {
        *diag = Diagnostics{};
    }

    auto getx = [&](int i, int j) -> const Prim& {
        i = std::max(0, std::min(nx - 1, i));
        return W[idx(i, j, nx)];
    };

    auto gety = [&](int i, int j) -> const Prim& {
        j = std::max(0, std::min(ny - 1, j));
        return W[idx(i, j, nx)];
    };

    for (int j = 0; j < ny; ++j) {
        for (int iface = 0; iface <= nx; ++iface) {
            Prim WL, WR;
            bool use_llf = false;
            bool used_weno = false;
            bool wall = false;
            bool strong = false;
            bool badrec = false;

            const int fi = iface * ny + j;

            if (iface == 0) {
                if (!fd.fR_x[fi]) {
                    continue;
                }
                WL = Winf;
                WR = W[idx(0, j, nx)];
            }
            else if (iface == nx) {
                if (!fd.fL_x[fi]) {
                    continue;
                }
                WL = W[idx(nx - 1, j, nx)];
                WR = WL;
            }
            else {
                const bool fL = fd.fL_x[fi];
                const bool fR = fd.fR_x[fi];

                if (!fL && !fR) {
                    continue;
                }

                if (fd.wall_x[fi]) {
                    wall = true;
                    if (fL) {
                        WL = W[idx(iface - 1, j, nx)];
                        WR = reflected_state(WL, fd.nrm_x[fi]);
                    }
                    else {
                        WR = W[idx(iface, j, nx)];
                        WL = reflected_state(WR, fd.nrm_x[fi]);
                    }
                }
                else if (fd.ok_x[fi]) {
                    const Prim& Lm2 = getx(iface - 3, j);
                    const Prim& Lm1 = getx(iface - 2, j);
                    const Prim& Li = getx(iface - 1, j);
                    const Prim& Lp1 = getx(iface, j);
                    const Prim& Lp2 = getx(iface + 1, j);

                    const Prim& Rm1 = getx(iface - 2, j);
                    const Prim& Ri = getx(iface - 1, j);
                    const Prim& Rp1 = getx(iface, j);
                    const Prim& Rp2 = getx(iface + 1, j);
                    const Prim& Rp3 = getx(iface + 2, j);

                    const double br = beta_ratio(Lm2.p, Lm1.p, Li.p, Lp1.p, Lp2.p, c.weno_eps);

                    if (br > c.trouble_ratio) {
                        WL = W[idx(iface - 1, j, nx)];
                        WR = W[idx(iface, j, nx)];
                        use_llf = true;
                        troubled_x[fi] = true;
                    }
                    else {
                        used_weno = true;
                        WL = {
                            weno5_left(Lm2.rho, Lm1.rho, Li.rho, Lp1.rho, Lp2.rho, c.weno_eps),
                            weno5_left(Lm2.u, Lm1.u, Li.u, Lp1.u, Lp2.u, c.weno_eps),
                            weno5_left(Lm2.v, Lm1.v, Li.v, Lp1.v, Lp2.v, c.weno_eps),
                            weno5_left(Lm2.p, Lm1.p, Li.p, Lp1.p, Lp2.p, c.weno_eps)
                        };

                        WR = {
                            weno5_right(Rm1.rho, Ri.rho, Rp1.rho, Rp2.rho, Rp3.rho, c.weno_eps),
                            weno5_right(Rm1.u, Ri.u, Rp1.u, Rp2.u, Rp3.u, c.weno_eps),
                            weno5_right(Rm1.v, Ri.v, Rp1.v, Rp2.v, Rp3.v, c.weno_eps),
                            weno5_right(Rm1.p, Ri.p, Rp1.p, Rp2.p, Rp3.p, c.weno_eps)
                        };

                        if (!physical(WL, c) || !physical(WR, c)) {
                            WL = W[idx(iface - 1, j, nx)];
                            WR = W[idx(iface, j, nx)];
                            use_llf = true;
                            troubled_x[fi] = true;
                            badrec = true;
                        }
                    }
                }
                else {
                    WL = W[idx(iface - 1, j, nx)];
                    WR = W[idx(iface, j, nx)];
                }
            }

            if (!use_llf && strong_shock_face(WL, WR, c)) {
                use_llf = true;
                strong = true;
            }

            if (diag) {
                diag->n_faces_x++;
                if (use_llf) diag->n_llf_x++; else diag->n_hllc_x++;
                if (used_weno) diag->n_weno_x++; else diag->n_first_order_x++;
                if (wall) diag->n_wall_x++;
                if (troubled_x[fi]) diag->n_troubled_x++;
                if (strong) diag->n_strong_shock_x++;
                if (badrec) diag->n_bad_reconstruction_x++;
            }

            const Cons F = (use_llf ? llf_flux(WL, WR, c.gamma, 1.0, 0.0) : hllc_flux(WL, WR, c.gamma, 1.0, 0.0)) * (1.0 / dx);

            if (iface > 0 && fd.fL_x[fi]) {
                RHS[idx(iface - 1, j, nx)] += F * (-1.0);
            }
            if (iface < nx && fd.fR_x[fi]) {
                RHS[idx(iface, j, nx)] += F;
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int jf = 0; jf <= ny; ++jf) {
            Prim WB, WT;
            bool use_llf = false;
            bool used_weno = false;
            bool wall = false;
            bool strong = false;
            bool badrec = false;

            const int fi = i * (ny + 1) + jf;

            if (jf == 0) {
                if (!fd.fT_y[fi]) {
                    continue;
                }
                WB = Winf;
                WT = W[idx(i, 0, nx)];
            }
            else if (jf == ny) {
                if (!fd.fB_y[fi]) {
                    continue;
                }
                WB = W[idx(i, ny - 1, nx)];
                WT = Winf;
            }
            else {
                const bool fB = fd.fB_y[fi];
                const bool fT = fd.fT_y[fi];

                if (!fB && !fT) {
                    continue;
                }

                if (fd.wall_y[fi]) {
                    wall = true;
                    if (fB) {
                        WB = W[idx(i, jf - 1, nx)];
                        WT = reflected_state(WB, fd.nrm_y[fi]);
                    }
                    else {
                        WT = W[idx(i, jf, nx)];
                        WB = reflected_state(WT, fd.nrm_y[fi]);
                    }
                }
                else if (fd.ok_y[fi]) {
                    const Prim& Bm2 = gety(i, jf - 3);
                    const Prim& Bm1 = gety(i, jf - 2);
                    const Prim& Bi = gety(i, jf - 1);
                    const Prim& Bp1 = gety(i, jf);
                    const Prim& Bp2 = gety(i, jf + 1);

                    const Prim& Tm1 = gety(i, jf - 2);
                    const Prim& Ti = gety(i, jf - 1);
                    const Prim& Tp1 = gety(i, jf);
                    const Prim& Tp2 = gety(i, jf + 1);
                    const Prim& Tp3 = gety(i, jf + 2);

                    const double br = beta_ratio(Bm2.p, Bm1.p, Bi.p, Bp1.p, Bp2.p, c.weno_eps);

                    if (br > c.trouble_ratio) {
                        WB = W[idx(i, jf - 1, nx)];
                        WT = W[idx(i, jf, nx)];
                        use_llf = true;
                        troubled_y[fi] = true;
                    }
                    else {
                        used_weno = true;
                        WB = {
                            weno5_left(Bm2.rho, Bm1.rho, Bi.rho, Bp1.rho, Bp2.rho, c.weno_eps),
                            weno5_left(Bm2.u, Bm1.u, Bi.u, Bp1.u, Bp2.u, c.weno_eps),
                            weno5_left(Bm2.v, Bm1.v, Bi.v, Bp1.v, Bp2.v, c.weno_eps),
                            weno5_left(Bm2.p, Bm1.p, Bi.p, Bp1.p, Bp2.p, c.weno_eps)
                        };

                        WT = {
                            weno5_right(Tm1.rho, Ti.rho, Tp1.rho, Tp2.rho, Tp3.rho, c.weno_eps),
                            weno5_right(Tm1.u, Ti.u, Tp1.u, Tp2.u, Tp3.u, c.weno_eps),
                            weno5_right(Tm1.v, Ti.v, Tp1.v, Tp2.v, Tp3.v, c.weno_eps),
                            weno5_right(Tm1.p, Ti.p, Tp1.p, Tp2.p, Tp3.p, c.weno_eps)
                        };

                        if (!physical(WB, c) || !physical(WT, c)) {
                            WB = W[idx(i, jf - 1, nx)];
                            WT = W[idx(i, jf, nx)];
                            use_llf = true;
                            troubled_y[fi] = true;
                            badrec = true;
                        }
                    }
                }
                else {
                    WB = W[idx(i, jf - 1, nx)];
                    WT = W[idx(i, jf, nx)];
                }
            }

            if (!use_llf && strong_shock_face(WB, WT, c)) {
                use_llf = true;
                strong = true;
            }

            if (diag) {
                diag->n_faces_y++;
                if (use_llf) diag->n_llf_y++; else diag->n_hllc_y++;
                if (used_weno) diag->n_weno_y++; else diag->n_first_order_y++;
                if (wall) diag->n_wall_y++;
                if (troubled_y[fi]) diag->n_troubled_y++;
                if (strong) diag->n_strong_shock_y++;
                if (badrec) diag->n_bad_reconstruction_y++;
            }

            const Cons F = (use_llf ? llf_flux(WB, WT, c.gamma, 0.0, 1.0) : hllc_flux(WB, WT, c.gamma, 0.0, 1.0)) * (1.0 / dy);

            if (jf > 0 && fd.fB_y[fi]) {
                RHS[idx(i, jf - 1, nx)] += F * (-1.0);
            }
            if (jf < ny && fd.fT_y[fi]) {
                RHS[idx(i, jf, nx)] += F;
            }
        }
    }

    ResidualInfo info{};
    double sum2 = 0.0;
    int count = 0;

    for (int k = 0; k < nx * ny; ++k) {
        if (geom[k].solid) {
            continue;
        }

        const double v = std::sqrt(
            RHS[k].r * RHS[k].r
            + RHS[k].ru * RHS[k].ru
            + RHS[k].rv * RHS[k].rv
            + RHS[k].rE * RHS[k].rE
        );

        info.linf = std::max(info.linf, v);
        sum2 += v * v;
        ++count;
    }

    info.l2 = std::sqrt(sum2 / std::max(count, 1));
    return info;
}

// -----------------------------------------------------------------------------
// Time step
// -----------------------------------------------------------------------------

static double compute_dt(const std::vector<Prim>& W, const std::vector<CellGeom>& geom, const Config& c, double dx, double dy) {
    double dt = HUGE_VAL_NUM;

    for (int j = 0; j < c.ny; ++j) {
        for (int i = 0; i < c.nx; ++i) {
            const int k = idx(i, j, c.nx);
            if (geom[k].solid) {
                continue;
            }

            const double a = sound_speed(W[k], c.gamma);
            const double spectral = (std::abs(W[k].u) + a) / dx + (std::abs(W[k].v) + a) / dy;
            dt = std::min(dt, c.cfl / std::max(spectral, 1.0e-30));
        }
    }

    return dt;
}

// -----------------------------------------------------------------------------
// Far-field damping
// -----------------------------------------------------------------------------

static void apply_sponge(std::vector<Cons>& U, const std::vector<CellGeom>& geom, const Cons& Uinf, const Config& c, double dx, double dy, double dt) {
    if (!c.use_sponge) {
        return;
    }

    const double xs = c.xhi - (c.xhi - c.xlo) * c.sponge_fraction;
    const double xlen = (c.xhi - c.xlo) * c.sponge_fraction;
    const double ylen = (c.yhi - c.ylo) * c.sponge_fraction;
    const double yb = c.ylo + ylen;
    const double yt = c.yhi - ylen;

    for (int j = 0; j < c.ny; ++j) {
        for (int i = 0; i < c.nx; ++i) {
            const int k = idx(i, j, c.nx);
            if (geom[k].solid) {
                continue;
            }

            const double x = c.xlo + (i + 0.5) * dx;
            const double y = c.ylo + (j + 0.5) * dy;
            double sigma = 0.0;

            if (x > xs) {
                const double xi = (x - xs) / xlen;
                sigma += c.sponge_sigma * xi * xi;
            }

            if (y < yb) {
                const double eta = (yb - y) / ylen;
                sigma += 0.5 * c.sponge_sigma * eta * eta;
            }

            if (y > yt) {
                const double eta = (y - yt) / ylen;
                sigma += 0.5 * c.sponge_sigma * eta * eta;
            }

            if (sigma <= 0.0) {
                continue;
            }

            const double f = std::exp(-sigma * dt);
            U[k].r = Uinf.r + (U[k].r - Uinf.r) * f;
            U[k].ru = Uinf.ru + (U[k].ru - Uinf.ru) * f;
            U[k].rv = Uinf.rv + (U[k].rv - Uinf.rv) * f;
            U[k].rE = Uinf.rE + (U[k].rE - Uinf.rE) * f;
        }
    }
}

// -----------------------------------------------------------------------------
// Wall pressure sampling
// -----------------------------------------------------------------------------

static PressureSample sample_wall_pressure(
    Vec2 surface_point,
    Vec2 normal,
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const Config& c,
    double dx,
    double dy
) {
    const Prim Winf = make_freestream(c);
    const double Vinf = velocity_magnitude(Winf);
    const double qinf = 0.5 * Winf.rho * Vinf * Vinf;

    const double h = std::min(dx, dy);
    const double first_offset = c.wall_pressure_first_offset_cells * h;
    const double step = c.wall_pressure_step_cells * h;
    const double max_offset = c.wall_pressure_max_offset_cells * h;

    const double p_reject = c.freestream_pressure_reject_rel * std::max(c.p_inf, 1.0);

    struct Candidate {
        double p = 0.0;
        double distance = 0.0;
        double weight = 0.0;
        bool freestream_like = false;
    };

    std::vector<Candidate> candidates;
    candidates.reserve(128);

    const Vec2 tangent = { -normal.y, normal.x };

    for (int ray = 0; ray < c.wall_pressure_normal_rays; ++ray) {
        const double s = first_offset + static_cast<double>(ray) * step;

        if (s > max_offset) {
            break;
        }

        const Vec2 ray_point = {
            surface_point.x + s * normal.x,
            surface_point.y + s * normal.y
        };

        for (int tj = -c.wall_pressure_tangent_radius; tj <= c.wall_pressure_tangent_radius; ++tj) {
            const Vec2 shifted = {
                ray_point.x + static_cast<double>(tj) * h * tangent.x,
                ray_point.y + static_cast<double>(tj) * h * tangent.y
            };

            const int ii0 = static_cast<int>((shifted.x - c.xlo) / dx);
            const int jj0 = static_cast<int>((shifted.y - c.ylo) / dy);

            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    const int ii = ii0 + di;
                    const int jj = jj0 + dj;

                    if (ii < 0 || ii >= c.nx || jj < 0 || jj >= c.ny) {
                        continue;
                    }

                    const int k = idx(ii, jj, c.nx);

                    if (geom[k].solid) {
                        continue;
                    }

                    const Prim Wk = cons_to_prim(U[k], c.gamma);

                    if (!physical(Wk, c)) {
                        continue;
                    }

                    const double rx = geom[k].x - surface_point.x;
                    const double ry = geom[k].y - surface_point.y;
                    const double dist = std::sqrt(rx * rx + ry * ry);

                    if (!(dist > 1.0e-14)) {
                        continue;
                    }

                    const double normal_distance = rx * normal.x + ry * normal.y;

                    if (normal_distance <= 0.15 * h) {
                        continue;
                    }

                    const bool freestream_like = std::abs(Wk.p - c.p_inf) <= p_reject;

                    const double w = 1.0 / std::pow(
                        std::max(dist, 0.25 * h),
                        c.wall_pressure_idw_power
                    );

                    candidates.push_back({ Wk.p, dist, w, freestream_like });
                }
            }
        }
    }

    std::vector<Candidate> good;
    good.reserve(candidates.size());

    for (const auto& q : candidates) {
        if (!q.freestream_like) {
            good.push_back(q);
        }
    }

    const bool used_candidate_fallback =
        static_cast<int>(good.size()) < c.wall_pressure_min_valid_samples;

    const std::vector<Candidate>* use = used_candidate_fallback ? &candidates : &good;

    if (use->empty()) {
        return { c.p_inf, 0.0, 0.0, 0, 1 };
    }

    std::vector<Candidate> selected = *use;

    std::sort(
        selected.begin(),
        selected.end(),
        [](const Candidate& a, const Candidate& b) {
            return a.distance < b.distance;
        }
    );

    if (static_cast<int>(selected.size()) > c.wall_pressure_max_valid_samples) {
        selected.resize(c.wall_pressure_max_valid_samples);
    }

    double wsum = 0.0;
    double psum = 0.0;
    double dsum = 0.0;

    for (const auto& q : selected) {
        wsum += q.weight;
        psum += q.weight * q.p;
        dsum += q.weight * q.distance;
    }

    if (!(wsum > 0.0)) {
        return { c.p_inf, 0.0, 0.0, 0, 1 };
    }

    const double p_wall = psum / wsum;
    const double d_avg = dsum / wsum;
    const double Cp = (p_wall - c.p_inf) / qinf;

    const int fallback = used_candidate_fallback ? 1 : 0;

    return { p_wall, Cp, d_avg, static_cast<int>(selected.size()), fallback };
}

// -----------------------------------------------------------------------------
// Panel force calculation
// -----------------------------------------------------------------------------

static double sample_panel_pressure(
    Vec2 mid,
    Vec2 normal,
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const Config& c,
    double dx,
    double dy
) {
    const double step = std::min(dx, dy);

    std::vector<double> samples;
    samples.reserve(c.panel_force_collect_cells);

    for (int s = 0;
         s < c.panel_force_skip_cells + c.panel_force_collect_cells
             && static_cast<int>(samples.size()) < c.panel_force_collect_cells;
         ++s) {
        const Vec2 pt = {
            mid.x + (static_cast<double>(s) + 0.5) * step * normal.x,
            mid.y + (static_cast<double>(s) + 0.5) * step * normal.y
        };

        const int ic = static_cast<int>((pt.x - c.xlo) / dx);
        const int jc = static_cast<int>((pt.y - c.ylo) / dy);

        double best_dist = HUGE_VAL_NUM;
        double best_p = -1.0;

        for (int dj = -1; dj <= 1; ++dj) {
            for (int di = -1; di <= 1; ++di) {
                const int ii = ic + di;
                const int jj = jc + dj;

                if (ii < 0 || ii >= c.nx || jj < 0 || jj >= c.ny) {
                    continue;
                }

                const int k = idx(ii, jj, c.nx);

                if (geom[k].solid) {
                    continue;
                }

                const double ex = geom[k].x - pt.x;
                const double ey = geom[k].y - pt.y;
                const double d2 = ex * ex + ey * ey;

                if (d2 < best_dist) {
                    const Prim Wk = cons_to_prim(U[k], c.gamma);

                    if (physical(Wk, c)) {
                        best_dist = d2;
                        best_p = Wk.p;
                    }
                }
            }
        }

        if (best_p > 0.0 && s >= c.panel_force_skip_cells) {
            samples.push_back(best_p);
        }
    }

    if (samples.empty()) {
        return c.p_inf;
    }

    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

static Forces compute_forces_panel_sampling(
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const std::vector<Panel>& panels,
    const Config& c,
    double dx,
    double dy,
    std::vector<PanelForce>* pfo = nullptr
) {
    const Prim Winf = make_freestream(c);
    const double Vinf = velocity_magnitude(Winf);
    const double q = 0.5 * Winf.rho * Vinf * Vinf;
    const double alpha = aoa_rad(c);

    double Fx = 0.0;
    double Fy = 0.0;
    double Mle = 0.0;
    double Mc4 = 0.0;

    if (pfo) {
        pfo->clear();
        pfo->reserve(panels.size());
    }

    for (std::size_t pid = 0; pid < panels.size(); ++pid) {
        const Panel& pan = panels[pid];

        const Vec2 mid = {
            0.5 * (pan.a.x + pan.b.x),
            0.5 * (pan.a.y + pan.b.y)
        };

        const double p = sample_panel_pressure(mid, pan.n, U, geom, c, dx, dy);

        const double pFx = -p * pan.n.x * pan.length;
        const double pFy = -p * pan.n.y * pan.length;

        const double pMle = (mid.x - c.x_le) * pFy - (mid.y - c.y_centre) * pFx;
        const double pMc4 = (mid.x - (c.x_le + 0.25 * c.chord)) * pFy
                          - (mid.y - c.y_centre) * pFx;

        Fx += pFx;
        Fy += pFy;
        Mle += pMle;
        Mc4 += pMc4;

        if (pfo) {
            const double Dp = pFx * std::cos(alpha) + pFy * std::sin(alpha);
            const double Lp = -pFx * std::sin(alpha) + pFy * std::cos(alpha);

            pfo->push_back({
                static_cast<int>(pid),
                panel_surface_name(static_cast<int>(pid)),
                mid.x,
                mid.y,
                pan.n.x,
                pan.n.y,
                pan.length,
                p,
                (p - c.p_inf) / q,
                pFx,
                pFy,
                Dp,
                Lp,
                Dp / (q * c.chord),
                Lp / (q * c.chord),
                pMle,
                pMc4,
                pMle / (q * c.chord * c.chord),
                pMc4 / (q * c.chord * c.chord),
                static_cast<double>(c.panel_force_skip_cells),
                static_cast<double>(c.panel_force_collect_cells),
                0
            });
        }
    }

    const double D = Fx * std::cos(alpha) + Fy * std::sin(alpha);
    const double L = -Fx * std::sin(alpha) + Fy * std::cos(alpha);

    return {
        Fx,
        Fy,
        D,
        L,
        D / (q * c.chord),
        L / (q * c.chord),
        Mle,
        Mc4,
        Mle / (q * c.chord * c.chord),
        Mc4 / (q * c.chord * c.chord),
        std::abs(D) > 1.0e-30 ? L / D : 0.0
    };
}

// -----------------------------------------------------------------------------
// Wall pressure force calculation
// -----------------------------------------------------------------------------

static Forces compute_forces_wall_pressure(
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const std::vector<Panel>& panels,
    const Config& c,
    double dx,
    double dy,
    std::vector<PanelForce>* pfo = nullptr
) {
    const Prim Winf = make_freestream(c);
    const double Vinf = velocity_magnitude(Winf);
    const double q = 0.5 * Winf.rho * Vinf * Vinf;
    const double alpha = aoa_rad(c);

    double Fx = 0.0;
    double Fy = 0.0;
    double Mle = 0.0;
    double Mc4 = 0.0;

    const int ns = std::max(1, c.force_samples_per_panel);

    if (pfo) {
        pfo->clear();
        pfo->reserve(panels.size());
    }

    for (std::size_t pid = 0; pid < panels.size(); ++pid) {
        const Panel& pan = panels[pid];
        const Vec2 ab = sub(pan.b, pan.a);
        const double dL = pan.length / static_cast<double>(ns);

        double pFx = 0.0;
        double pFy = 0.0;
        double pMle = 0.0;
        double pMc4 = 0.0;
        double pInt = 0.0;
        double sample_distance_sum = 0.0;
        double valid_samples_sum = 0.0;
        int fallback_count = 0;

        for (int s = 0; s < ns; ++s) {
            const double eta = (static_cast<double>(s) + 0.5) / static_cast<double>(ns);
            const Vec2 sp = { pan.a.x + eta * ab.x, pan.a.y + eta * ab.y };

            const PressureSample ps = sample_wall_pressure(sp, pan.n, U, geom, c, dx, dy);
            const double p = ps.pressure;

            const double dFx = -p * pan.n.x * dL;
            const double dFy = -p * pan.n.y * dL;

            pFx += dFx;
            pFy += dFy;
            pInt += p * dL;

            pMle += (sp.x - c.x_le) * dFy - (sp.y - c.y_centre) * dFx;
            pMc4 += (sp.x - (c.x_le + 0.25 * c.chord)) * dFy - (sp.y - c.y_centre) * dFx;

            sample_distance_sum += ps.sample_distance;
            valid_samples_sum += ps.valid_samples;
            fallback_count += ps.used_fallback;
        }

        Fx += pFx;
        Fy += pFy;
        Mle += pMle;
        Mc4 += pMc4;

        if (pfo) {
            const Vec2 mid = { 0.5 * (pan.a.x + pan.b.x), 0.5 * (pan.a.y + pan.b.y) };
            const double pavg = pInt / pan.length;
            const double D = pFx * std::cos(alpha) + pFy * std::sin(alpha);
            const double L = -pFx * std::sin(alpha) + pFy * std::cos(alpha);

            pfo->push_back({
                static_cast<int>(pid),
                panel_surface_name(static_cast<int>(pid)),
                mid.x,
                mid.y,
                pan.n.x,
                pan.n.y,
                pan.length,
                pavg,
                (pavg - c.p_inf) / q,
                pFx,
                pFy,
                D,
                L,
                D / (q * c.chord),
                L / (q * c.chord),
                pMle,
                pMc4,
                pMle / (q * c.chord * c.chord),
                pMc4 / (q * c.chord * c.chord),
                sample_distance_sum / static_cast<double>(ns),
                valid_samples_sum / static_cast<double>(ns),
                fallback_count
            });
        }
    }

    const double D = Fx * std::cos(alpha) + Fy * std::sin(alpha);
    const double L = -Fx * std::sin(alpha) + Fy * std::cos(alpha);

    return {
        Fx,
        Fy,
        D,
        L,
        D / (q * c.chord),
        L / (q * c.chord),
        Mle,
        Mc4,
        Mle / (q * c.chord * c.chord),
        Mc4 / (q * c.chord * c.chord),
        std::abs(D) > 1.0e-30 ? L / D : 0.0
    };
}

// -----------------------------------------------------------------------------
// Force method
// -----------------------------------------------------------------------------

static Forces compute_forces(
    const std::vector<Cons>& U,
    const std::vector<CellGeom>& geom,
    const std::vector<Panel>& panels,
    const Config& c,
    double dx,
    double dy,
    std::vector<PanelForce>* pfo = nullptr
) {
    if (c.use_panel_force_for_primary_coefficients) {
        return compute_forces_panel_sampling(U, geom, panels, c, dx, dy, pfo);
    }

    return compute_forces_wall_pressure(U, geom, panels, c, dx, dy, pfo);
}

// -----------------------------------------------------------------------------
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
