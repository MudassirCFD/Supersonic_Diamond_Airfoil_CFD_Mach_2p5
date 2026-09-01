#pragma once

#include "types.hpp"

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
