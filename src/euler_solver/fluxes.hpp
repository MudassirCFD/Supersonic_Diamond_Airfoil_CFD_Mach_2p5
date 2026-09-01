#pragma once

#include "flow.hpp"

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
