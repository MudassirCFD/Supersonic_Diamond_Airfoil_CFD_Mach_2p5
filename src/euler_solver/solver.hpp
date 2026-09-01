#pragma once

#include "fluxes.hpp"
#include "immersed_boundary.hpp"

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
