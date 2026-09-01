#pragma once

#include "geometry.hpp"

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
