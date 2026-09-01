#pragma once

#include "solver.hpp"

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
