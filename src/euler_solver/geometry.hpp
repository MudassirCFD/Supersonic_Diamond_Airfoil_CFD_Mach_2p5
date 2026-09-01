#pragma once

#include "flow.hpp"

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
