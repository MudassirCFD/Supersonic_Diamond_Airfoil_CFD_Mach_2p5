"""Analytical references used to check the Euler solution."""

import math

from config import (
    AOA_RAD, BETA_ACK, CHORD, GAMMA, M_INF, P_INF, Q_INF,
    THETA_ACK, T_OVER_C,
)

def ackeret():
    cd = 4.0 * (THETA_ACK * THETA_ACK + AOA_RAD * AOA_RAD) / BETA_ACK
    cl = 4.0 * AOA_RAD / BETA_ACK

    cp_uf = 2.0 * (THETA_ACK - AOA_RAD) / BETA_ACK
    cp_ur = -2.0 * (THETA_ACK + AOA_RAD) / BETA_ACK
    cp_lf = 2.0 * (THETA_ACK + AOA_RAD) / BETA_ACK
    cp_lr = -2.0 * (THETA_ACK - AOA_RAD) / BETA_ACK

    return {
        "Cd": cd,
        "Cl": cl,
        "Cp_uf": cp_uf,
        "Cp_ur": cp_ur,
        "Cp_lf": cp_lf,
        "Cp_lr": cp_lr,
    }

def theta_beta_mach(beta, mach, theta, gamma=GAMMA):
    s = math.sin(beta)
    mn2 = mach * mach * s * s

    if mn2 <= 1.0:
        return 1.0e30

    return (
        math.atan(
            2.0 / math.tan(beta)
            * (mn2 - 1.0)
            / (mach * mach * (gamma + math.cos(2.0 * beta)) + 2.0)
        )
        - theta
    )

def solve_oblique_shock(mach, theta):
    if theta <= 0.0:
        return None

    mu = math.asin(1.0 / mach)
    lo = None
    hi = None

    previous_x = mu + 1.0e-8
    previous_value = theta_beta_mach(previous_x, mach, theta)

    for k in range(1, 8000):
        x = mu + 1.0e-8 + k * (math.pi / 2.0 - mu - 2.0e-8) / 8000.0
        value = theta_beta_mach(x, mach, theta)

        if previous_value * value <= 0.0:
            lo = previous_x
            hi = x
            break

        previous_x = x
        previous_value = value

    if lo is None:
        return None

    for _ in range(240):
        mid = 0.5 * (lo + hi)
        if theta_beta_mach(lo, mach, theta) * theta_beta_mach(mid, mach, theta) <= 0.0:
            hi = mid
        else:
            lo = mid

        if hi - lo < 1.0e-13:
            break

    beta = 0.5 * (lo + hi)
    mn1_sq = mach * mach * math.sin(beta) ** 2

    p2p1 = 1.0 + 2.0 * GAMMA * (mn1_sq - 1.0) / (GAMMA + 1.0)
    rho2rho1 = (GAMMA + 1.0) * mn1_sq / ((GAMMA - 1.0) * mn1_sq + 2.0)

    mn2_sq = ((GAMMA - 1.0) * mn1_sq + 2.0) / (
        2.0 * GAMMA * mn1_sq - (GAMMA - 1.0)
    )

    mach2 = math.sqrt(mn2_sq) / math.sin(beta - theta)
    cp = 2.0 * (p2p1 - 1.0) / (GAMMA * mach * mach)

    return {
        "beta": beta,
        "beta_deg": math.degrees(beta),
        "theta_deg": math.degrees(theta),
        "M2": mach2,
        "p2p1": p2p1,
        "rho2rho1": rho2rho1,
        "Cp": cp,
    }

def nu_pm(mach, gamma=GAMMA):
    mach = max(float(mach), 1.000001)
    return (
        math.sqrt((gamma + 1.0) / (gamma - 1.0))
        * math.atan(math.sqrt((gamma - 1.0) / (gamma + 1.0) * (mach * mach - 1.0)))
        - math.atan(math.sqrt(mach * mach - 1.0))
    )

def inv_nu_pm(target, guess=2.5, gamma=GAMMA):
    mach = max(1.01, float(guess))

    for _ in range(300):
        nu = nu_pm(mach, gamma)
        derivative = math.sqrt(mach * mach - 1.0) / (
            mach * (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
        )

        delta = (target - nu) / max(derivative, 1.0e-14)
        mach += delta

        if mach < 1.001:
            mach = 1.001

        if abs(delta) < 1.0e-12:
            break

    return mach

def prandtl_meyer_expansion(mach1, turn_rad):
    nu1 = nu_pm(mach1)
    nu2 = nu1 + abs(turn_rad)
    mach2 = inv_nu_pm(nu2, mach1 + 0.3)

    r1 = 1.0 + 0.5 * (GAMMA - 1.0) * mach1 * mach1
    r2 = 1.0 + 0.5 * (GAMMA - 1.0) * mach2 * mach2
    p2p1 = (r1 / r2) ** (GAMMA / (GAMMA - 1.0))

    return {
        "M1": mach1,
        "M2": mach2,
        "nu1_deg": math.degrees(nu1),
        "nu2_deg": math.degrees(nu2),
        "turn_deg": math.degrees(abs(turn_rad)),
        "p2p1": p2p1,
    }

def nonlinear_theory():
    theta_exact = math.atan(T_OVER_C)

    deflection_upper_front = theta_exact - AOA_RAD
    deflection_lower_front = theta_exact + AOA_RAD

    shock_upper_front = solve_oblique_shock(M_INF, deflection_upper_front)
    shock_lower_front = solve_oblique_shock(M_INF, deflection_lower_front)

    mach_upper_front = shock_upper_front["M2"] if shock_upper_front else M_INF
    mach_lower_front = shock_lower_front["M2"] if shock_lower_front else M_INF

    cp_upper_front = shock_upper_front["Cp"] if shock_upper_front else 0.0
    cp_lower_front = shock_lower_front["Cp"] if shock_lower_front else 0.0

    # Each shoulder changes the surface direction by twice the diamond half-angle.
    turn_upper_rear = 2.0 * theta_exact
    turn_lower_rear = 2.0 * theta_exact

    pm_upper_rear = prandtl_meyer_expansion(mach_upper_front, turn_upper_rear)
    pm_lower_rear = prandtl_meyer_expansion(mach_lower_front, turn_lower_rear)

    pressure_upper_rear = (
        (shock_upper_front["p2p1"] if shock_upper_front else 1.0)
        * P_INF
        * pm_upper_rear["p2p1"]
    )

    pressure_lower_rear = (
        (shock_lower_front["p2p1"] if shock_lower_front else 1.0)
        * P_INF
        * pm_lower_rear["p2p1"]
    )

    cp_upper_rear = (pressure_upper_rear - P_INF) / Q_INF
    cp_lower_rear = (pressure_lower_rear - P_INF) / Q_INF

    half_chord = 0.5 * CHORD
    half_thickness = 0.5 * T_OVER_C * CHORD
    panel_length = math.sqrt(half_chord * half_chord + half_thickness * half_thickness)

    sin_theta = half_thickness / panel_length
    cos_theta = half_chord / panel_length

    normals = [
        [-sin_theta, +cos_theta],  # upper front
        [+sin_theta, +cos_theta],  # upper rear
        [+sin_theta, -cos_theta],  # lower rear
        [-sin_theta, -cos_theta],  # lower front
    ]

    cp_faces = [cp_upper_front, cp_upper_rear, cp_lower_rear, cp_lower_front]

    fx = 0.0
    fy = 0.0

    for normal, cp_value in zip(normals, cp_faces):
        pressure = P_INF + cp_value * Q_INF
        fx += -pressure * normal[0] * panel_length
        fy += -pressure * normal[1] * panel_length

    drag = fx * math.cos(AOA_RAD) + fy * math.sin(AOA_RAD)
    lift = -fx * math.sin(AOA_RAD) + fy * math.cos(AOA_RAD)

    return {
        "Cd": drag / (Q_INF * CHORD),
        "Cl": lift / (Q_INF * CHORD),
        "Cp_uf": cp_upper_front,
        "Cp_ur": cp_upper_rear,
        "Cp_lf": cp_lower_front,
        "Cp_lr": cp_lower_rear,
        "os_uf": shock_upper_front,
        "os_lf": shock_lower_front,
        "pm_ur": pm_upper_rear,
        "pm_lr": pm_lower_rear,
        "turn_ur_deg": math.degrees(turn_upper_rear),
        "turn_lr_deg": math.degrees(turn_lower_rear),
        "theta_exact_deg": math.degrees(theta_exact),
        "theta_exact_rad": theta_exact,
    }

