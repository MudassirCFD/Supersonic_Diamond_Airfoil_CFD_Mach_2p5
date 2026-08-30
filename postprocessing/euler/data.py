"""Read solver CSV files and prepare data used by the plots."""

import glob
import os

import numpy as np
import pandas as pd

from config import (
    CHORD, GAS_R, GAMMA, M_INF, P0_INF, P_INF, Q_INF, RHO_INF,
    T_INF, X_LE, numeric_key,
)
from theory import ackeret, nonlinear_theory

def find_existing_file(folder, names):
    for name in names:
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            return path
    return None

def detect_result_files(folder):
    candidates = [
        {
            "name": "weno5",
            "solution": ["weno5_solution.csv", "weno5_solution"],
            "residual": ["weno5_residual.csv", "weno5_residual"],
            "forces": ["weno5_forces.csv", "weno5_forces"],
            "force_history": ["weno5_force_history.csv", "weno5_force_history"],
        },
        {
            "name": "hllc",
            "solution": ["hllc_solution.csv", "hllc_solution"],
            "residual": ["hllc_residual.csv", "hllc_residual"],
            "forces": ["hllc_forces.csv", "hllc_forces"],
            "force_history": ["hllc_force_history.csv", "hllc_force_history"],
        },
        {
            "name": "default",
            "solution": ["solution.csv"],
            "residual": ["residual.csv"],
            "forces": ["forces.csv"],
            "force_history": ["force_history.csv"],
        },
    ]

    for case in candidates:
        paths = {}
        ok = True

        for key in ["solution", "residual", "forces", "force_history"]:
            path = find_existing_file(folder, case[key])
            if path is None:
                ok = False
                break
            paths[key] = path

        if ok:
            print("detected case:", case["name"])
            for key, value in paths.items():
                print(f"{key:14s}: {value}")
            return paths

    raise FileNotFoundError(
        "Could not find a complete CSV set. Expected files such as "
        "weno5_solution.csv, weno5_residual.csv, weno5_forces.csv, "
        "and weno5_force_history.csv."
    )

def find_snapshot_files(folder, stride=1):
    patterns = [
        "weno5_snap_*.csv",
        "weno5_snap_*",
        "*snap_*.csv",
        "*snap_*",
    ]

    files = []

    for pattern in patterns:
        found = sorted(glob.glob(os.path.join(folder, pattern)), key=numeric_key)
        found = [
            file for file in found
            if os.path.isfile(file)
            and file.lower().endswith(".csv")
        ]

        if found:
            files = found
            break

    if stride > 1:
        files = files[::stride]

    print("snapshot files found:", len(files))
    if files:
        print("first snapshot:", os.path.basename(files[0]))
        print("last snapshot :", os.path.basename(files[-1]))

    return files

def grid_from_dataframe(df, fields):
    if "is_solid" not in df.columns:
        raise ValueError("CSV file must contain the column 'is_solid'.")

    fluid = df[df["is_solid"] == 0].copy()

    xs = np.sort(df["x"].unique())
    ys = np.sort(df["y"].unique())

    nx = len(xs)
    ny = len(ys)

    if nx < 2 or ny < 2:
        raise ValueError("Grid appears too small or empty.")

    x_values = fluid["x"].to_numpy()
    y_values = fluid["y"].to_numpy()

    ii = np.searchsorted(xs, x_values)
    jj = np.searchsorted(ys, y_values)

    grids = {}

    solid_grid = np.zeros((ny, nx), dtype=bool)
    solid_rows = df[df["is_solid"] == 1]
    if not solid_rows.empty:
        si = np.searchsorted(xs, solid_rows["x"].to_numpy())
        sj = np.searchsorted(ys, solid_rows["y"].to_numpy())
        solid_grid[sj, si] = True

    grids["solid_mask"] = solid_grid

    for field in fields:
        grid = np.full((ny, nx), np.nan)
        if field in fluid.columns:
            grid[jj, ii] = fluid[field].to_numpy()
        grids[field] = grid

    return xs, ys, grids

def load_solution(path):
    print("loading:", path)

    df = pd.read_csv(path)

    required = ["x", "y", "is_solid", "rho", "u", "v", "p", "Cp", "T", "a", "Mach"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"solution file missing required columns: {missing}")

    fields = ["rho", "u", "v", "p", "Cp", "T", "a", "Mach"]
    xs, ys, grids = grid_from_dataframe(df, fields)

    grids["V"] = np.sqrt(grids["u"] ** 2 + grids["v"] ** 2)
    grids["q"] = 0.5 * grids["rho"] * grids["V"] ** 2
    grids["q_over_qinf"] = grids["q"] / Q_INF

    grids["T0"] = grids["T"] * (1.0 + 0.5 * (GAMMA - 1.0) * grids["Mach"] ** 2)
    grids["p0_isentropic"] = grids["p"] * (1.0 + 0.5 * (GAMMA - 1.0) * grids["Mach"] ** 2) ** (GAMMA / (GAMMA - 1.0))
    grids["p0_ratio"] = grids["p0_isentropic"] / P0_INF

    grids["entropy_proxy"] = np.log(np.maximum(grids["p"], 1.0e-30) / np.maximum(grids["rho"], 1.0e-30) ** GAMMA)

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    pressure = np.nan_to_num(grids["p"], nan=P_INF)
    rho = np.nan_to_num(grids["rho"], nan=RHO_INF)

    dpdx = np.gradient(pressure, dx, axis=1)
    dpdy = np.gradient(pressure, dy, axis=0)

    drdx = np.gradient(rho, dx, axis=1)
    drdy = np.gradient(rho, dy, axis=0)

    grids["gradp"] = np.sqrt(dpdx ** 2 + dpdy ** 2)
    grids["gradrho"] = np.sqrt(drdx ** 2 + drdy ** 2)

    grids["gradp"][np.isnan(grids["p"])] = np.nan
    grids["gradrho"][np.isnan(grids["rho"])] = np.nan

    print(f"grid: {len(xs)} x {len(ys)}")
    print(f"fluid cells: {np.count_nonzero(~np.isnan(grids['Mach']))}")

    return xs, ys, grids, df

def load_snapshot_grid(path, field_name):
    df = pd.read_csv(path)
    xs, ys, grids = grid_from_dataframe(df, [field_name])
    return xs, ys, grids[field_name]

def load_snapshot_schlieren(path):
    df = pd.read_csv(path)
    xs, ys, grids = grid_from_dataframe(df, ["rho", "p"])

    if "rho" in grids and not np.all(np.isnan(grids["rho"])):
        base = np.nan_to_num(grids["rho"], nan=RHO_INF)
        label = "rho"
    else:
        base = np.nan_to_num(grids["p"], nan=P_INF)
        label = "p"

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    dfdx = np.gradient(base, dx, axis=1)
    dfdy = np.gradient(base, dy, axis=0)
    grad = np.sqrt(dfdx ** 2 + dfdy ** 2)

    upper = np.nanpercentile(grad, 99.8)
    grad = np.clip(grad, 1.0, upper)
    schlieren = np.log10(grad)

    return xs, ys, schlieren, label

def sample_line_at_y(xs, ys, field, y_target):
    j = int(np.argmin(np.abs(ys - y_target)))
    return field[j, :], ys[j]

def extract_surface_cp(df):
    fluid = df[df["is_solid"] == 0].copy()
    solid = df[df["is_solid"] == 1].copy()

    if solid.empty:
        return pd.DataFrame(columns=["x_over_c", "Cp_upper", "Cp_lower"])

    xs = np.sort(solid["x"].unique())

    rows = []

    for x in xs:
        if x < X_LE or x > X_LE + CHORD:
            continue

        fluid_column = fluid[np.abs(fluid["x"] - x) < 1.0e-10]
        solid_column = solid[np.abs(solid["x"] - x) < 1.0e-10]

        if fluid_column.empty or solid_column.empty:
            continue

        y_top = solid_column["y"].max()
        y_bottom = solid_column["y"].min()

        upper = fluid_column[fluid_column["y"] > y_top].sort_values("y")
        lower = fluid_column[fluid_column["y"] < y_bottom].sort_values("y", ascending=False)

        if upper.empty or lower.empty:
            continue

        rows.append({
            "x_over_c": (x - X_LE) / CHORD,
            "x": x,
            "y_top_solid": y_top,
            "y_bottom_solid": y_bottom,
            "Cp_upper": float(upper.iloc[0]["Cp"]),
            "Cp_lower": float(lower.iloc[0]["Cp"]),
            "p_upper": float(upper.iloc[0]["p"]),
            "p_lower": float(lower.iloc[0]["p"]),
            "Mach_upper": float(upper.iloc[0]["Mach"]),
            "Mach_lower": float(lower.iloc[0]["Mach"]),
        })

    return pd.DataFrame(rows)

def panel_averages(surface_df):
    if surface_df.empty:
        return pd.DataFrame()

    regions = {
        "upper_front": (surface_df["x_over_c"] >= 0.06) & (surface_df["x_over_c"] < 0.48),
        "upper_rear": (surface_df["x_over_c"] > 0.52) & (surface_df["x_over_c"] <= 0.96),
        "lower_front": (surface_df["x_over_c"] >= 0.06) & (surface_df["x_over_c"] < 0.48),
        "lower_rear": (surface_df["x_over_c"] > 0.52) & (surface_df["x_over_c"] <= 0.96),
    }

    rows = []
    for name, mask in regions.items():
        if name.startswith("upper"):
            col = "Cp_upper"
        else:
            col = "Cp_lower"

        values = surface_df.loc[mask, col].to_numpy()
        if values.size == 0:
            mean = np.nan
            std = np.nan
        else:
            mean = float(np.nanmean(values))
            std = float(np.nanstd(values))

        rows.append({
            "panel": name,
            "Cp_mean_CFD": mean,
            "Cp_std_CFD": std,
            "samples": int(values.size),
        })

    return pd.DataFrame(rows)

def extract_field_statistics(grids):
    rows = []
    keys = [
        "rho", "u", "v", "p", "Cp", "T", "a", "Mach", "V", "q_over_qinf",
        "T0", "p0_ratio", "entropy_proxy", "gradp", "gradrho",
    ]

    for key in keys:
        if key not in grids:
            continue

        data = grids[key]
        clean = data[np.isfinite(data)]

        if clean.size == 0:
            continue

        rows.append({
            "field": key,
            "min": float(np.nanmin(clean)),
            "p01": float(np.nanpercentile(clean, 1)),
            "mean": float(np.nanmean(clean)),
            "p99": float(np.nanpercentile(clean, 99)),
            "max": float(np.nanmax(clean)),
        })

    return pd.DataFrame(rows)

def extract_centerline_tables(xs, ys, grids):
    cuts = [
        ("centre", 0.0),
        ("upper_near_body", +0.08 * CHORD),
        ("lower_near_body", -0.08 * CHORD),
        ("upper_far", +0.25 * CHORD),
        ("lower_far", -0.25 * CHORD),
    ]

    rows = []

    for name, y_target in cuts:
        _, y_actual = sample_line_at_y(xs, ys, grids["Mach"], y_target)
        for i, x in enumerate(xs):
            rows.append({
                "cut": name,
                "x_over_c": (x - X_LE) / CHORD,
                "x": x,
                "y": y_actual,
                "Mach": grids["Mach"][int(np.argmin(np.abs(ys - y_actual))), i],
                "p": grids["p"][int(np.argmin(np.abs(ys - y_actual))), i],
                "Cp": grids["Cp"][int(np.argmin(np.abs(ys - y_actual))), i],
                "rho": grids["rho"][int(np.argmin(np.abs(ys - y_actual))), i],
                "T": grids["T"][int(np.argmin(np.abs(ys - y_actual))), i],
                "V": grids["V"][int(np.argmin(np.abs(ys - y_actual))), i],
            })

    return pd.DataFrame(rows)

def lower_shock_profile_table(xs, ys, grids, y_target=-0.35 * CHORD):
    """Extract a horizontal cut through the lower leading-edge oblique shock."""
    j = int(np.argmin(np.abs(ys - y_target)))
    y_actual = float(ys[j])
    nonlinear = nonlinear_theory()
    shock = nonlinear["os_lf"]

    p_rh = P_INF * shock["p2p1"] if shock else np.nan
    rho_rh = RHO_INF * shock["rho2rho1"] if shock else np.nan
    t_rh = p_rh / (rho_rh * GAS_R) if shock else np.nan
    mach_rh = shock["M2"] if shock else np.nan

    rows = []
    for i, x in enumerate(xs):
        rows.append({
            "x_over_c": (float(x) - X_LE) / CHORD,
            "x": float(x),
            "y": y_actual,
            "Mach_CFD": float(grids["Mach"][j, i]),
            "pressure_CFD_Pa": float(grids["p"][j, i]),
            "density_CFD_kg_m3": float(grids["rho"][j, i]),
            "temperature_CFD_K": float(grids["T"][j, i]),
            "Mach_freestream": M_INF,
            "pressure_freestream_Pa": P_INF,
            "density_freestream_kg_m3": RHO_INF,
            "temperature_freestream_K": T_INF,
            "Mach_RH_lower_shock": mach_rh,
            "pressure_RH_lower_shock_Pa": p_rh,
            "density_RH_lower_shock_kg_m3": rho_rh,
            "temperature_RH_lower_shock_K": t_rh,
        })

    return pd.DataFrame(rows)

def write_data_extracts(xs, ys, grids, df, force_row, residual_df, force_history_df, folders):
    ack = ackeret()
    nonlinear = nonlinear_theory()

    field_stats = extract_field_statistics(grids)
    field_stats.to_csv(os.path.join(folders.data, "field_statistics.csv"), index=False)

    surface = extract_surface_cp(df)
    surface.to_csv(os.path.join(folders.data, "surface_cp_extract.csv"), index=False)

    panel_avg = panel_averages(surface)
    if not panel_avg.empty:
        theory_map = {
            "upper_front": (ack["Cp_uf"], nonlinear["Cp_uf"]),
            "upper_rear": (ack["Cp_ur"], nonlinear["Cp_ur"]),
            "lower_front": (ack["Cp_lf"], nonlinear["Cp_lf"]),
            "lower_rear": (ack["Cp_lr"], nonlinear["Cp_lr"]),
        }
        panel_avg["Cp_Ackeret"] = panel_avg["panel"].map(lambda s: theory_map[s][0])
        panel_avg["Cp_nonlinear_OS_PM"] = panel_avg["panel"].map(lambda s: theory_map[s][1])
        panel_avg["error_vs_nonlinear_percent"] = (
            100.0 * (panel_avg["Cp_mean_CFD"] - panel_avg["Cp_nonlinear_OS_PM"])
            / panel_avg["Cp_nonlinear_OS_PM"].replace(0.0, np.nan)
        )
    panel_avg.to_csv(os.path.join(folders.data, "surface_cp_panel_averages.csv"), index=False)

    centerline = extract_centerline_tables(xs, ys, grids)
    centerline.to_csv(os.path.join(folders.data, "axial_cuts_extract.csv"), index=False)

    shock_profile = lower_shock_profile_table(xs, ys, grids)
    shock_profile.to_csv(os.path.join(folders.data, "lower_oblique_shock_profile_extract.csv"), index=False)

    theory_rows = [
        {"source": "Ackeret", "Cd": ack["Cd"], "Cl": ack["Cl"],
         "Cp_upper_front": ack["Cp_uf"], "Cp_upper_rear": ack["Cp_ur"],
         "Cp_lower_front": ack["Cp_lf"], "Cp_lower_rear": ack["Cp_lr"]},
        {"source": "Nonlinear OS plus PM", "Cd": nonlinear["Cd"], "Cl": nonlinear["Cl"],
         "Cp_upper_front": nonlinear["Cp_uf"], "Cp_upper_rear": nonlinear["Cp_ur"],
         "Cp_lower_front": nonlinear["Cp_lf"], "Cp_lower_rear": nonlinear["Cp_lr"]},
        {"source": "CFD", "Cd": float(force_row["Cd"]), "Cl": float(force_row["Cl"]),
         "Cp_upper_front": np.nan, "Cp_upper_rear": np.nan,
         "Cp_lower_front": np.nan, "Cp_lower_rear": np.nan},
    ]
    theory_df = pd.DataFrame(theory_rows)
    theory_df.to_csv(os.path.join(folders.data, "theory_and_cfd_force_summary.csv"), index=False)

    shock_rows = []
    for label, key in [("upper_front_shock", "os_uf"), ("lower_front_shock", "os_lf")]:
        shock = nonlinear[key]
        if shock:
            shock_rows.append({
                "state": label,
                "process": "oblique shock",
                "theta_deg": shock["theta_deg"],
                "beta_deg": shock["beta_deg"],
                "M_downstream": shock["M2"],
                "p2p1": shock["p2p1"],
                "rho2rho1": shock["rho2rho1"],
                "Cp": shock["Cp"],
            })

    for label, key in [("upper_rear_expansion", "pm_ur"), ("lower_rear_expansion", "pm_lr")]:
        pm = nonlinear[key]
        shock_rows.append({
            "state": label,
            "process": "Prandtl Meyer expansion",
            "theta_deg": pm["turn_deg"],
            "beta_deg": np.nan,
            "M_downstream": pm["M2"],
            "p2p1": pm["p2p1"],
            "rho2rho1": np.nan,
            "Cp": np.nan,
        })

    pd.DataFrame(shock_rows).to_csv(os.path.join(folders.data, "shock_and_prandtl_meyer_states.csv"), index=False)

    residual_summary = {
        "iterations": int(residual_df["iteration"].iloc[-1]),
        "initial_residual": float(residual_df["residual"].iloc[0]) if "residual" in residual_df else np.nan,
        "final_residual": float(residual_df["residual"].iloc[-1]) if "residual" in residual_df else np.nan,
        "final_normalized_residual": float(residual_df["normalized_residual"].iloc[-1]) if "normalized_residual" in residual_df else np.nan,
        "final_solution_change": float(residual_df["solution_change"].iloc[-1]) if "solution_change" in residual_df else np.nan,
        "final_Cd_history": float(force_history_df["Cd"].iloc[-1]) if "Cd" in force_history_df else np.nan,
        "final_Cl_history": float(force_history_df["Cl"].iloc[-1]) if "Cl" in force_history_df else np.nan,
        "final_Cd_force_file": float(force_row["Cd"]),
        "final_Cl_force_file": float(force_row["Cl"]),
        "Cd_error_vs_Ackeret_percent": 100.0 * (float(force_row["Cd"]) - ack["Cd"]) / ack["Cd"],
        "Cl_error_vs_Ackeret_percent": 100.0 * (float(force_row["Cl"]) - ack["Cl"]) / ack["Cl"],
        "Cd_error_vs_nonlinear_percent": 100.0 * (float(force_row["Cd"]) - nonlinear["Cd"]) / nonlinear["Cd"],
        "Cl_error_vs_nonlinear_percent": 100.0 * (float(force_row["Cl"]) - nonlinear["Cl"]) / nonlinear["Cl"],
    }
    pd.DataFrame([residual_summary]).to_csv(os.path.join(folders.data, "run_summary.csv"), index=False)

    print("saved data extracts in:", folders.data)

