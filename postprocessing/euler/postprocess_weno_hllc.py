#!/usr/bin/env python3
"""Post-process the WENO5-HLLC reference calculation."""

import argparse
import os

import pandas as pd

from animations import make_animations
from config import make_output_folders
from convergence_plots import plot_cl_vs_cd, plot_force_history, plot_residual
from data import detect_result_files, load_solution, write_data_extracts
from field_plots import (
    plot_all_contours,
    plot_schlieren,
    plot_schlieren_shock_angles,
    plot_streamlines,
)
from verification_plots import (
    plot_force_comparison,
    plot_lower_oblique_shock_profiles,
    plot_near_body_cuts,
    plot_prandtl_meyer,
    plot_riemann_fan_comparison,
    plot_surface_cp,
    plot_theory_cp_bar,
)


def main():
    parser = argparse.ArgumentParser(
        description="Post-process WENO5-HLLC diamond-airfoil results."
    )
    parser.add_argument(
        "--dir", default=".",
        help="Folder containing the solver CSV files."
    )
    parser.add_argument(
        "--out", default=None,
        help="Output folder. Default: <input>/post"
    )
    parser.add_argument(
        "--animate", action="store_true",
        help="Create GIFs when snapshot CSV files are available."
    )
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--stride", type=int, default=1)
    args = parser.parse_args()

    input_folder = os.path.abspath(args.dir)
    output_folder = (
        os.path.abspath(args.out)
        if args.out
        else os.path.join(input_folder, "post")
    )

    folders = make_output_folders(output_folder)
    paths = detect_result_files(input_folder)

    xs, ys, grids, solution = load_solution(paths["solution"])
    residual = pd.read_csv(paths["residual"])
    force_row = pd.read_csv(paths["forces"]).iloc[0]
    force_history = pd.read_csv(paths["force_history"])

    print("Writing data extracts")
    write_data_extracts(
        xs, ys, grids, solution, force_row, residual, force_history, folders
    )

    print("Generating flow-field plots")
    plot_all_contours(xs, ys, grids, folders)
    plot_schlieren(xs, ys, grids, folders)
    plot_schlieren_shock_angles(xs, ys, grids, folders)
    plot_streamlines(xs, ys, grids, folders)

    print("Generating verification plots")
    plot_near_body_cuts(xs, ys, grids, folders)
    plot_lower_oblique_shock_profiles(xs, ys, grids, folders)
    plot_surface_cp(solution, folders)
    plot_force_comparison(force_row, folders)
    plot_prandtl_meyer(folders)
    plot_riemann_fan_comparison(folders)
    plot_theory_cp_bar(folders)

    print("Generating convergence plots")
    plot_residual(residual, folders)
    plot_force_history(force_history, force_row, folders)
    plot_cl_vs_cd(force_history, force_row, folders)

    if args.animate:
        print("Generating animations")
        make_animations(
            input_folder,
            folders,
            force_history,
            force_row,
            fps=args.fps,
            stride=max(1, args.stride),
        )

    print("Done")
    print("Output:", output_folder)


if __name__ == "__main__":
    main()
