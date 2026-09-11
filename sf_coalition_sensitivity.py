"""Certified sensitivity scan for objective heterogeneity and governance scale."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sf_revision_common import (
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    solve_benchmarks,
    solve_profile,
    write_json,
)


RAW_PATH = Path("sf_coalition_sensitivity_results.json")
SUMMARY_PATH = Path("sf_coalition_sensitivity_summary.json")
ALPHAS = (0.5, 0.9)
GAMMAS = (0.0, 0.25, 0.5, 0.75, 1.0)
TYPE_PAIRS = {
    "wide": (0.25, 1.75),
    "baseline": (0.50, 1.50),
    "narrow": (0.75, 1.25),
}
SLOPE_SCALES = (0.5, 1.0, 2.0)
SHARES = np.full(4, 0.25)
PARTITIONS = {
    "assortative": ((0, 1), (2, 3)),
    "mixed": ((0, 2), (1, 3)),
}


def record(alpha, gamma, type_name, types, slope_scale, partition_name,
           partition, sol, game, j_ue, j_so):
    return {
        "alpha": float(alpha),
        "gamma": float(gamma),
        "type_design": type_name,
        "objective_types": list(map(float, types)),
        "objective_gap": float(types[2] - types[0]),
        "slope_scale": float(slope_scale),
        "partition": partition_name,
        "coalitions": [list(block) for block in partition],
        "company_HHI": 0.25,
        "coalition_HHI": 0.5,
        "J": float(sol["J"]),
        "recovery": float((j_ue - sol["J"]) / (j_ue - j_so)),
        "VI_gap": float(sol["gap"]),
        "omitted_route_slack": float(sol["max_reduced_cost"]),
        "route_cost_scale": float(sol["route_cost_scale"]),
        "certified": bool(sol["certified"]),
        "iterations": int(sol["iters"]),
        "cg_rounds": int(sol["cg_rounds"]),
        "paths": int(sol["n_paths"]),
    }


def summarize(raw):
    rows = raw["profiles"]
    groups = {}
    for row in rows:
        key = (row["alpha"], row["type_design"], row["slope_scale"], row["gamma"])
        groups.setdefault(key, {})[row["partition"]] = row
    contrasts = []
    for key, values in sorted(groups.items()):
        if set(values) != {"assortative", "mixed"}:
            continue
        assortative, mixed = values["assortative"], values["mixed"]
        contrasts.append({
            "alpha": key[0],
            "type_design": key[1],
            "objective_gap": assortative["objective_gap"],
            "slope_scale": key[2],
            "gamma": key[3],
            "delta_J_assortative_minus_mixed": assortative["J"] - mixed["J"],
            "delta_J_over_J_UE": (assortative["J"] - mixed["J"]) / raw["J_UE"],
            "delta_recovery_percentage_points":
                100.0 * (assortative["recovery"] - mixed["recovery"]),
        })
    by_design = {}
    for row in contrasts:
        by_design.setdefault((row["alpha"], row["type_design"], row["slope_scale"]), []).append(row)
    design_summary = []
    for key, items in sorted(by_design.items()):
        signs = {np.sign(item["delta_recovery_percentage_points"]) for item in items
                 if abs(item["delta_recovery_percentage_points"]) > 1e-8}
        design_summary.append({
            "alpha": key[0],
            "type_design": key[1],
            "objective_gap": items[0]["objective_gap"],
            "slope_scale": key[2],
            "max_abs_recovery_contrast_pp": max(
                abs(item["delta_recovery_percentage_points"]) for item in items),
            "sign_reversal_over_gamma": len(signs) > 1,
            "min_recovery_contrast_pp": min(item["delta_recovery_percentage_points"] for item in items),
            "max_recovery_contrast_pp": max(item["delta_recovery_percentage_points"] for item in items),
        })
    return {
        "protocol": raw["protocol"],
        "rows": len(rows),
        "certified_rows": sum(bool(row["certified"]) for row in rows),
        "max_VI_gap": max(row["VI_gap"] for row in rows),
        "max_relative_omitted_route_slack": max(
            row["omitted_route_slack"] / row["route_cost_scale"] for row in rows),
        "assortative_minus_mixed": contrasts,
        "design_summary": design_summary,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game, so = solve_benchmarks(sf)
    base_slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])

    prior = {}
    if args.resume and RAW_PATH.exists():
        old = json.loads(RAW_PATH.read_text())
        prior = {(r["alpha"], r["gamma"], r["type_design"],
                  r["slope_scale"], r["partition"]): r
                 for r in old.get("profiles", [])}
    raw = {
        "protocol": {
            "governance": "identity-preserving coordination",
            "objective_type_pairs": {k: list(v) for k, v in TYPE_PAIRS.items()},
            "slope_scales": list(SLOPE_SCALES),
            "alphas": list(ALPHAS),
            "gammas": list(GAMMAS),
            "partitions": {k: [list(b) for b in v] for k, v in PARTITIONS.items()},
            "total_AV_mass_fixed": True,
            "company_masses": "equal",
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }
    total_rows = len(ALPHAS) * len(GAMMAS) * len(TYPE_PAIRS) * len(SLOPE_SCALES) * len(PARTITIONS)
    done = 0
    for alpha in ALPHAS:
        for type_name, (low, high) in TYPE_PAIRS.items():
            types = np.array([low, low, high, high], dtype=float)
            for slope_scale in SLOPE_SCALES:
                slope = base_slope * slope_scale
                warm = {}
                for gamma in GAMMAS:
                    for partition_name, partition in PARTITIONS.items():
                        done += 1
                        key = (alpha, gamma, type_name, slope_scale, partition_name)
                        if key in prior:
                            raw["profiles"].append(prior[key])
                            print(f"[{done:03d}/{total_rows}] reused {key}", flush=True)
                            continue
                        sol, game = solve_profile(
                            sf, routes, oracle, alpha, SHARES, lam_beta,
                            coalitions=partition, objective_types=types,
                            coordination_gamma=gamma,
                            management_slope=slope,
                            warm=warm.get(partition_name), verbose=args.verbose,
                        )
                        row = record(alpha, gamma, type_name, types, slope_scale,
                                     partition_name, partition, sol, game,
                                     ue["J"], so["J"])
                        raw["profiles"].append(row)
                        warm[partition_name] = (sol["fH"], sol["fF"])
                        write_json(RAW_PATH, raw)
                        print(f"[{done:03d}/{total_rows}] alpha={alpha:.1f} type={type_name:<8} "
                              f"scale={slope_scale:g} gamma={gamma:.2f} {partition_name:<10} "
                              f"rho={row['recovery']:.5f} gap={row['VI_gap']:.1e} cert={row['certified']}",
                              flush=True)
    raw["profiles"].sort(key=lambda r: (r["alpha"], r["type_design"],
                                         r["slope_scale"], r["gamma"], r["partition"]))
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if summary["certified_rows"] != summary["rows"]:
        raise SystemExit("one or more sensitivity profiles failed certification")


if __name__ == "__main__":
    main()
