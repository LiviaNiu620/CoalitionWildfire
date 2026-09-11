"""Certified identity-preserving versus operational-pooling comparison."""
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


RAW_PATH = Path("sf_coalition_pooling_results.json")
SUMMARY_PATH = Path("sf_coalition_pooling_summary.json")
ALPHAS = (0.5, 0.9)
SHARES = np.full(4, 0.25)
OBJECTIVE_TYPES = np.array([0.5, 0.5, 1.5, 1.5])
PARTITION = ((0, 1), (2, 3))
GAMMA = 1.0


def record(alpha, regime, sol, game, j_ue, j_so):
    return {
        "alpha": float(alpha),
        "regime": regime,
        "operational_pooling": regime == "pooling",
        "gamma": GAMMA,
        "partition": [list(block) for block in PARTITION],
        "objective_types": OBJECTIVE_TYPES.tolist(),
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
    contrasts = []
    for alpha in ALPHAS:
        identity = next(r for r in rows if r["alpha"] == alpha and r["regime"] == "identity")
        pooling = next(r for r in rows if r["alpha"] == alpha and r["regime"] == "pooling")
        contrasts.append({
            "alpha": alpha,
            "delta_J_pooling_minus_identity": pooling["J"] - identity["J"],
            "delta_J_over_J_UE": (pooling["J"] - identity["J"]) / raw["J_UE"],
            "delta_recovery_percentage_points":
                100.0 * (pooling["recovery"] - identity["recovery"]),
        })
    return {
        "protocol": raw["protocol"],
        "rows": len(rows),
        "certified_rows": sum(bool(r["certified"]) for r in rows),
        "max_VI_gap": max(r["VI_gap"] for r in rows),
        "max_relative_omitted_route_slack": max(
            r["omitted_route_slack"] / r["route_cost_scale"] for r in rows),
        "pooling_minus_identity": contrasts,
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
    slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])
    prior = {}
    if args.resume and RAW_PATH.exists():
        old = json.loads(RAW_PATH.read_text())
        prior = {(r["alpha"], r["regime"]): r for r in old.get("profiles", [])}
    raw = {
        "protocol": {
            "governance": "same identity-preserving coalition objective versus operational pooling",
            "partition": [list(block) for block in PARTITION],
            "objective_types": OBJECTIVE_TYPES.tolist(),
            "gamma": GAMMA,
            "alphas": list(ALPHAS),
            "company_masses_fixed": True,
            "aggregate_OD_fixed": True,
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }
    total_rows = len(ALPHAS) * 2
    done = 0
    for alpha in ALPHAS:
        warm = None
        for regime, pooling in (("identity", False), ("pooling", True)):
            done += 1
            key = (alpha, regime)
            if key in prior:
                raw["profiles"].append(prior[key])
                print(f"[{done:02d}/{total_rows}] reused {key}", flush=True)
                continue
            sol, game = solve_profile(
                sf, routes, oracle, alpha, SHARES, lam_beta,
                coalitions=PARTITION, objective_types=OBJECTIVE_TYPES,
                coordination_gamma=GAMMA, management_slope=slope,
                operational_pooling=pooling, warm=warm, verbose=args.verbose,
            )
            row = record(alpha, regime, sol, game, ue["J"], so["J"])
            raw["profiles"].append(row)
            warm = (sol["fH"], sol["fF"])
            write_json(RAW_PATH, raw)
            print(f"[{done:02d}/{total_rows}] alpha={alpha:.1f} {regime:<8} "
                  f"rho={row['recovery']:.5f} gap={row['VI_gap']:.1e} cert={row['certified']}",
                  flush=True)
    raw["profiles"].sort(key=lambda r: (r["alpha"], r["regime"]))
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if summary["certified_rows"] != summary["rows"]:
        raise SystemExit("one or more pooling profiles failed certification")


if __name__ == "__main__":
    main()
