"""Certified n=4 Sioux Falls profiles with company accounting costs.

This runner expands the existing partition scan along the objective-
heterogeneity dimension needed by the coalition-formation analysis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coalition_private_costs import company_accounting_costs
from sf_coalition_partition_scan import PARTITIONS, PARTITION_IDS, coalition_hhi
from sf_revision_common import (
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    solve_benchmarks,
    solve_profile,
    write_json,
)


OUT = Path("sf_coalition_formation_profiles.json")
ALPHAS = (0.5, 0.9)
GAMMAS = (0.5, 1.0)
TYPE_DESIGNS = {
    "homogeneous": np.array([1.0, 1.0, 1.0, 1.0]),
    "baseline": np.array([0.5, 0.5, 1.5, 1.5]),
    "wide": np.array([0.25, 0.25, 1.75, 1.75]),
}
SHARES = np.full(4, 0.25)


def make_row(alpha, gamma, design, types, partition, sol, game, j_ue, j_so):
    return {
        "network": "sioux_falls",
        "alpha": float(alpha),
        "gamma": float(gamma),
        "type_design": design,
        "objective_types": types.tolist(),
        "objective_gap": float(np.max(types) - np.min(types)),
        "partition_id": PARTITION_IDS[partition],
        "partition": [list(block) for block in partition],
        "n_companies": 4,
        "n_coalitions": len(partition),
        "coalition_HHI": coalition_hhi(partition),
        "J": float(sol["J"]),
        "recovery": float((j_ue - sol["J"]) / (j_ue - j_so)),
        "VI_gap": float(sol["gap"]),
        "omitted_route_slack": float(sol["max_reduced_cost"]),
        "route_cost_scale": float(sol["route_cost_scale"]),
        "certified": bool(sol["certified"]),
        "iterations": int(sol["iters"]),
        "cg_rounds": int(sol["cg_rounds"]),
        "paths": int(sol["n_paths"]),
        **company_accounting_costs(sol, game),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total_demand = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total_demand)
    ue, ue_game, so = solve_benchmarks(sf)
    slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])

    prior = {}
    if args.resume and OUT.exists():
        old = json.loads(OUT.read_text())
        prior = {
            (r["alpha"], r["gamma"], r["type_design"], r["partition_id"]): r
            for r in old.get("profiles", [])
            if "company_accounting_cost" in r
        }

    raw = {
        "protocol": {
            "network": "Sioux Falls",
            "governance": "identity-preserving coordination",
            "alphas": list(ALPHAS),
            "gammas": list(GAMMAS),
            "type_designs": {k: v.tolist() for k, v in TYPE_DESIGNS.items()},
            "partitions": len(PARTITIONS),
            "company_cost": (
                "sum_e f_i,e*c_e(x) + 0.5*sum_e management_slope_e*"
                "tau_i*f_i,e^2; coalition cross terms are not allocated"
            ),
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
        },
        "benchmarks": {"J_UE": float(ue["J"]), "J_SO": float(so["J"])},
        "profiles": [],
    }
    total = len(ALPHAS) * len(GAMMAS) * len(TYPE_DESIGNS) * len(PARTITIONS)
    done = 0
    for design, types in TYPE_DESIGNS.items():
        warm = {}
        for alpha in ALPHAS:
            for gamma in GAMMAS:
                for partition in PARTITIONS:
                    done += 1
                    pid = PARTITION_IDS[partition]
                    key = (alpha, gamma, design, pid)
                    if key in prior:
                        raw["profiles"].append(prior[key])
                        print(f"[{done:03d}/{total}] reused {key}", flush=True)
                        continue
                    sol, game = solve_profile(
                        sf,
                        routes,
                        oracle,
                        alpha,
                        SHARES,
                        lam_beta,
                        coalitions=partition,
                        objective_types=types,
                        coordination_gamma=gamma,
                        management_slope=slope,
                        warm=warm.get((alpha, pid)),
                        verbose=args.verbose,
                    )
                    row = make_row(
                        alpha, gamma, design, types, partition, sol, game,
                        ue["J"], so["J"]
                    )
                    raw["profiles"].append(row)
                    warm[(alpha, pid)] = (sol["fH"], sol["fF"])
                    write_json(OUT, raw)
                    print(
                        f"[{done:03d}/{total}] {design:<11} alpha={alpha:.1f} "
                        f"gamma={gamma:.1f} {pid} K={len(partition)} "
                        f"gap={row['VI_gap']:.2e} cert={row['certified']}",
                        flush=True,
                    )
    raw["profiles"].sort(
        key=lambda r: (r["type_design"], r["alpha"], r["gamma"], r["partition_id"])
    )
    write_json(OUT, raw)
    certified = sum(r["certified"] for r in raw["profiles"])
    if certified != len(raw["profiles"]):
        raise SystemExit(f"formation profiles not fully certified: {certified}/{len(raw['profiles'])}")
    print(f"certified {certified}/{len(raw['profiles'])} profiles", flush=True)


if __name__ == "__main__":
    main()
