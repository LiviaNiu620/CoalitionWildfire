"""Pilot n=4 wildfire snapshot partition audit on reduced Sioux Falls OD data.

This writes planning/pilot evidence only. It is not the manuscript's full
wildfire result until the all-OD scan and its validator are complete.
"""
from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np

from coalition_private_costs import company_accounting_costs
from sf_coalition_partition_scan import PARTITIONS, PARTITION_IDS, coalition_hhi
from sf_revision_common import LAM_TIMES_M, SE, network_oracle, solve_benchmarks, solve_profile
from wildfire_scenarios import apply_snapshot, build_snapshot, hazard_oracle


OUT = Path("wildfire_sioux_falls_pilot_v2_results.json")
SHARES = np.full(4, 0.25)
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
ALPHAS = (0.5, 0.9)
GAMMAS = (0.5, 1.0)


def hazard_reference_slope(sf, routes, lam_beta):
    """Calibrate reference curvature from the same-state HDV-only UE."""
    game = SE.SFGame(
        sf, alpha=0.9, shares=[1.0], lam_beta=lam_beta,
        routes=routes, mode="ue",
    )
    sol = game.solve(tol=1e-8, max_iter=20000, verbose=False)
    if sol["gap"] > 1e-8:
        raise RuntimeError(f"hazard HDV-only reference did not converge: {sol['gap']}")
    return SE.bpr_deriv(sol["x"], game.t0, game.cap)


def realized_critical_edges(scenario, base_ue_x, n_edges=3):
    """Choose critical edges with realized baseline flow, disjoint from disruptions."""
    excluded = set(int(e) for e in scenario.get("closed_edges", []))
    excluded.update(int(e) for e in scenario.get("degraded_edges", []))
    ranked = np.argsort(-np.asarray(base_ue_x, dtype=float))
    selected = [int(e) for e in ranked if int(e) not in excluded and base_ue_x[int(e)] > 1e-9]
    if len(selected) < n_edges:
        selected.extend(int(e) for e in ranked if int(e) not in set(selected) and int(e) not in set(scenario.get("closed_edges", [])))
    return selected[:n_edges]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-mode", choices=("base", "hazard"), default="hazard")
    parser.add_argument("--output", type=Path, default=Path("wildfire_sioux_falls_pilot_v2_results.json"))
    args = parser.parse_args()
    global OUT
    OUT = args.output
    base, scenarios = build_snapshot(od_limit=20, demand_scale=1.0)
    total = sum(item["demand"] for item in base["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game, so = solve_benchmarks(base)
    base_slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    rows = []
    for scenario in scenarios:
        scenario = dict(scenario)
        scenario["critical_edges"] = realized_critical_edges(scenario, ue["x"])
        default_risk = 0.0 if scenario["scenario_id"] == "no_hazard" else 0.25
        scenario["risk_score"] = {
            str(e): float(scenario.get("risk_score", {}).get(str(e), default_risk))
            for e in scenario["critical_edges"]
        }
        sf = apply_snapshot(base, scenario)
        oracle = hazard_oracle(scenario)
        routes = SE.RouteSet(sf["od_list"])
        hazard_slope = (
            hazard_reference_slope(sf, routes, lam_beta)
            if args.reference_mode == "hazard" else base_slope
        )
        # The no-hazard profile is compared against the same reduced-OD base;
        # no committed full-network result is overwritten.
        for alpha in ALPHAS:
            for gamma in GAMMAS:
                warm = {}
                for partition in PARTITIONS:
                    pid = PARTITION_IDS[partition]
                    sol, game = solve_profile(
                        sf,
                        routes,
                        oracle,
                        alpha,
                        SHARES,
                        lam_beta,
                        coalitions=partition,
                        objective_types=TYPES,
                        coordination_gamma=gamma,
                        management_slope=hazard_slope,
                        warm=warm.get(pid),
                        max_rounds=12,
                    )
                    x = np.asarray(sol["x"], dtype=float)
                    risk = sum(
                        float(scenario.get("risk_score", {}).get(str(e), 0.0)) * x[int(e)]
                        for e in scenario.get("critical_edges", [])
                    )
                    critical = [int(e) for e in scenario.get("critical_edges", [])]
                    vc = x[critical] / np.maximum(game.cap[critical], 1e-12) if critical else np.zeros(1)
                    rows.append(
                        {
                            "scenario_id": scenario["scenario_id"],
                            "alpha": float(alpha),
                            "gamma": float(gamma),
                            "partition_id": pid,
                            "partition": [list(block) for block in partition],
                            "n_coalitions": len(partition),
                            "coalition_HHI": coalition_hhi(partition),
                            "J": float(sol["J"]),
                            "critical_risk_weighted_flow": float(risk),
                            "critical_vc_max": float(np.max(vc)),
                            "critical_vc_mean": float(np.mean(vc)),
                            "edge_flow": x.tolist(),
                            "capacity": game.cap.tolist(),
                            "total_av_travel_cost": float(np.sum(company_accounting_costs(sol, game)["company_travel_cost"])),
                            "VI_gap": float(sol["gap"]),
                            "omitted_route_slack": float(sol["max_reduced_cost"]),
                            "route_cost_scale": float(sol["route_cost_scale"]),
                            "certified": bool(sol["certified"]),
                            "iterations": int(sol["iters"]),
                            "cg_rounds": int(sol["cg_rounds"]),
                            "paths": int(sol["n_paths"]),
                            **company_accounting_costs(sol, game),
                        }
                    )
                    warm[pid] = (sol["fH"], sol["fF"])
    result = {
        "status": "pilot_only",
        "planning_data": True,
        "protocol": {
            "network": "Sioux Falls reduced OD pilot",
            "od_limit": 20,
            "demand_scale": 1.0,
            "ordinary_hdv": True,
            "commercial_av_orders_identity_preserving": True,
            "partitions": len(PARTITIONS),
            "alphas": list(ALPHAS),
            "gammas": list(GAMMAS),
            "reference_slope": (
                "same-hazard-state reduced-OD HDV-only UE derivative"
                if args.reference_mode == "hazard"
                else "base-network reduced-OD HDV-only UE derivative reused in every hazard state"
            ),
        },
        "scenarios": scenarios,
        "profiles": rows,
    }
    grouped = {}
    for row in rows:
        grouped.setdefault((row["scenario_id"], row["alpha"], row["gamma"]), []).append(row)
    summaries = []
    for (scenario_id, alpha, gamma), group in sorted(grouped.items()):
        best = min(group, key=lambda r: r["J"])
        singleton = next(r for r in group if r["n_coalitions"] == 4)
        grand = next(r for r in group if r["n_coalitions"] == 1)
        summaries.append(
            {
                "scenario_id": scenario_id,
                "alpha": alpha,
                "gamma": gamma,
                "best_partition_id": best["partition_id"],
                "best_K": best["n_coalitions"],
                "best_J": best["J"],
                "singleton_J": singleton["J"],
                "grand_J": grand["J"],
                "grand_minus_singleton_J": grand["J"] - singleton["J"],
                "best_critical_risk_weighted_flow": best["critical_risk_weighted_flow"],
                "best_regret_of_grand": grand["J"] - best["J"],
            }
        )
    result["summary"] = summaries
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"pilot profiles={len(rows)} certified={sum(r['certified'] for r in rows)}")


if __name__ == "__main__":
    main()
