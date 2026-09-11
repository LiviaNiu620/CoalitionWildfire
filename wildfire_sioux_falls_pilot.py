"""Pilot n=4 wildfire snapshot partition audit on reduced Sioux Falls OD data.

This writes planning/pilot evidence only. It is not the manuscript's full
wildfire result until the all-OD scan and its validator are complete.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from coalition_private_costs import company_accounting_costs
from sf_coalition_partition_scan import PARTITIONS, PARTITION_IDS, coalition_hhi
from sf_revision_common import LAM_TIMES_M, SE, network_oracle, solve_benchmarks, solve_profile
from wildfire_scenarios import apply_snapshot, build_snapshot, hazard_oracle


OUT = Path("wildfire_sioux_falls_pilot_results.json")
SHARES = np.full(4, 0.25)
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
ALPHAS = (0.5, 0.9)
GAMMAS = (0.5, 1.0)


def main() -> None:
    base, scenarios = build_snapshot(od_limit=20, demand_scale=1.0)
    total = sum(item["demand"] for item in base["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    base_oracle = network_oracle()
    ue, ue_game, so = solve_benchmarks(base)
    base_slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    rows = []
    for scenario in scenarios:
        sf = apply_snapshot(base, scenario)
        oracle = hazard_oracle(scenario)
        routes = SE.RouteSet(sf["od_list"])
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
                        management_slope=base_slope,
                        warm=warm.get(pid),
                        max_rounds=12,
                    )
                    x = np.asarray(sol["x"], dtype=float)
                    risk = sum(
                        float(scenario.get("risk_score", {}).get(str(e), 0.0)) * x[int(e)]
                        for e in scenario.get("critical_edges", [])
                    )
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
            "reference_slope": "base reduced-OD HDV UE derivative",
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
