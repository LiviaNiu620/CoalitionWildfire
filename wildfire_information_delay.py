"""Pilot separating perceived public risk information from true evaluation."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sf_coalition_partition_scan import PARTITIONS, PARTITION_IDS
from sf_revision_common import LAM_TIMES_M, SE, solve_benchmarks, solve_profile
from wildfire_scenarios import apply_snapshot, build_snapshot, hazard_oracle
from wildfire_sioux_falls_pilot import hazard_reference_slope, realized_critical_edges


OUT = Path("wildfire_information_delay_results.json")


def main() -> None:
    base, catalog = build_snapshot(od_limit=20, demand_scale=1.0)
    scenario0 = next(s for s in catalog if s["scenario_id"] == "critical_edge_closure")
    scenario = dict(scenario0)
    base_ue, _, _ = solve_benchmarks(base)
    scenario["critical_edges"] = realized_critical_edges(scenario, base_ue["x"])
    scenario["risk_score"] = {str(e): 1.0 for e in scenario["critical_edges"]}
    sf = apply_snapshot(base, scenario)
    routes = SE.RouteSet(sf["od_list"])
    oracle = hazard_oracle(scenario)
    total = sum(o["demand"] for o in sf["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    slope = hazard_reference_slope(sf, routes, lam_beta)
    true_penalty = np.zeros(sf["edges"])
    for edge in scenario["critical_edges"]:
        true_penalty[int(edge)] = 2.0
    rows = []
    for delay_name, perceived_scale in (("perfect", 1.0), ("stale_half", 0.5), ("stale_none", 0.0)):
        perceived = perceived_scale * true_penalty
        for alpha in (0.5, 0.9):
            for gamma in (0.5, 1.0):
                for partition in PARTITIONS:
                    sol, game = solve_profile(
                        sf, routes, oracle, alpha, np.full(4, 0.25), lam_beta,
                        coalitions=partition, objective_types=np.array([0.5, 0.5, 1.5, 1.5]),
                        coordination_gamma=gamma, management_slope=slope,
                        perceived_edge_penalty=perceived,
                        solver_method="potential", solver_max_iter=12000,
                    )
                    x = np.asarray(sol["x"], dtype=float)
                    rows.append({
                        "information_case": delay_name,
                        "perceived_penalty_scale": perceived_scale,
                        "alpha": alpha,
                        "gamma": gamma,
                        "partition_id": PARTITION_IDS[partition],
                        "n_coalitions": len(partition),
                        "J_physical": float(sol["J"]),
                        "true_risk_cost": float(np.dot(true_penalty, x)),
                        "critical_vc_max": float(np.max(x[scenario["critical_edges"]] / np.maximum(game.cap[scenario["critical_edges"]], 1e-12))),
                        "VI_gap": float(sol["gap"]),
                        "certified": bool(sol["certified"]),
                    })
    result = {
        "status": "information_delay_pilot",
        "planning_data": True,
        "network": "Sioux Falls reduced OD 20",
        "true_state": "critical-edge closure with true public risk penalty",
        "note": "Only the perceived common risk penalty is delayed; closed-edge discovery is a later extension.",
        "profiles": rows,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"information-delay profiles={len(rows)} certified={sum(r['certified'] for r in rows)}")


if __name__ == "__main__":
    main()
