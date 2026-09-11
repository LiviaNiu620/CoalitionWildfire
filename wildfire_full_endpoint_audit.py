"""Checkpointed full-OD hazard audit over representative coalition endpoints."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coalition_private_costs import company_accounting_costs
from sf_coalition_partition_scan import PARTITIONS, PARTITION_IDS
from sf_revision_common import LAM_TIMES_M, SE, solve_benchmarks, solve_profile
from wildfire_scenarios import apply_snapshot, build_snapshot, hazard_oracle
from wildfire_sioux_falls_pilot import hazard_reference_slope, realized_critical_edges


OUT = Path("wildfire_full_endpoint_audit_results.json")
SHARES = np.full(4, 0.25)
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
PIDS = ("P01", "P03", "P06", "P15")
PART_BY_ID = {PARTITION_IDS[p]: p for p in PARTITIONS}


def checkpoint(rows, target):
    OUT.write_text(json.dumps({
        "status": "running" if len(rows) < target else "complete",
        "planning_data": True,
        "network": "Sioux Falls full OD",
        "target_profiles": target,
        "profiles": rows,
    }, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-profiles", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    base, catalog = build_snapshot(None, 1.0)
    total = sum(o["demand"] for o in base["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game, so = solve_benchmarks(base)
    scenarios = [s for s in catalog if s["scenario_id"] in {"mild_capacity_degradation", "critical_edge_closure"}]
    rows = []
    if args.resume and OUT.exists():
        prior = json.loads(OUT.read_text())
        rows = prior.get("profiles", [])
    done = {(r["scenario_id"], r["alpha"], r["gamma"], r["partition_id"]) for r in rows}
    target = min(args.max_profiles, len(scenarios) * 2 * 2 * len(PIDS))
    for scenario0 in scenarios:
        scenario = dict(scenario0)
        scenario["critical_edges"] = realized_critical_edges(scenario, ue["x"])
        sf = apply_snapshot(base, scenario)
        oracle = hazard_oracle(scenario)
        reference_routes = SE.RouteSet(sf["od_list"])
        slope = hazard_reference_slope(sf, reference_routes, lam_beta)
        routes = SE.RouteSet(sf["od_list"])
        warm, warm_keys = {}, {}
        for alpha in (0.5, 0.9):
            for gamma in (0.5, 1.0):
                for pid in PIDS:
                    key = (scenario["scenario_id"], alpha, gamma, pid)
                    if key in done:
                        continue
                    partition = PART_BY_ID[pid]
                    sol, game = solve_profile(
                        sf, routes, oracle, alpha, SHARES, lam_beta,
                        coalitions=partition, objective_types=TYPES,
                        coordination_gamma=gamma, management_slope=slope,
                        warm=warm.get(pid), warm_keys=warm_keys.get(pid),
                        solver_method="potential", solver_max_iter=12000,
                    )
                    accounting = company_accounting_costs(sol, game)
                    x = np.asarray(sol["x"], dtype=float)
                    critical = scenario["critical_edges"]
                    vc = x[critical] / np.maximum(game.cap[critical], 1e-12)
                    rows.append({
                        "scenario_id": scenario["scenario_id"],
                        "alpha": alpha,
                        "gamma": gamma,
                        "partition_id": pid,
                        "n_coalitions": len(partition),
                        "J": float(sol["J"]),
                        "total_av_travel_cost": float(accounting["total_av_travel_cost"]),
                        "critical_vc_max": float(np.max(vc)),
                        "critical_vc_mean": float(np.mean(vc)),
                        "VI_gap": float(sol["gap"]),
                        "omitted_route_slack": float(sol["max_reduced_cost"]),
                        "route_cost_scale": float(sol["route_cost_scale"]),
                        "certified": bool(sol["certified"]),
                        "iterations": int(sol["iters"]),
                        "paths": int(sol["n_paths"]),
                    })
                    done.add(key)
                    warm[pid] = (sol["fH"], sol["fF"])
                    warm_keys[pid] = routes.keys()
                    checkpoint(rows, target)
                    print(f"[{len(rows)}/{target}] {key} gap={sol['gap']:.2e} cert={sol['certified']}", flush=True)
                    if len(rows) >= target:
                        return
    checkpoint(rows, target)


if __name__ == "__main__":
    main()
