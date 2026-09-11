"""Checkpointed full n=4 partition audit for full-OD wildfire snapshots."""
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


OUT = Path("wildfire_full_partition_audit_results.json")
SHARES = np.full(4, 0.25)
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
PART_BY_ID = {PARTITION_IDS[p]: p for p in PARTITIONS}


def write_checkpoint(meta, rows, status):
    payload = dict(meta)
    payload.update({"status": status, "profiles": rows})
    OUT.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenarios", nargs="+", default=["mild_capacity_degradation", "critical_edge_closure"])
    parser.add_argument("--alphas", nargs="+", type=float, default=[0.5, 0.9])
    parser.add_argument("--gammas", nargs="+", type=float, default=[0.5, 1.0])
    parser.add_argument("--max-profiles", type=int, default=10_000)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    base, catalog = build_snapshot(None, 1.0)
    total = sum(o["demand"] for o in base["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game, so = solve_benchmarks(base)
    scenario_map = {s["scenario_id"]: s for s in catalog}
    missing = sorted(set(args.scenarios) - set(scenario_map))
    if missing:
        raise SystemExit(f"unknown scenarios: {missing}")
    rows = []
    if args.resume and OUT.exists():
        old = json.loads(OUT.read_text())
        rows = list(old.get("profiles", []))
    done = {(r["scenario_id"], r["alpha"], r["gamma"], r["partition_id"]) for r in rows}
    target = min(args.max_profiles, len(args.scenarios) * len(args.alphas) * len(args.gammas) * len(PARTITIONS))
    meta = {
        "status": "running",
        "planning_data": True,
        "network": "Sioux Falls full OD",
        "scenarios": args.scenarios,
        "alphas": args.alphas,
        "gammas": args.gammas,
        "partition_count": len(PARTITIONS),
        "reference_slope": "same-hazard-state HDV-only UE derivative",
        "metrics": ["TSTT", "total_av_travel_cost", "critical_vc_max", "critical_vc_mean", "critical_risk_weighted_flow"],
    }
    for scenario_id in args.scenarios:
        scenario0 = scenario_map[scenario_id]
        scenario = dict(scenario0)
        scenario["critical_edges"] = realized_critical_edges(scenario, ue["x"])
        default_risk = 0.0 if scenario_id == "no_hazard" else 0.25
        scenario["risk_score"] = {
            str(e): float(scenario.get("risk_score", {}).get(str(e), default_risk))
            for e in scenario["critical_edges"]
        }
        sf = apply_snapshot(base, scenario)
        oracle = hazard_oracle(scenario)
        ref_routes = SE.RouteSet(sf["od_list"])
        slope = hazard_reference_slope(sf, ref_routes, lam_beta)
        routes = SE.RouteSet(sf["od_list"])
        for alpha in args.alphas:
            for gamma in args.gammas:
                warm, warm_keys = {}, {}
                for partition in PARTITIONS:
                    pid = PARTITION_IDS[partition]
                    key = (scenario_id, alpha, gamma, pid)
                    if key in done:
                        continue
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
                        "scenario_id": scenario_id,
                        "alpha": float(alpha),
                        "gamma": float(gamma),
                        "partition_id": pid,
                        "partition": [list(block) for block in partition],
                        "n_coalitions": len(partition),
                        "J": float(sol["J"]),
                        "total_av_travel_cost": float(accounting["total_av_travel_cost"]),
                        "critical_vc_max": float(np.max(vc)),
                        "critical_vc_mean": float(np.mean(vc)),
                        "critical_risk_weighted_flow": float(sum(scenario["risk_score"][str(e)] * x[int(e)] for e in critical)),
                        "edge_flow": x.tolist(),
                        "capacity": game.cap.tolist(),
                        "VI_gap": float(sol["gap"]),
                        "omitted_route_slack": float(sol["max_reduced_cost"]),
                        "route_cost_scale": float(sol["route_cost_scale"]),
                        "certified": bool(sol["certified"]),
                        "iterations": int(sol["iters"]),
                        "cg_rounds": int(sol["cg_rounds"]),
                        "paths": int(sol["n_paths"]),
                    })
                    done.add(key)
                    warm[pid] = (sol["fH"], sol["fF"])
                    warm_keys[pid] = routes.keys()
                    write_checkpoint(meta, rows, "running")
                    print(f"[{len(rows)}/{target}] {key} gap={sol['gap']:.2e} cert={sol['certified']}", flush=True)
                    if len(rows) >= target:
                        write_checkpoint(meta, rows, "complete")
                        return
    write_checkpoint(meta, rows, "complete")


if __name__ == "__main__":
    main()
