"""Pareto and epsilon-constraint diagnostics for full-OD endpoint profiles."""
from __future__ import annotations

import json
from pathlib import Path


RAW = Path("wildfire_full_endpoint_audit_results.json")
OUT = Path("wildfire_endpoint_policy_analysis.json")


def dominates(a, b):
    keys = ("J", "total_av_travel_cost", "critical_vc_max")
    return all(a[k] <= b[k] for k in keys) and any(a[k] < b[k] for k in keys)


def main() -> None:
    data = json.loads(RAW.read_text())
    groups = {}
    for row in data["profiles"]:
        groups.setdefault((row["scenario_id"], row["alpha"], row["gamma"]), []).append(row)
    states = []
    for key, rows in sorted(groups.items()):
        front = [r for r in rows if not any(dominates(o, r) for o in rows if o is not r)]
        t_best = min(rows, key=lambda r: r["J"])
        # Threshold grid is a diagnostic; no threshold is claimed as universal.
        epsilon_choices = []
        for vc_limit in (1.70, 1.73, 1.75, 1.80):
            feasible = [r for r in rows if r["critical_vc_max"] <= vc_limit]
            choice = min(feasible, key=lambda r: r["J"]) if feasible else None
            epsilon_choices.append({
                "critical_vc_limit": vc_limit,
                "feasible": choice is not None,
                "selected_partition_id": choice["partition_id"] if choice else None,
                "selected_K": choice["n_coalitions"] if choice else None,
                "selected_J": choice["J"] if choice else None,
            })
        states.append({
            "scenario_id": key[0],
            "alpha": key[1],
            "gamma": key[2],
            "tstt_best_partition_id": t_best["partition_id"],
            "tstt_best_K": t_best["n_coalitions"],
            "pareto_partition_ids": sorted(r["partition_id"] for r in front),
            "pareto_size": len(front),
            "epsilon_choices": epsilon_choices,
        })
    result = {
        "status": "endpoint_policy_analysis",
        "planning_data": True,
        "note": "Only four representative partitions are audited; epsilon limits are diagnostics, not universal policy thresholds.",
        "states": states,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
