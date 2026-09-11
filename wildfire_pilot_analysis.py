"""Outcome-vector, Pareto, and epsilon-constraint analysis for pilot profiles."""
from __future__ import annotations

import json
from pathlib import Path


RAW = Path("wildfire_sioux_falls_pilot_v2_results.json")
OUT = Path("wildfire_pilot_outcome_analysis.json")


def dominates(a, b):
    """All metrics are minimized; strict improvement in at least one."""
    keys = ("J", "total_av_travel_cost", "critical_vc_max", "critical_risk_weighted_flow")
    return all(float(a[k]) <= float(b[k]) for k in keys) and any(
        float(a[k]) < float(b[k]) for k in keys
    )


def pareto(rows):
    return [r for r in rows if not any(dominates(other, r) for other in rows if other is not r)]


def main() -> None:
    data = json.loads(RAW.read_text())
    groups = {}
    for row in data["profiles"]:
        groups.setdefault((row["scenario_id"], row["alpha"], row["gamma"]), []).append(row)
    summaries = []
    for key, rows in sorted(groups.items()):
        front = pareto(rows)
        tmin = min(rows, key=lambda r: r["J"])
        summaries.append({
            "scenario_id": key[0],
            "alpha": key[1],
            "gamma": key[2],
            "tstt_best_partition_id": tmin["partition_id"],
            "tstt_best_K": tmin["n_coalitions"],
            "pareto_partition_ids": sorted(r["partition_id"] for r in front),
            "pareto_size": len(front),
            "critical_vc_min_partition_id": min(rows, key=lambda r: r["critical_vc_max"])["partition_id"],
            "critical_vc_min": min(r["critical_vc_max"] for r in rows),
            "tstt_range": [min(r["J"] for r in rows), max(r["J"] for r in rows)],
            "critical_vc_range": [min(r["critical_vc_max"] for r in rows), max(r["critical_vc_max"] for r in rows)],
        })
    result = {
        "status": "pilot_outcome_analysis",
        "planning_data": True,
        "objective_direction": "minimize TSTT, AV order cost, critical v/c, and risk diagnostic",
        "note": "Pareto diagnostics do not select a single policy without declared epsilon constraints.",
        "states": summaries,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
