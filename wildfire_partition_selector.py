"""Formal epsilon-constrained selector over complete partition-audit rows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


RAW = Path("wildfire_full_partition_audit_results.json")
OUT = Path("wildfire_partition_selection.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=RAW)
    parser.add_argument("--critical-vc-limit", type=float, default=None)
    parser.add_argument("--order-cost-limit", type=float, default=None)
    args = parser.parse_args()
    raw_path = args.raw
    data = json.loads(raw_path.read_text())
    rows = data["profiles"]
    groups = {}
    for row in rows:
        groups.setdefault((row["scenario_id"], row["alpha"], row["gamma"]), []).append(row)
    states = []
    for key, group in sorted(groups.items()):
        feasible = list(group)
        if args.critical_vc_limit is not None:
            feasible = [r for r in feasible if r["critical_vc_max"] <= args.critical_vc_limit]
        if args.order_cost_limit is not None:
            feasible = [r for r in feasible if r["total_av_travel_cost"] <= args.order_cost_limit]
        choice = min(feasible, key=lambda r: r["J"]) if feasible else None
        states.append({
            "scenario_id": key[0],
            "alpha": key[1],
            "gamma": key[2],
            "critical_vc_limit": args.critical_vc_limit,
            "order_cost_limit": args.order_cost_limit,
            "feasible_count": len(feasible),
            "selected_partition_id": choice["partition_id"] if choice else None,
            "selected_K": choice["n_coalitions"] if choice else None,
            "selected_J": choice["J"] if choice else None,
            "selected_critical_vc_max": choice["critical_vc_max"] if choice else None,
            "selected_total_av_travel_cost": choice["total_av_travel_cost"] if choice else None,
        })
    result = {
        "status": "epsilon_selection",
        "planning_data": True,
        "source": str(raw_path),
        "partition_count_required": 15,
        "states": states,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
