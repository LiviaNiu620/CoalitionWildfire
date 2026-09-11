"""Compare V0 no-hazard and V2 same-hazard reference-slope pilots."""
from __future__ import annotations

import json
from pathlib import Path


V0 = Path("wildfire_sioux_falls_pilot_v0_spatial_results.json")
V2 = Path("wildfire_sioux_falls_pilot_v2_results.json")
OUT = Path("wildfire_slope_robustness_pilot.json")


def main() -> None:
    old = json.loads(V0.read_text())["summary"]
    new = json.loads(V2.read_text())["summary"]
    old_by = {(r["scenario_id"], r["alpha"], r["gamma"]): r for r in old}
    new_by = {(r["scenario_id"], r["alpha"], r["gamma"]): r for r in new}
    rows = []
    for key in sorted(new_by):
        a, b = old_by[key], new_by[key]
        rows.append({
            "scenario_id": key[0],
            "alpha": key[1],
            "gamma": key[2],
            "v0_best_K": a["best_K"],
            "v0_best_partition_id": a["best_partition_id"],
            "v2_best_K": b["best_K"],
            "v2_best_partition_id": b["best_partition_id"],
            "best_K_changed": a["best_K"] != b["best_K"],
            "partition_changed": a["best_partition_id"] != b["best_partition_id"],
            "v0_grand_regret": a["best_regret_of_grand"],
            "v2_grand_regret": b["best_regret_of_grand"],
            "v2_minus_v0_best_J": b["best_J"] - a["best_J"],
        })
    result = {
        "status": "pilot_comparison",
        "planning_data": True,
        "v0": "base-network HDV-only UE derivative reused in every hazard state",
        "v2": "same-hazard-state HDV-only UE derivative",
        "rows": rows,
        "best_K_change_count": sum(r["best_K_changed"] for r in rows),
        "partition_change_count": sum(r["partition_changed"] for r in rows),
        "v0_max_grand_regret": max(r["v0_grand_regret"] for r in rows),
        "v2_max_grand_regret": max(r["v2_grand_regret"] for r in rows),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
