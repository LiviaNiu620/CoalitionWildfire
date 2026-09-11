"""Exact full-OD identity certificate for the no-hazard adapter."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from wildfire_scenarios import apply_snapshot, build_snapshot, hazard_oracle


OUT = Path("wildfire_nohazard_identity_certificate.json")


def main() -> None:
    base, scenarios = build_snapshot(od_limit=None, demand_scale=1.0)
    no_hazard = next(s for s in scenarios if s["scenario_id"] == "no_hazard")
    adapted = apply_snapshot(base, no_hazard)
    graph, edge_idx = hazard_oracle(no_hazard)
    checks = {
        "edge_count_equal": base["edges"] == adapted["edges"] == len(edge_idx),
        "t0_exact": bool(np.array_equal(base["t0"], adapted["t0"])),
        "capacity_exact": bool(np.array_equal(base["cap"], adapted["cap"])),
        "od_count_equal": len(base["od_list"]) == len(adapted["od_list"]) == 528,
        "od_demand_exact": all(a["demand"] == b["demand"] for a, b in zip(base["od_list"], adapted["od_list"])),
        "od_routes_exact": all(a["paths"] == b["paths"] for a, b in zip(base["od_list"], adapted["od_list"])),
        "no_closed_edges": no_hazard["closed_edges"] == [],
        "oracle_edge_count": graph.number_of_edges() == 76,
    }
    result = {
        "status": "full_od_input_identity",
        "network": "Sioux Falls",
        "od_count": 528,
        "checks": checks,
        "passed": all(checks.values()),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("no-hazard adapter identity failed")


if __name__ == "__main__":
    main()
