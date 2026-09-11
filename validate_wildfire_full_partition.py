"""Fail-closed validator for the complete 120-profile hazard audit."""
from __future__ import annotations

import json
import math
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_full_partition_audit_results.json").read_text())
rows = data["profiles"]
assert data["status"] == "complete"
assert data["network"] == "Sioux Falls full OD"
assert len(rows) == 120
assert len({(r["scenario_id"], r["alpha"], r["gamma"]) for r in rows}) == 8
assert all(sum(r["scenario_id"] == s for r in rows) == 60 for s in ("mild_capacity_degradation", "critical_edge_closure"))
assert all(r["certified"] for r in rows)
assert all(r["VI_gap"] <= 1e-6 for r in rows)
assert all(r["omitted_route_slack"] / max(r["route_cost_scale"], 1e-12) <= 1e-6 for r in rows)
for row in rows:
    for key in ("J", "total_av_travel_cost", "critical_vc_max", "critical_vc_mean", "critical_risk_weighted_flow"):
        assert math.isfinite(float(row[key])), (row["scenario_id"], row["alpha"], row["gamma"], row["partition_id"], key)
print("wildfire full 15-partition hazard audit validation passed: 120/120")
