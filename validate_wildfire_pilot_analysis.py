"""Fail-closed validator for pilot outcome-vector analysis."""
from __future__ import annotations

import json
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_pilot_outcome_analysis.json").read_text())
assert data["status"] == "pilot_outcome_analysis"
assert data["planning_data"] is True
assert len(data["states"]) == 16
assert all(s["pareto_size"] >= 1 for s in data["states"])
assert all(s["critical_vc_range"][0] <= s["critical_vc_range"][1] for s in data["states"])
print("wildfire pilot outcome analysis validation passed")
