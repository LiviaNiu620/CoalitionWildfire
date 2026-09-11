"""Validator for full-OD endpoint Pareto/epsilon diagnostics."""
from __future__ import annotations

import json
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_endpoint_policy_analysis.json").read_text())
assert data["status"] == "endpoint_policy_analysis"
assert data["planning_data"] is True
assert len(data["states"]) == 8
assert all(s["pareto_size"] >= 1 for s in data["states"])
assert all(len(s["epsilon_choices"]) == 4 for s in data["states"])
print("wildfire endpoint policy analysis validation passed")
