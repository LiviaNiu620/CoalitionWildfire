"""Validator for the perceived-state/true-state information-delay pilot."""
from __future__ import annotations

import json
import math
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_information_delay_results.json").read_text())
rows = data["profiles"]
assert data["status"] == "information_delay_pilot"
assert data["planning_data"] is True
assert len(rows) == 180
assert {r["information_case"] for r in rows} == {"perfect", "stale_half", "stale_none"}
assert all(r["certified"] for r in rows)
assert all(math.isfinite(r["J_physical"]) and math.isfinite(r["true_risk_cost"]) for r in rows)
print("wildfire information-delay pilot validation passed")
