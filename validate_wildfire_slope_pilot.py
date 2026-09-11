"""Fail-closed validator for the V0/V2 reference-slope comparison."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / "wildfire_slope_robustness_pilot.json").read_text())
assert data["status"] == "pilot_comparison"
assert data["planning_data"] is True
assert len(data["rows"]) == 16
assert data["v0_max_grand_regret"] >= data["v2_max_grand_regret"]
assert data["best_K_change_count"] >= 0
assert data["partition_change_count"] >= 0
print("wildfire slope robustness pilot validation passed")
