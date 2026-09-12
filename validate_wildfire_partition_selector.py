"""Validator for epsilon-constrained partition selection."""
from __future__ import annotations

import json
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_partition_selection.json").read_text())
assert data["status"] == "epsilon_selection"
assert data["planning_data"] is True
assert data["partition_count_required"] == 15
assert len(data["states"]) in {4, 8}
assert all(s["feasible_count"] >= 0 for s in data["states"])
print("wildfire epsilon-constrained partition selector validation passed")
