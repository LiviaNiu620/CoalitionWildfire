"""Validator for the multiple-corridor full partition audit."""
from __future__ import annotations

import json
import math
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_full_partition_multiple_results.json").read_text())
rows = data["profiles"]
assert data["status"] == "complete"
assert len(rows) == 60
assert all(r["scenario_id"] == "multiple_corridor_degradation" for r in rows)
assert all(r["certified"] for r in rows)
assert all(r["VI_gap"] <= 1e-6 for r in rows)
assert all(math.isfinite(r["J"]) and math.isfinite(r["critical_vc_max"]) for r in rows)
print("wildfire multiple-corridor full partition audit validation passed: 60/60")
