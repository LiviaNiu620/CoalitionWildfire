"""Fail-closed validator for the same-state-slope wildfire pilot."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "wildfire_sioux_falls_pilot_v2_results.json"
data = json.loads(RAW.read_text())
assert data["status"] == "pilot_only"
assert data["planning_data"] is True
assert data["protocol"]["reference_slope"].startswith("same-hazard-state")
profiles = data["profiles"]
assert len(profiles) == 240
assert len(data["summary"]) == 16
assert all(row["certified"] for row in profiles)
assert all(row["VI_gap"] <= 1e-6 for row in profiles)
assert all(math.isfinite(float(row["J"])) for row in profiles)
print("wildfire same-state-slope pilot validation passed")
