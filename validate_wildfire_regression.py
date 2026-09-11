"""Fail-closed validator for the full-OD no-hazard regression."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / "wildfire_nohazard_regression_results.json").read_text())
assert data["status"] == "regression"
assert data["network"] == "Sioux Falls full OD" or data["network"].startswith("Sioux Falls OD limit ")
assert len(data["profiles"]) >= 1
assert data["certified_rows"] == len(data["profiles"])
assert data["max_VI_gap"] <= 1e-6
if data["baseline_comparable"]:
    assert data["max_relative_delta_J"] <= 1e-8
assert data["passed"] is True
print("wildfire full-OD no-hazard regression passed")
