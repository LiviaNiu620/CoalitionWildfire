"""Validator for completed rows in the full-OD endpoint hazard audit."""
from __future__ import annotations

import json
import math
from pathlib import Path


data = json.loads((Path(__file__).resolve().parent / "wildfire_full_endpoint_audit_results.json").read_text())
rows = data["profiles"]
assert data["network"] == "Sioux Falls full OD"
assert len(rows) >= 1
assert all(r["certified"] for r in rows)
assert all(r["VI_gap"] <= 1e-6 for r in rows)
assert all(r["omitted_route_slack"] / max(r["route_cost_scale"], 1e-12) <= 1e-6 for r in rows)
assert all(math.isfinite(r["J"]) and math.isfinite(r["critical_vc_max"]) for r in rows)
print(f"wildfire full-OD endpoint rows validated: {len(rows)}/{data['target_profiles']}")
