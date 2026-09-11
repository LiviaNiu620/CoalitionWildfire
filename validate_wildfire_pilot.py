"""Fail-closed validator for the reduced-OD wildfire Sioux Falls pilot."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "wildfire_sioux_falls_pilot_results.json"


def main() -> None:
    data = json.loads(RAW.read_text())
    assert data["status"] == "pilot_only"
    assert data["planning_data"] is True
    assert len(data["scenarios"]) == 4
    profiles = data["profiles"]
    assert len(profiles) == 240
    assert len(data["summary"]) == 16
    assert all(row["certified"] for row in profiles)
    assert all(row["VI_gap"] <= 1e-6 for row in profiles)
    assert all(
        row["omitted_route_slack"] / max(row["route_cost_scale"], 1e-12) <= 1e-6
        for row in profiles
    )
    for row in profiles:
        for key in ("J", "critical_risk_weighted_flow", "VI_gap", "omitted_route_slack"):
            assert math.isfinite(float(row[key])), (row["scenario_id"], key)
    assert {row["n_coalitions"] for row in profiles} == {1, 2, 3, 4}
    print("wildfire Sioux Falls pilot validation passed")


if __name__ == "__main__":
    main()
