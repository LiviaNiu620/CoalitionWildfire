"""Fail-closed validation for objective/governance sensitivity experiments."""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def finite(v):
    if isinstance(v, dict):
        for x in v.values(): finite(x)
    elif isinstance(v, list):
        for x in v: finite(x)
    elif isinstance(v, float) and not math.isfinite(v):
        raise AssertionError("nonfinite result")

def main():
    raw = json.loads((ROOT / "sf_coalition_sensitivity_results.json").read_text())
    summary = json.loads((ROOT / "sf_coalition_sensitivity_summary.json").read_text())
    rows = raw["profiles"]
    expected = 2 * 5 * 3 * 3 * 2
    if len(rows) != expected: raise AssertionError(f"expected {expected} rows, found {len(rows)}")
    if not all(r["certified"] for r in rows): raise AssertionError("uncertified sensitivity profile")
    if max(r["VI_gap"] for r in rows) > 1.0001e-6: raise AssertionError("VI tolerance exceeded")
    if max(r["omitted_route_slack"] / r["route_cost_scale"] for r in rows) > 1.0001e-6:
        raise AssertionError("omitted-route tolerance exceeded")
    if {r["type_design"] for r in rows} != {"wide", "baseline", "narrow"}: raise AssertionError("type grid incomplete")
    if {r["slope_scale"] for r in rows} != {0.5, 1.0, 2.0}: raise AssertionError("slope grid incomplete")
    finite(raw); finite(summary)
    print("coalition sensitivity package validation passed")

if __name__ == "__main__": main()
