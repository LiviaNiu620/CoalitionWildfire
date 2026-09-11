"""Validation for the company-count and coalition-degree experiment."""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def finite(v):
    if isinstance(v, dict):
        for x in v.values(): finite(x)
    elif isinstance(v, list):
        for x in v: finite(x)
    elif isinstance(v, float) and not math.isfinite(v):
        raise AssertionError("nonfinite result")


def main():
    raw = read("sf_coalition_degree_results.json")
    summary = read("sf_coalition_degree_summary.json")
    rows = raw["profiles"]
    expected = sum({2: 2, 4: 3, 8: 4}[n] for n in (2, 4, 8)) * 2 * 3
    if len(rows) != expected:
        raise AssertionError(f"expected {expected} rows, found {len(rows)}")
    if not all(r["certified"] for r in rows):
        raise AssertionError("uncertified coalition-degree profile")
    if max(r["VI_gap"] for r in rows) > 1.0001e-6:
        raise AssertionError("VI tolerance exceeded")
    if max(r["omitted_route_slack"] / r["route_cost_scale"] for r in rows) > 1.0001e-6:
        raise AssertionError("omitted-route tolerance exceeded")
    if {r["n_companies"] for r in rows} != {2, 4, 8}:
        raise AssertionError("company-count grid incomplete")
    if any(r["n_coalitions"] != len(r["coalitions"]) for r in rows):
        raise AssertionError("coalition count inconsistent")
    if any(abs(r["company_HHI"] - 1.0 / r["n_companies"]) > 1e-14 for r in rows):
        raise AssertionError("company HHI mismatch")
    finite(raw)
    finite(summary)
    print("coalition degree package validation passed")


if __name__ == "__main__":
    main()
