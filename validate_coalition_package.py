"""Fail-closed validation for the coalition-objective experiment package."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def finite(value, path="root"):
    if isinstance(value, dict):
        for key, item in value.items():
            finite(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            finite(item, f"{path}[{index}]")
    elif isinstance(value, float) and not math.isfinite(value):
        raise AssertionError(f"nonfinite value at {path}")


def main():
    raw = read("sf_coalition_composition_results.json")
    summary = read("sf_coalition_composition_summary.json")
    rows = raw["profiles"]
    if len(rows) != 40:
        raise AssertionError(f"expected 40 rows, found {len(rows)}")
    expected = {
        (alpha, gamma, partition)
        for alpha in (0.5, 0.9)
        for gamma in (0.0, 0.25, 0.5, 0.75, 1.0)
        for partition in ("singletons", "assortative", "mixed", "grand")
    }
    keys = {(row["alpha"], row["gamma"], row["partition"]) for row in rows}
    if keys != expected:
        raise AssertionError("coalition profile grid is incomplete")
    if not all(row["certified"] for row in rows):
        raise AssertionError("one or more coalition equilibria are uncertified")
    if max(row["VI_gap"] for row in rows) > 1.0001e-6:
        raise AssertionError("VI tolerance exceeded")
    if max(row["omitted_route_slack"] / row["route_cost_scale"] for row in rows) > 1.0001e-6:
        raise AssertionError("omitted-route tolerance exceeded")
    if any(abs(row["company_HHI"] - 0.25) > 1e-14 for row in rows):
        raise AssertionError("company HHI changed")
    for row in rows:
        target = 0.25 if row["partition"] == "singletons" else (
            1.0 if row["partition"] == "grand" else 0.5
        )
        if abs(row["coalition_HHI"] - target) > 1e-14:
            raise AssertionError("coalition HHI does not match the declared partition")
    if any(item["max_relative_TSTT_difference"] > 1e-8
           for item in summary["partition_invariance_gamma0"]):
        raise AssertionError("gamma=0 partition invariance failed in TSTT")
    if any(item["max_abs_edge_flow_difference"] > 1e-3
           for item in summary["partition_invariance_gamma0"]):
        raise AssertionError("gamma=0 partition invariance failed in edge flow")
    finite(raw, "raw")
    finite(summary, "summary")
    for name in (
        "response_paper/figures/fig_sf_coalition_composition.pdf",
        "response_paper/figures/fig_sf_coalition_composition.png",
    ):
        path = ROOT / name
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing figure {name}")
    print("coalition package validation passed")


if __name__ == "__main__":
    main()
