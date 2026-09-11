"""Fail-closed consistency checks for the major-revision evidence package."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def assert_finite(value, path="root"):
    if isinstance(value, dict):
        for key, item in value.items():
            assert_finite(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            assert_finite(item, f"{path}[{index}]")
    elif isinstance(value, float) and not math.isfinite(value):
        raise AssertionError(f"nonfinite value at {path}")


def certified(rows, expected):
    if len(rows) != expected:
        raise AssertionError(f"expected {expected} rows, found {len(rows)}")
    if not all(row["certified"] for row in rows):
        raise AssertionError("one or more equilibrium rows are uncertified")
    if max(row["VI_gap"] for row in rows) > 1.0001e-6:
        raise AssertionError("VI tolerance exceeded")
    if max(row["omitted_route_slack"] / row["route_cost_scale"] for row in rows) > 1.0001e-6:
        raise AssertionError("omitted-route tolerance exceeded")


def main():
    decomposition = read("sf_response_decomposition_results.json")
    alpha = read("sf_alpha_screening_results.json")
    spatial = read("sf_spatial_od_results.json")
    anchors = read("sf_bpr_anchor_audit_results.json")
    coalition_sensitivity = read("sf_coalition_sensitivity_results.json")
    coalition_partition = read("sf_coalition_partition_scan_results.json")
    coalition_pooling = read("sf_coalition_pooling_results.json")
    certified(decomposition["profiles"], 600)
    certified(alpha["profiles"], 500)
    certified(spatial["profiles"], 82)
    certified(anchors["profiles"], 10)
    certified(coalition_sensitivity["profiles"], 180)
    certified(coalition_partition["profiles"], 60)
    certified(coalition_pooling["profiles"], 4)

    regimes = {row["regime"] for row in decomposition["profiles"]}
    if regimes != {"mass_linked_prior", "full_information", "common_half"}:
        raise AssertionError("decomposition regimes are incomplete")
    if {row["alpha"] for row in alpha["profiles"]} != {0.1, 0.3, 0.5, 0.7, 0.9}:
        raise AssertionError("alpha design is incomplete")
    if any(abs(row["HHI"] - 0.25) > 1e-10 for row in spatial["profiles"]):
        raise AssertionError("fixed-HHI constraint failed")
    if any(row["max_fleet_mass_error"] > 1e-6 for row in spatial["profiles"]):
        raise AssertionError("equal fleet-mass constraint failed")

    for name in (
        "sf_response_decomposition_summary.json",
        "sf_alpha_screening_summary.json",
        "sf_spatial_od_summary.json",
        "sf_bpr_anchor_audit_summary.json",
        "sf_coalition_sensitivity_summary.json",
        "sf_coalition_partition_scan_summary.json",
        "sf_coalition_pooling_summary.json",
        "response_paper/figures/fig_braess_response_matrices_data.json",
    ):
        assert_finite(read(name), name)

    required_figures = (
        "response_paper/figures/fig_braess_response_matrices.pdf",
        "response_paper/figures/fig_sf_response_decomposition.pdf",
        "response_paper/figures/fig_sf_spatial_od.pdf",
        "response_paper/figures/fig_sf_alpha_screening.pdf",
        "response_paper/figures/fig_coalition_sensitivity.pdf",
        "response_paper/figures/fig_coalition_partition_scan.pdf",
    )
    for name in required_figures:
        path = ROOT / name
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing figure {name}")

    text_files = [ROOT / "paper_response_interoperability.tex"]
    text_files.extend((ROOT / "response_paper/chapters").glob("*.tex"))
    text_files.append(ROOT / "response_to_reviewers_2026-08-01.md")
    for path in text_files:
        text = path.read_text()
        if "[[" in text or "AFTER_MERGE" in text:
            raise AssertionError(f"unresolved result placeholder in {path}")
    print("revision package validation passed")


if __name__ == "__main__":
    main()
