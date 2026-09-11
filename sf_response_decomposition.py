"""Mechanism-separated Sioux Falls ownership scan for the major revision."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sf_revision_common import (
    K,
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    profile_record,
    share_design,
    solve_benchmarks,
    solve_profile,
    summarize_group,
    write_json,
)


ALPHA = 0.9
RAW_PATH = Path("sf_response_decomposition_results.json")
SUMMARY_PATH = Path("sf_response_decomposition_summary.json")
REGIME_ORDER = ("mass_linked_prior", "full_information", "common_half")


def chunk_path(index):
    return Path(f"sf_response_decomposition_chunk_{index:02d}.json")


def regime_parameters(name, theta_bar):
    if name == "mass_linked_prior":
        return dict(theta_bar=theta_bar, beta_override=None, pure_fidelity=False)
    if name == "full_information":
        return dict(theta_bar=None, beta_override=1.0, pure_fidelity=True)
    if name == "common_half":
        return dict(theta_bar=None, beta_override=0.5, pure_fidelity=True)
    raise ValueError(name)


def summarize(raw):
    return {
        "protocol": raw["protocol"],
        "design_profiles": 200,
        "rows": len(raw["profiles"]),
        "by_regime": {
            regime: summarize_group(
                [row for row in raw["profiles"] if row["regime"] == regime],
                raw["J_UE"],
            )
            for regime in REGIME_ORDER
        },
    }


def merge_chunks(count):
    chunks = [json.loads(chunk_path(i).read_text()) for i in range(count)]
    rows = [row for chunk in chunks for row in chunk["profiles"]]
    keys = [(row["id"], row["regime"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate profile/regime keys")
    expected = 200 * len(REGIME_ORDER)
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} rows, found {len(rows)}")
    first = chunks[0]
    for chunk in chunks[1:]:
        if not np.isclose(chunk["J_UE"], first["J_UE"], rtol=0, atol=1e-6):
            raise RuntimeError("inconsistent J_UE across chunks")
        if not np.isclose(chunk["J_SO"], first["J_SO"], rtol=0, atol=1e-6):
            raise RuntimeError("inconsistent J_SO across chunks")
    by_id = {}
    for row in rows:
        by_id.setdefault(row["id"], []).append(row)
    if set(map(len, by_id.values())) != {len(REGIME_ORDER)}:
        raise RuntimeError("each ownership profile must appear in all regimes")
    for profile_rows in by_id.values():
        reference = profile_rows[0]
        for row in profile_rows[1:]:
            if not np.allclose(row["shares"], reference["shares"], rtol=0, atol=1e-14):
                raise RuntimeError("regime rows do not share an ownership vector")
            if not np.isclose(row["HHI"], reference["HHI"], rtol=0, atol=1e-14):
                raise RuntimeError("regime rows do not share an HHI")
    protocol = dict(first["protocol"])
    protocol.pop("chunk_index", None)
    protocol.pop("chunks", None)
    protocol["parallel_chunks"] = count
    raw = {
        "protocol": protocol,
        "J_UE": first["J_UE"],
        "J_SO": first["J_SO"],
        "profiles": sorted(rows, key=lambda row: (row["regime"], row["id"])),
    }
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if any(value["certified_n"] != value["n"] for value in summary["by_regime"].values()):
        raise SystemExit("one or more decomposition profiles failed certification")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunks", type=int, default=1)
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--merge-chunks", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    if args.merge_chunks:
        merge_chunks(args.chunks)
        return
    if not 0 <= args.chunk_index < args.chunks:
        parser.error("invalid chunk index")

    output = RAW_PATH if args.chunks == 1 else chunk_path(args.chunk_index)
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)
    ue, ue_game, so = solve_benchmarks(sf)
    theta_bar = 0.5 * SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])
    design = sorted(share_design(), key=lambda row: row["HHI"])
    selected = [row for i, row in enumerate(design) if i % args.chunks == args.chunk_index]

    prior = {}
    if args.resume and output.exists():
        old = json.loads(output.read_text())
        prior = {(row["id"], row["regime"]): row for row in old.get("profiles", [])}

    raw = {
        "protocol": {
            "alpha": ALPHA,
            "share_seed": 20260731,
            "profiles": 200,
            "regimes": list(REGIME_ORDER),
            "mass_linked_kappa": 0.5,
            "constant_beta": 0.5,
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
            "chunks": args.chunks,
            "chunk_index": args.chunk_index,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }

    extremes = [design[0]["shares"], design[-1]["shares"]]
    anchors = [[0.25] * K, [0.7, 0.1, 0.1, 0.1], *extremes]
    print("Building a shared route basis across three response regimes", flush=True)
    for shares in anchors:
        for regime in REGIME_ORDER:
            params = regime_parameters(regime, theta_bar)
            sol, _ = solve_profile(
                sf, routes, oracle, ALPHA, shares, lam_beta,
                theta_bar=params["theta_bar"],
                beta_override=params["beta_override"],
                verbose=args.verbose,
            )
            if not sol["certified"]:
                raise RuntimeError("route-basis anchor failed")
    print(f"Shared basis contains {routes.total()} paths", flush=True)

    warm = {regime: None for regime in REGIME_ORDER}
    for index, spec in enumerate(selected, 1):
        for regime in REGIME_ORDER:
            key = (spec["id"], regime)
            if key in prior:
                raw["profiles"].append(prior[key])
                continue
            params = regime_parameters(regime, theta_bar)
            sol, game = solve_profile(
                sf, routes, oracle, ALPHA, spec["shares"], lam_beta,
                theta_bar=params["theta_bar"],
                beta_override=params["beta_override"],
                warm=warm[regime], verbose=args.verbose,
            )
            record = profile_record(
                {**spec, "regime": regime}, sol, game, ue["J"], so["J"],
                params["pure_fidelity"],
            )
            raw["profiles"].append(record)
            warm[regime] = (sol["fH"], sol["fF"])
            write_json(output, raw)
            print(
                f"[{index:03d}/{len(selected)}] {spec['id']} {regime:<18} "
                f"rho={record['recovery']:.4f} gap={record['VI_gap']:.1e} "
                f"PSD={record['profile_psd_certified']} cert={record['certified']}",
                flush=True,
            )

    raw["profiles"].sort(key=lambda row: (row["regime"], row["id"]))
    write_json(output, raw)
    if args.chunks == 1:
        summary = summarize(raw)
        write_json(SUMMARY_PATH, summary)
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
