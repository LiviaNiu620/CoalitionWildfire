"""Full-information ownership scans across fleet-controlled demand shares."""
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


ALPHAS = (0.1, 0.3, 0.5, 0.7, 0.9)
RAW_PATH = Path("sf_alpha_screening_results.json")
SUMMARY_PATH = Path("sf_alpha_screening_summary.json")
DECOMPOSITION_PATH = Path("sf_response_decomposition_results.json")


def chunk_path(index):
    return Path(f"sf_alpha_screening_chunk_{index:02d}.json")


def fixed_subset():
    rows = []
    for row in share_design():
        within_group = int(row["id"].rsplit("-", 1)[1])
        if within_group < 20:
            rows.append(row)
    if len(rows) != 100:
        raise RuntimeError(f"expected 100 profiles, found {len(rows)}")
    return sorted(rows, key=lambda row: row["HHI"])


def summarize(raw):
    return {
        "protocol": raw["protocol"],
        "rows": len(raw["profiles"]),
        "by_alpha": {
            f"{alpha:g}": summarize_group(
                [row for row in raw["profiles"] if row["alpha"] == alpha],
                raw["J_UE"],
            )
            for alpha in ALPHAS
        },
    }


def merge_chunks(count):
    chunks = [json.loads(chunk_path(i).read_text()) for i in range(count)]
    rows = [row for chunk in chunks for row in chunk["profiles"]]
    keys = [(row["id"], row["alpha"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate profile/alpha keys")
    expected = 100 * len(ALPHAS)
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} rows, found {len(rows)}")
    first = chunks[0]
    for chunk in chunks[1:]:
        if not np.isclose(chunk["J_UE"], first["J_UE"], rtol=0, atol=1e-6):
            raise RuntimeError("inconsistent J_UE across alpha chunks")
        if not np.isclose(chunk["J_SO"], first["J_SO"], rtol=0, atol=1e-6):
            raise RuntimeError("inconsistent J_SO across alpha chunks")
    by_id = {}
    for row in rows:
        by_id.setdefault(row["id"], []).append(row)
    if set(map(len, by_id.values())) != {len(ALPHAS)}:
        raise RuntimeError("each ownership profile must appear at every alpha")
    for profile_rows in by_id.values():
        reference = profile_rows[0]
        if any(not np.allclose(row["shares"], reference["shares"], rtol=0, atol=1e-14)
               for row in profile_rows[1:]):
            raise RuntimeError("alpha rows do not share an ownership vector")
    protocol = dict(first["protocol"])
    protocol.pop("chunk_index", None)
    protocol.pop("chunks", None)
    protocol["parallel_chunks"] = count
    raw = {
        "protocol": protocol,
        "J_UE": first["J_UE"],
        "J_SO": first["J_SO"],
        "profiles": sorted(rows, key=lambda row: (row["alpha"], row["id"])),
    }
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if any(value["certified_n"] != value["n"] for value in summary["by_alpha"].values()):
        raise SystemExit("one or more alpha profiles failed certification")


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
    ue, _, so = solve_benchmarks(sf)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])
    design = fixed_subset()
    selected = [row for i, row in enumerate(design) if i % args.chunks == args.chunk_index]

    prior = {}
    reuse_alpha_nine = False
    # Alpha=0.9 with beta=1 is exactly the full-information regime in the
    # mechanism-decomposition scan.  Reuse those certified rows instead of
    # solving the same profile a second time.
    if DECOMPOSITION_PATH.exists():
        decomposition = json.loads(DECOMPOSITION_PATH.read_text())
        subset_ids = {row["id"] for row in design}
        for row in decomposition.get("profiles", []):
            if row.get("regime") != "full_information" or row["id"] not in subset_ids:
                continue
            reused = dict(row)
            reused.pop("regime", None)
            reused["alpha"] = 0.9
            reused["source"] = "reused certified full-information decomposition row"
            prior[(reused["id"], 0.9)] = reused
        reuse_alpha_nine = sum(alpha == 0.9 for _, alpha in prior) == len(design)
        if not reuse_alpha_nine:
            prior = {key: value for key, value in prior.items() if key[1] != 0.9}
    if args.resume and output.exists():
        old = json.loads(output.read_text())
        prior.update({
            (row["id"], row["alpha"]): row for row in old.get("profiles", [])
        })

    raw = {
        "protocol": {
            "alphas": list(ALPHAS),
            "profiles_per_alpha": 100,
            "subset_rule": "first 20 seeded draws within each Dirichlet concentration",
            "beta_override": 1.0,
            "alpha_0.9_source": (
                "certified full-information decomposition rows"
                if reuse_alpha_nine else "direct solve"
            ),
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
            "chunks": args.chunks,
            "chunk_index": args.chunk_index,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }

    anchors = [[0.25] * K, [0.7, 0.1, 0.1, 0.1], design[0]["shares"], design[-1]["shares"]]
    print("Building a shared route basis across alpha values", flush=True)
    anchor_alphas = ALPHAS[:-1] if reuse_alpha_nine else ALPHAS
    for alpha in anchor_alphas:
        lam_beta = LAM_TIMES_M / (alpha * total)
        for anchor_index, shares in enumerate(anchors):
            sol, _ = solve_profile(
                sf, routes, oracle, alpha, shares, lam_beta,
                beta_override=1.0, verbose=args.verbose,
            )
            if not sol["certified"]:
                raise RuntimeError(
                    "alpha route-basis anchor failed: "
                    f"alpha={alpha:g}, anchor={anchor_index}, "
                    f"gap={sol['gap']:.3e}, "
                    f"relative_omitted={sol['max_reduced_cost'] / sol['route_cost_scale']:.3e}, "
                    f"cg_rounds={sol['cg_rounds']}"
                )
    print(f"Shared basis contains {routes.total()} paths", flush=True)

    warm = {alpha: None for alpha in ALPHAS}
    for index, spec in enumerate(selected, 1):
        for alpha in ALPHAS:
            key = (spec["id"], alpha)
            if key in prior:
                raw["profiles"].append(prior[key])
                continue
            lam_beta = LAM_TIMES_M / (alpha * total)
            sol, game = solve_profile(
                sf, routes, oracle, alpha, spec["shares"], lam_beta,
                beta_override=1.0, warm=warm[alpha], verbose=args.verbose,
            )
            record = profile_record(
                {**spec, "alpha": alpha}, sol, game, ue["J"], so["J"], True
            )
            raw["profiles"].append(record)
            warm[alpha] = (sol["fH"], sol["fF"])
            write_json(output, raw)
            print(
                f"[{index:03d}/{len(selected)}] {spec['id']} alpha={alpha:.1f} "
                f"rho={record['recovery']:.4f} PSD={record['profile_psd_certified']} "
                f"cert={record['certified']}",
                flush=True,
            )

    raw["profiles"].sort(key=lambda row: (row["alpha"], row["id"]))
    write_json(output, raw)
    if args.chunks == 1:
        summary = summarize(raw)
        write_json(SUMMARY_PATH, summary)
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
