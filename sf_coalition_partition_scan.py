"""Certified Sioux Falls scan over all set partitions of four companies."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coalition_private_costs import company_accounting_costs
from sf_revision_common import (
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    solve_benchmarks,
    solve_profile,
    write_json,
)
from sf_simplex_scan import kendall_tau_b, spearman_rho


RAW_PATH = Path("sf_coalition_partition_scan_results.json")
SUMMARY_PATH = Path("sf_coalition_partition_scan_summary.json")
ALPHAS = (0.5, 0.9)
GAMMAS = (0.5, 1.0)
N = 4
SHARES = np.full(N, 0.25)
OBJECTIVE_TYPES = np.array([0.5, 0.5, 1.5, 1.5])


def all_partitions(n):
    blocks = []
    def extend(i):
        if i == n:
            yield tuple(tuple(block) for block in blocks)
            return
        for j in range(len(blocks)):
            blocks[j].append(i)
            yield from extend(i + 1)
            blocks[j].pop()
        blocks.append([i])
        yield from extend(i + 1)
        blocks.pop()
    parts = list(extend(0))
    return sorted(parts, key=lambda p: (len(p), tuple(tuple(b) for b in p)))


PARTITIONS = all_partitions(N)
PARTITION_IDS = {
    partition: f"P{index:02d}" for index, partition in enumerate(PARTITIONS, start=1)
}


def type_class(partition):
    signatures = []
    for block in partition:
        low = sum(OBJECTIVE_TYPES[i] < 1.0 for i in block)
        high = sum(OBJECTIVE_TYPES[i] > 1.0 for i in block)
        signatures.append((low, high))
    signatures.sort()
    if signatures == [(0, 2), (2, 0)]:
        return "assortative"
    if signatures == [(1, 1), (1, 1)]:
        return "mixed_pairs"
    if signatures == [(0, 1), (1, 0), (1, 0), (1, 1)]:
        return "mixed_singletons"
    return "other"


def coalition_hhi(partition):
    shares = np.array([len(block) / N for block in partition], dtype=float)
    return float(np.dot(shares, shares))


def record(alpha, gamma, partition, sol, game, j_ue, j_so):
    return {
        "alpha": float(alpha),
        "gamma": float(gamma),
        "partition_id": PARTITION_IDS[partition],
        "partition": [list(block) for block in partition],
        "n_companies": N,
        "n_coalitions": len(partition),
        "objective_types": OBJECTIVE_TYPES.tolist(),
        "type_class": type_class(partition),
        "company_HHI": 0.25,
        "coalition_HHI": coalition_hhi(partition),
        "J": float(sol["J"]),
        "recovery": float((j_ue - sol["J"]) / (j_ue - j_so)),
        "VI_gap": float(sol["gap"]),
        "omitted_route_slack": float(sol["max_reduced_cost"]),
        "route_cost_scale": float(sol["route_cost_scale"]),
        "certified": bool(sol["certified"]),
        "iterations": int(sol["iters"]),
        "cg_rounds": int(sol["cg_rounds"]),
        "paths": int(sol["n_paths"]),
        **company_accounting_costs(sol, game),
    }


def pairwise_inversion_rate(hhi, recovery):
    total = inversions = 0
    for i in range(len(hhi)):
        for j in range(i + 1, len(hhi)):
            dh = hhi[i] - hhi[j]
            dr = recovery[i] - recovery[j]
            if abs(dh) < 1e-12:
                continue
            total += 1
            if dh * dr < 0:
                inversions += 1
    return float(inversions / total) if total else 0.0


def summarize(raw):
    groups = {}
    for row in raw["profiles"]:
        groups.setdefault((row["alpha"], row["gamma"]), []).append(row)
    out_groups = []
    for (alpha, gamma), rows in sorted(groups.items()):
        hhi = np.array([r["coalition_HHI"] for r in rows])
        recovery = np.array([r["recovery"] for r in rows])
        out_groups.append({
            "alpha": alpha,
            "gamma": gamma,
            "rows": len(rows),
            "spearman_coalition_HHI_recovery": spearman_rho(hhi, recovery),
            "kendall_coalition_HHI_recovery": kendall_tau_b(hhi, recovery),
            "pairwise_inversion_rate": pairwise_inversion_rate(hhi, recovery),
            "recovery_range": [float(recovery.min()), float(recovery.max())],
            "recovery_sd": float(np.std(recovery, ddof=1)),
            "max_J_range_over_J_UE": float((max(r["J"] for r in rows)
                                              - min(r["J"] for r in rows)) / raw["J_UE"]),
            "class_means": {
                cls: float(np.mean([r["recovery"] for r in rows if r["type_class"] == cls]))
                for cls in sorted({r["type_class"] for r in rows})
            },
        })
    return {
        "protocol": raw["protocol"],
        "partition_count": len(PARTITIONS),
        "rows": len(raw["profiles"]),
        "certified_rows": sum(bool(r["certified"]) for r in raw["profiles"]),
        "max_VI_gap": max(r["VI_gap"] for r in raw["profiles"]),
        "max_relative_omitted_route_slack": max(
            r["omitted_route_slack"] / r["route_cost_scale"] for r in raw["profiles"]),
        "by_group": out_groups,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game, so = solve_benchmarks(sf)
    slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])
    prior = {}
    if args.resume and RAW_PATH.exists():
        old = json.loads(RAW_PATH.read_text())
        prior = {(r["alpha"], r["gamma"], r["partition_id"]): r
                 for r in old.get("profiles", [])}
    raw = {
        "protocol": {
            "governance": "identity-preserving coordination",
            "n_companies": N,
            "partition_count": len(PARTITIONS),
            "alphas": list(ALPHAS),
            "gammas": list(GAMMAS),
            "partition_ids": {PARTITION_IDS[p]: [list(b) for b in p]
                              for p in PARTITIONS},
            "objective_types": OBJECTIVE_TYPES.tolist(),
            "total_AV_mass_fixed": True,
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }
    total_rows = len(ALPHAS) * len(GAMMAS) * len(PARTITIONS)
    done = 0
    for alpha in ALPHAS:
        warm = {}
        for gamma in GAMMAS:
            for partition in PARTITIONS:
                done += 1
                pid = PARTITION_IDS[partition]
                key = (alpha, gamma, pid)
                if key in prior:
                    raw["profiles"].append(prior[key])
                    print(f"[{done:02d}/{total_rows}] reused {key}", flush=True)
                    continue
                sol, game = solve_profile(
                    sf, routes, oracle, alpha, SHARES, lam_beta,
                    coalitions=partition, objective_types=OBJECTIVE_TYPES,
                    coordination_gamma=gamma, management_slope=slope,
                    warm=warm.get(pid), verbose=args.verbose,
                )
                row = record(
                    alpha, gamma, partition, sol, game, ue["J"], so["J"]
                )
                raw["profiles"].append(row)
                warm[pid] = (sol["fH"], sol["fF"])
                write_json(RAW_PATH, raw)
                print(f"[{done:02d}/{total_rows}] alpha={alpha:.1f} gamma={gamma:.1f} "
                      f"{pid} K={len(partition)} rho={row['recovery']:.5f} "
                      f"gap={row['VI_gap']:.1e} cert={row['certified']}", flush=True)
    raw["profiles"].sort(key=lambda r: (r["alpha"], r["gamma"], r["partition_id"]))
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if summary["certified_rows"] != summary["rows"]:
        raise SystemExit("one or more partition-scan profiles failed certification")


if __name__ == "__main__":
    main()
