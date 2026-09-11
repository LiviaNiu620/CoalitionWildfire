"""Certified Sioux Falls scan over company count and coalition degree."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sf_revision_common import LAM_TIMES_M, SE, build_sioux_falls, network_oracle, solve_benchmarks, solve_profile, write_json


RAW_PATH = Path("sf_coalition_degree_results.json")
SUMMARY_PATH = Path("sf_coalition_degree_summary.json")
ALPHAS = (0.5, 0.9)
GAMMAS = (0.0, 0.5, 1.0)
COMPANY_COUNTS = (2, 4, 8)


def partitions(n):
    blocks = {"singletons": tuple((i,) for i in range(n)),
              "grand": (tuple(range(n)),)}
    if n >= 4:
        blocks["balanced_pairs"] = tuple((i, i + 1) for i in range(0, n, 2))
    if n >= 8:
        blocks["balanced_quads"] = tuple(
            tuple(range(i, i + 4)) for i in range(0, n, 4)
        )
    return blocks


def objective_types(n):
    return np.array([0.5 if i % 2 == 0 else 1.5 for i in range(n)], dtype=float)


def company_hhi(n):
    return 1.0 / n


def coalition_hhi(partition):
    shares = np.array([len(block) / sum(map(len, partition)) for block in partition])
    return float(np.dot(shares, shares))


def record(alpha, gamma, n, name, partition, sol, game, j_ue, j_so):
    return {
        "alpha": float(alpha),
        "gamma": float(gamma),
        "n_companies": int(n),
        "n_coalitions": int(len(partition)),
        "partition": name,
        "coalitions": [list(block) for block in partition],
        "objective_types": objective_types(n).tolist(),
        "company_HHI": company_hhi(n),
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
    }


def summarize(raw):
    rows = raw["profiles"]
    key = {(r["alpha"], r["gamma"], r["n_companies"], r["partition"]): r for r in rows}
    degree = []
    for alpha in ALPHAS:
        for gamma in GAMMAS:
            for n in COMPANY_COUNTS:
                names = list(partitions(n))
                singleton = key[(alpha, gamma, n, "singletons")]
                grand = key[(alpha, gamma, n, "grand")]
                degree.append({
                    "alpha": alpha,
                    "gamma": gamma,
                    "n_companies": n,
                    "singleton_recovery": singleton["recovery"],
                    "grand_recovery": grand["recovery"],
                    "grand_minus_singleton_recovery_pp": 100.0 * (grand["recovery"] - singleton["recovery"]),
                    "grand_minus_singleton_J_over_J_UE": (grand["J"] - singleton["J"]) / raw["J_UE"],
                    "partitions": names,
                })
    composition = []
    for alpha in ALPHAS:
        for gamma in GAMMAS:
            for n in (4, 8):
                if "balanced_pairs" not in partitions(n):
                    continue
                pair = key[(alpha, gamma, n, "balanced_pairs")]
                singleton = key[(alpha, gamma, n, "singletons")]
                composition.append({
                    "alpha": alpha,
                    "gamma": gamma,
                    "n_companies": n,
                    "pairs_minus_singletons_recovery_pp": 100.0 * (pair["recovery"] - singleton["recovery"]),
                    "pairs_minus_singletons_J_over_J_UE": (pair["J"] - singleton["J"]) / raw["J_UE"],
                })
    return {
        "protocol": raw["protocol"],
        "rows": len(rows),
        "certified_rows": sum(r["certified"] for r in rows),
        "max_VI_gap": max(r["VI_gap"] for r in rows),
        "max_relative_omitted_route_slack": max(r["omitted_route_slack"] / r["route_cost_scale"] for r in rows),
        "degree_contrasts": degree,
        "pair_contrasts": composition,
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
    management_slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle()
    prior = {}
    if args.resume and RAW_PATH.exists():
        old = json.loads(RAW_PATH.read_text())
        prior = {(r["alpha"], r["gamma"], r["n_companies"], r["partition"]): r for r in old.get("profiles", [])}
    raw = {
        "protocol": {
            "governance": "identity-preserving coordination",
            "company_counts": list(COMPANY_COUNTS),
            "alphas": list(ALPHAS),
            "gammas": list(GAMMAS),
            "total_AV_mass_fixed": True,
            "objective_types": "alternating 0.5 and 1.5, mean one",
            "company_masses": "equal within n",
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }
    total_rows = sum(len(partitions(n)) for n in COMPANY_COUNTS) * len(ALPHAS) * len(GAMMAS)
    done = 0
    routes = SE.RouteSet(sf["od_list"])
    for alpha in ALPHAS:
        for n in COMPANY_COUNTS:
            shares = np.full(n, 1.0 / n)
            types = objective_types(n)
            for gamma in GAMMAS:
                for name, partition in partitions(n).items():
                    done += 1
                    k = (alpha, gamma, n, name)
                    if k in prior:
                        raw["profiles"].append(prior[k])
                        print(f"[{done:02d}/{total_rows}] reused {k}", flush=True)
                        continue
                    if name == "singletons" and gamma != 0.0:
                        base_key = (alpha, 0.0, n, name)
                        base = prior.get(base_key)
                        if base is None:
                            base = next((item for item in raw["profiles"]
                                         if (item["alpha"], item["gamma"],
                                             item["n_companies"], item["partition"])
                                         == base_key), None)
                        if base is not None:
                            derived = dict(base)
                            derived["gamma"] = float(gamma)
                            derived["derived_by_partition_invariance"] = True
                            raw["profiles"].append(derived)
                            write_json(RAW_PATH, raw)
                            print(f"[{done:02d}/{total_rows}] derived {k} from gamma=0 invariance", flush=True)
                            continue
                    sol, game = solve_profile(
                        sf, routes, oracle,
                        alpha, shares, lam_beta,
                        coalitions=partition,
                        objective_types=types,
                        coordination_gamma=gamma,
                        management_slope=management_slope,
                        verbose=args.verbose,
                    )
                    row = record(alpha, gamma, n, name, partition, sol, game, ue["J"], so["J"])
                    raw["profiles"].append(row)
                    write_json(RAW_PATH, raw)
                    print(f"[{done:02d}/{total_rows}] alpha={alpha:.1f} n={n} gamma={gamma:.1f} {name:<15} rho={row['recovery']:.5f} gap={row['VI_gap']:.1e} cert={row['certified']}", flush=True)
    raw["profiles"].sort(key=lambda r: (r["alpha"], r["n_companies"], r["gamma"], r["partition"]))
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if summary["certified_rows"] != summary["rows"]:
        raise SystemExit("one or more coalition-degree profiles failed certification")


if __name__ == "__main__":
    main()
