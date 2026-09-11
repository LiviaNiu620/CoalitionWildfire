"""Certified Sioux Falls experiment for heterogeneous AV coalition objectives.

The experiment uses identity-preserving coordination: every company retains
its own OD-demand constraints, while companies in the same coalition jointly
internalize a convex quadratic management penalty.  It isolates objective
composition because the primary assortative and mixed partitions have equal
company masses, equal coalition sizes, identical OD obligations, and the same
company- and coalition-level HHI.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sf_revision_common import (
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    solve_benchmarks,
    solve_profile,
    write_json,
)


RAW_PATH = Path("sf_coalition_composition_results.json")
SUMMARY_PATH = Path("sf_coalition_composition_summary.json")
ALPHAS = (0.5, 0.9)
GAMMAS = (0.0, 0.25, 0.5, 0.75, 1.0)
OBJECTIVE_TYPES = np.array([0.5, 0.5, 1.5, 1.5])
SHARES = np.full(4, 0.25)
PARTITIONS = {
    "singletons": ((0,), (1,), (2,), (3,)),
    "assortative": ((0, 1), (2, 3)),
    "mixed": ((0, 2), (1, 3)),
    "grand": ((0, 1, 2, 3),),
}


def coalition_hhi(partition, masses):
    total = float(np.sum(masses))
    shares = [float(np.sum(masses[list(block)])) / total for block in partition]
    return float(np.dot(shares, shares))


def coalition_objective_values(sol, game):
    """Evaluate the declared convex continuation objective by coalition."""
    x = np.asarray(sol["x"])
    own = np.asarray(sol["fF"]) @ game.A.T
    values = []
    for block in game.coalitions:
        idx = np.asarray(block, dtype=int)
        coalition_flow = np.sum(own[idx], axis=0)
        external = np.maximum(x - coalition_flow, 0.0)
        physical = game.t0 * (
            coalition_flow
            + SE.BPR_B
            * ((external + coalition_flow) ** (SE.BPR_ETA + 1.0)
               - external ** (SE.BPR_ETA + 1.0))
            / ((SE.BPR_ETA + 1.0) * game.cap ** SE.BPR_ETA)
        )
        governance = game.coalition_governance_matrix(block)
        coupled = governance @ own[idx]
        penalty = 0.5 * float(np.sum(
            game.management_slope[None, :] * own[idx] * coupled
        ))
        values.append({
            "members": list(block),
            "physical_integral": float(np.sum(physical)),
            "management_penalty": penalty,
            "objective": float(np.sum(physical) + penalty),
        })
    return values


def record_profile(alpha, gamma, partition_name, partition, sol, game, j_ue, j_so):
    route_support = {
        "hdv": int(np.sum(sol["fH"] > 1e-9)),
        "companies": [int(np.sum(row > 1e-9)) for row in sol["fF"]],
    }
    return {
        "alpha": float(alpha),
        "gamma": float(gamma),
        "partition": partition_name,
        "coalitions": [list(block) for block in partition],
        "objective_types": OBJECTIVE_TYPES.tolist(),
        "company_HHI": float(np.dot(SHARES, SHARES)),
        "coalition_HHI": coalition_hhi(partition, game.masses),
        "J": float(sol["J"]),
        "recovery": float((j_ue - sol["J"]) / (j_ue - j_so)),
        "VI_gap": float(sol["gap"]),
        "omitted_route_slack": float(sol["max_reduced_cost"]),
        "route_cost_scale": float(sol["route_cost_scale"]),
        "certified": bool(sol["certified"]),
        "iterations": int(sol["iters"]),
        "cg_rounds": int(sol["cg_rounds"]),
        "paths": int(sol["n_paths"]),
        "edge_flow": np.asarray(sol["x"]).tolist(),
        "route_support": route_support,
        "coalition_objectives": coalition_objective_values(sol, game),
    }


def summarize(raw):
    rows = raw["profiles"]
    lookup = {
        (row["alpha"], row["gamma"], row["partition"]): row
        for row in rows
    }
    contrasts = []
    for alpha in ALPHAS:
        for gamma in GAMMAS:
            assortative = lookup[(alpha, gamma, "assortative")]
            mixed = lookup[(alpha, gamma, "mixed")]
            delta_j = assortative["J"] - mixed["J"]
            contrasts.append({
                "alpha": alpha,
                "gamma": gamma,
                "J_assortative": assortative["J"],
                "J_mixed": mixed["J"],
                "delta_J_assortative_minus_mixed": delta_j,
                "delta_J_over_J_UE": delta_j / raw["J_UE"],
                "delta_recovery_percentage_points":
                    100.0 * (assortative["recovery"] - mixed["recovery"]),
            })

    invariance = []
    for alpha in ALPHAS:
        base = lookup[(alpha, 0.0, "singletons")]
        base_x = np.asarray(base["edge_flow"])
        candidates = [lookup[(alpha, 0.0, name)] for name in PARTITIONS]
        invariance.append({
            "alpha": alpha,
            "max_abs_TSTT_difference": float(max(
                abs(row["J"] - base["J"]) for row in candidates
            )),
            "max_abs_edge_flow_difference": float(max(
                np.max(np.abs(np.asarray(row["edge_flow"]) - base_x))
                for row in candidates
            )),
            "max_relative_TSTT_difference": float(max(
                abs(row["J"] - base["J"]) / raw["J_UE"]
                for row in candidates
            )),
        })

    coordination = []
    for alpha in ALPHAS:
        for name in PARTITIONS:
            base = lookup[(alpha, 0.0, name)]
            full = lookup[(alpha, 1.0, name)]
            coordination.append({
                "alpha": alpha,
                "partition": name,
                "delta_J_gamma1_minus_gamma0": full["J"] - base["J"],
                "delta_J_over_J_UE": (full["J"] - base["J"]) / raw["J_UE"],
                "delta_recovery_percentage_points":
                    100.0 * (full["recovery"] - base["recovery"]),
            })

    return {
        "protocol": raw["protocol"],
        "rows": len(rows),
        "certified_rows": sum(row["certified"] for row in rows),
        "max_VI_gap": max(row["VI_gap"] for row in rows),
        "max_relative_omitted_route_slack": max(
            row["omitted_route_slack"] / row["route_cost_scale"] for row in rows
        ),
        "partition_invariance_gamma0": invariance,
        "assortative_minus_mixed": contrasts,
        "coordination_gamma1_minus_gamma0": coordination,
        "max_abs_composition_delta_J_over_J_UE": max(
            abs(row["delta_J_over_J_UE"]) for row in contrasts
        ),
    }


def deterministic_operator_check(sf, management_slope, lam_beta):
    """Check the theorem that gamma=0 makes partition labels inert."""
    routes = SE.RouteSet(sf["od_list"])
    games = [
        SE.SFGame(
            sf, 0.5, SHARES, lam_beta, routes, mode="eq",
            coalitions=partition, objective_types=OBJECTIVE_TYPES,
            coordination_gamma=0.0, management_slope=management_slope,
        )
        for partition in PARTITIONS.values()
    ]
    f_h, f_f = games[0].initial(seed=17)
    base_h, base_f, *_ = games[0].operator(f_h, f_f)
    for game in games[1:]:
        grad_h, grad_f, *_ = game.operator(f_h, f_f)
        if not np.allclose(grad_h, base_h, rtol=0, atol=1e-12):
            raise AssertionError("gamma=0 changed the HDV operator")
        if not np.allclose(grad_f, base_f, rtol=0, atol=1e-12):
            raise AssertionError("gamma=0 changed the company operator")


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
    deterministic_operator_check(sf, management_slope, lam_beta)

    prior = {}
    if args.resume and RAW_PATH.exists():
        old = json.loads(RAW_PATH.read_text())
        prior = {
            (row["alpha"], row["gamma"], row["partition"]): row
            for row in old.get("profiles", [])
        }

    raw = {
        "protocol": {
            "governance": "identity-preserving coordination",
            "objective_types": OBJECTIVE_TYPES.tolist(),
            "company_shares": SHARES.tolist(),
            "company_HHI": float(np.dot(SHARES, SHARES)),
            "alphas": list(ALPHAS),
            "gammas": list(GAMMAS),
            "partitions": {
                name: [list(block) for block in partition]
                for name, partition in PARTITIONS.items()
            },
            "management_slope_calibration": "BPR derivative at certified all-HDV UE",
            "management_slope": management_slope.tolist(),
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
            "pooling": False,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }

    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])
    warm_by_partition = {}
    gamma0_warm = {}
    total_profiles = len(ALPHAS) * len(GAMMAS) * len(PARTITIONS)
    completed = 0

    for alpha in ALPHAS:
        for gamma in GAMMAS:
            for name, partition in PARTITIONS.items():
                completed += 1
                key = (alpha, gamma, name)
                if key in prior:
                    raw["profiles"].append(prior[key])
                    print(f"[{completed:02d}/{total_profiles}] reused {key}", flush=True)
                    continue
                if gamma == 0.0:
                    warm = gamma0_warm.get(alpha)
                else:
                    warm = warm_by_partition.get((alpha, name))
                sol, game = solve_profile(
                    sf, routes, oracle, alpha, SHARES, lam_beta,
                    coalitions=partition,
                    objective_types=OBJECTIVE_TYPES,
                    coordination_gamma=gamma,
                    management_slope=management_slope,
                    warm=warm,
                    verbose=args.verbose,
                )
                row = record_profile(
                    alpha, gamma, name, partition, sol, game,
                    ue["J"], so["J"],
                )
                raw["profiles"].append(row)
                state = (sol["fH"], sol["fF"])
                warm_by_partition[(alpha, name)] = state
                if gamma == 0.0 and alpha not in gamma0_warm:
                    gamma0_warm[alpha] = state
                raw["profiles"].sort(
                    key=lambda item: (item["alpha"], item["gamma"], item["partition"])
                )
                write_json(RAW_PATH, raw)
                print(
                    f"[{completed:02d}/{total_profiles}] alpha={alpha:.1f} "
                    f"gamma={gamma:.2f} {name:<11} rho={row['recovery']:.5f} "
                    f"gap={row['VI_gap']:.1e} cert={row['certified']}",
                    flush=True,
                )

    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if summary["certified_rows"] != summary["rows"]:
        raise SystemExit("one or more coalition profiles failed certification")


if __name__ == "__main__":
    main()
