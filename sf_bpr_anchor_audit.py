"""Pointwise BPR monotonicity audit for the reported Sioux Falls anchors."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sf_revision_common import (
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    profile_record,
    solve_benchmarks,
    write_json,
)


ALPHA = 0.9
CONFIGS = (
    ("8 equal", [1 / 8] * 8),
    ("4 equal", [1 / 4] * 4),
    ("2 equal", [1 / 2] * 2),
    ("skewed", [0.7, 0.1, 0.1, 0.1]),
    ("monopoly", [1.0]),
)
REGIMES = ("mass_linked_prior", "full_information")
RAW_PATH = Path("sf_bpr_anchor_audit_results.json")
SUMMARY_PATH = Path("sf_bpr_anchor_audit_summary.json")


def remap(old_keys, new_keys, sol, players):
    positions = {key: index for index, key in enumerate(old_keys)}
    f_h = np.zeros(len(new_keys))
    f_f = np.zeros((players, len(new_keys)))
    for new_index, key in enumerate(new_keys):
        old_index = positions.get(key)
        if old_index is not None:
            f_h[new_index] = sol["fH"][old_index]
            f_f[:, new_index] = sol["fF"][:, old_index]
    return f_h, f_f


def solve(sf, routes, oracle, shares, lam_beta, theta_bar, beta_override):
    graph, edge_idx = oracle
    game = SE.SFGame(
        sf, ALPHA, shares, lam_beta, routes, mode="eq",
        theta_bar=theta_bar, beta_override=beta_override,
    )
    warm = None
    for round_index in range(25):
        sol = game.solve(tol=1e-6, warm=warm, verbose=False)
        weights_h, _ = game.edge_weights(sol["fH"], sol["fF"])
        scale = max(float(np.max(game.A.T @ weights_h)), 1e-12)
        violations, worst = SE.certificate(
            game, sol, graph, edge_idx, 1e-6 * scale
        )
        if not violations:
            sol.update(
                certified=bool(sol["gap"] <= 1e-6),
                cg_rounds=round_index + 1,
                n_paths=routes.total(),
                max_reduced_cost=float(worst),
                route_cost_scale=scale,
            )
            return sol, game
        old_keys = routes.keys()
        if sum(routes.add(od, path) for od, path, _ in violations) == 0:
            break
        game.rebuild()
        warm = remap(old_keys, routes.keys(), sol, game.K)
    sol.update(
        certified=False,
        cg_rounds=25,
        n_paths=routes.total(),
        max_reduced_cost=float(worst),
        route_cost_scale=scale,
    )
    return sol, game


def summarize(rows):
    by_regime = {}
    for regime in REGIMES:
        selected = [row for row in rows if row["regime"] == regime]
        item = {
            "n": len(selected),
            "certified_n": sum(row["certified"] for row in selected),
            "profile_exact_pass_n": sum(
                row["profile_psd_certified"] for row in selected
            ),
            "max_exact_violation_ratio": max(
                row["max_exact_violation_ratio"] for row in selected
            ),
            "anchors": {
                row["name"]: {
                    "recovery": row["recovery"],
                    "profile_exact_pass": row["profile_psd_certified"],
                    "max_exact_violation_ratio": row["max_exact_violation_ratio"],
                    "profile_coarse_pass": row["profile_coarse_certified"],
                    "max_coarse_violation_ratio": row["max_coarse_violation_ratio"],
                }
                for row in selected
            },
        }
        coarse = [row for row in selected if row["coarse_bound_available"]]
        item["profile_coarse_pass_n"] = (
            sum(row["profile_coarse_certified"] for row in coarse) if coarse else None
        )
        item["max_coarse_violation_ratio"] = (
            max(row["max_coarse_violation_ratio"] for row in coarse)
            if coarse else None
        )
        by_regime[regime] = item
    return {"by_regime": by_regime}


def main():
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)
    ue, ue_game, so = solve_benchmarks(sf)
    theta_bar = 0.5 * SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    routes = SE.RouteSet(sf["od_list"])
    oracle = network_oracle()
    rows = []
    for regime in REGIMES:
        for name, shares in CONFIGS:
            full_information = regime == "full_information"
            sol, game = solve(
                sf, routes, oracle, shares, lam_beta,
                theta_bar=None if full_information else theta_bar,
                beta_override=1.0 if full_information else None,
            )
            record = profile_record(
                {"name": name, "shares": shares, "regime": regime},
                sol, game, ue["J"], so["J"], full_information,
            )
            rows.append(record)
            print(
                f"{regime:<18} {name:<10} recovery={record['recovery']:.4f} "
                f"exact={record['profile_psd_certified']} "
                f"coarse={record['profile_coarse_certified']} "
                f"cert={record['certified']}",
                flush=True,
            )
    raw = {
        "protocol": {
            "alpha": ALPHA,
            "positive_prior_kappa": 0.5,
            "regimes": list(REGIMES),
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": rows,
    }
    write_json(RAW_PATH, raw)
    summary = summarize(rows)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2))
    if any(not row["certified"] for row in rows):
        raise SystemExit("one or more anchor profiles failed equilibrium certification")


if __name__ == "__main__":
    main()
