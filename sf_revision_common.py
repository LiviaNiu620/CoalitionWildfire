"""Shared certified-solver utilities for the 2026-08 major-revision experiments."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

import sf_equilibrium as SE
from sf_simplex_scan import (
    DIRICHLET_A,
    SHARE_SEED,
    average_ranks,
    kendall_tau_b,
    pairwise_summary,
    share_design,
    spearman_rho,
)
from sioux_falls_loader import DATA_DIR, build_sioux_falls, parse_net


VI_TOL = 1e-6
CG_TOL_REL = 1e-6
K = 4
LAM_TIMES_M = 6.78


def network_oracle():
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    return SE.build_graph(links)


def remap_warm(old_keys, new_keys, sol, players):
    pos = {key: i for i, key in enumerate(old_keys)}
    f_h = np.zeros(len(new_keys))
    f_f = np.zeros((players, len(new_keys)))
    for j, key in enumerate(new_keys):
        i = pos.get(key)
        if i is not None:
            f_h[j] = sol["fH"][i]
            f_f[:, j] = sol["fF"][:, i]
    return f_h, f_f


def solve_profile(sf, routes, oracle, alpha, shares, lam_beta,
                  theta_bar=None, beta_override=None,
                  dem_firm_override=None, coalitions=None,
                  objective_types=None, coordination_gamma=0.0,
                  management_slope=None, management_slope_mode="reference",
                  warm=None, max_rounds=25,
                  operational_pooling=False,
                  verbose=False):
    """Solve and certify one profile, expanding the shared route set."""
    graph, edge_idx = oracle
    game = SE.SFGame(
        sf,
        alpha,
        shares,
        lam_beta,
        routes,
        mode="eq",
        theta_bar=theta_bar,
        beta_override=beta_override,
        dem_firm_override=dem_firm_override,
        coalitions=coalitions,
        objective_types=objective_types,
        coordination_gamma=coordination_gamma,
        management_slope=management_slope,
        management_slope_mode=management_slope_mode,
        operational_pooling=operational_pooling,
    )
    player_count = len(shares)
    current_warm = (warm if warm is not None
                    and warm[1].shape == (player_count, game.P) else None)

    for rnd in range(max_rounds):
        sol = game.solve(tol=VI_TOL, warm=current_warm, verbose=False)
        weights_h, _ = game.edge_weights(sol["fH"], sol["fF"])
        scale = max(float(np.max(game.A.T @ weights_h)), 1e-12)
        violations, worst = SE.certificate(
            game, sol, graph, edge_idx, CG_TOL_REL * scale
        )
        if verbose:
            print(
                f"      round={rnd + 1} paths={routes.total()} "
                f"gap={sol['gap']:.2e} omitted={worst:+.2e} "
                f"violations={len(violations)}",
                flush=True,
            )
        if not violations and sol["gap"] <= VI_TOL:
            sol.update(
                certified=True,
                cg_rounds=rnd + 1,
                n_paths=routes.total(),
                max_reduced_cost=float(worst),
                route_cost_scale=scale,
            )
            return sol, game

        # A column-generation round can exhaust the 20,000-iteration VI
        # budget just above tolerance after the oracle has already found no
        # omitted route.  Continue from that iterate on the unchanged route
        # set rather than returning an uncertified anchor.  This preserves the
        # declared 1e-6 threshold; it only gives the solver additional warm
        # restarts to attain it.
        if not violations:
            current_warm = (sol["fH"], sol["fF"])
            continue

        old_keys = routes.keys()
        added = sum(routes.add(oi, path) for oi, path, _ in violations)
        if added == 0 and sol["gap"] > VI_TOL:
            current_warm = (sol["fH"], sol["fF"])
            continue
        if added == 0:
            sol.update(
                certified=False,
                cg_rounds=rnd + 1,
                n_paths=routes.total(),
                max_reduced_cost=float(worst),
                route_cost_scale=scale,
            )
            return sol, game
        game.rebuild()
        current_warm = remap_warm(old_keys, routes.keys(), sol, player_count)

    sol.update(
        certified=False,
        cg_rounds=max_rounds,
        n_paths=routes.total(),
        max_reduced_cost=float(worst),
        route_cost_scale=scale,
    )
    return sol, game


def solve_benchmarks(sf):
    """UE and SO depend on total OD demand, not the population partition."""
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game = SE.solve_with_column_generation(
        sf, 0.9, [1.0], lam_beta, mode="ue", tol=1e-7, verbose=False
    )
    so, _ = SE.solve_with_column_generation(
        sf, 0.9, [1.0], lam_beta, mode="so", tol=1e-7, verbose=False
    )
    if not ue["certified"] or not so["certified"]:
        raise RuntimeError("UE or SO benchmark failed certification")
    return ue, ue_game, so


def bpr_monotonicity_diagnostics(sol, game, pure_fidelity):
    """Return exact edgewise PSD and coarse BPR sufficient-condition checks."""
    x = np.asarray(sol["x"])
    own = np.asarray(sol["fF"]) @ game.A.T
    cp = SE.bpr_deriv(x, game.t0, game.cap)
    cpp = SE.bpr_second_deriv(x, game.t0, game.cap)
    theta = game.perceived_slope(cp, own)
    mask = (x > 1e-9) & (cp > 1e-14)
    if not np.any(mask):
        return {
            "positive_slope_edges": 0,
            "edge_psd_pass_fraction": 1.0,
            "profile_psd_certified": True,
            "max_exact_mu": 0.0,
            "max_exact_violation_ratio": 0.0,
            "coarse_bound_available": bool(pure_fidelity),
            "coarse_bound_pass_fraction": 1.0 if pure_fidelity else None,
            "profile_coarse_certified": True if pure_fidelity else None,
            "max_beta_sigma": 0.0 if pure_fidelity else None,
            "max_coarse_violation_ratio": 0.0 if pure_fidelity else None,
        }

    kappa = cpp[mask] / cp[mask]
    if pure_fidelity:
        mu = kappa * np.sqrt(
            np.sum(game.betas[:, None] * own[:, mask] ** 2, axis=0)
        )
    else:
        a = theta[:, mask] / cp[None, mask]
        weighted = (game.betas[:, None] * own[:, mask]) ** 2 / a
        mu = kappa * np.sqrt(np.sum(weighted, axis=0))

    exact_pass = mu <= 2.0 + 1e-12
    out = {
        "positive_slope_edges": int(mask.sum()),
        "edge_psd_pass_fraction": float(np.mean(exact_pass)),
        "profile_psd_certified": bool(np.all(exact_pass)),
        "max_exact_mu": float(np.max(mu)),
        "max_exact_violation_ratio": float(np.max(mu) / 2.0),
        "coarse_bound_available": bool(pure_fidelity),
    }
    if pure_fidelity:
        sigma = np.max(own[:, mask], axis=0) / x[mask]
        beta_sigma = float(np.max(game.betas)) * sigma
        threshold = 4.0 / (SE.BPR_ETA - 1.0) ** 2
        coarse_pass = beta_sigma <= threshold + 1e-12
        out.update(
            coarse_bound_pass_fraction=float(np.mean(coarse_pass)),
            profile_coarse_certified=bool(np.all(coarse_pass)),
            max_beta_sigma=float(np.max(beta_sigma)),
            max_coarse_violation_ratio=float(np.max(beta_sigma) / threshold),
        )
    else:
        out.update(
            coarse_bound_pass_fraction=None,
            profile_coarse_certified=None,
            max_beta_sigma=None,
            max_coarse_violation_ratio=None,
        )
    return out


def profile_record(spec, sol, game, j_ue, j_so, pure_fidelity):
    return {
        **spec,
        "J": float(sol["J"]),
        "recovery": float((j_ue - sol["J"]) / (j_ue - j_so)),
        "VI_gap": float(sol["gap"]),
        "omitted_route_slack": float(sol["max_reduced_cost"]),
        "route_cost_scale": float(sol["route_cost_scale"]),
        "certified": bool(sol["certified"]),
        "iterations": int(sol["iters"]),
        "cg_rounds": int(sol["cg_rounds"]),
        "paths": int(sol["n_paths"]),
        **bpr_monotonicity_diagnostics(sol, game, pure_fidelity),
    }


def summarize_group(rows, j_ue):
    certified = [row for row in rows if row["certified"]]
    hhi = np.array([row["HHI"] for row in certified])
    recovery = np.array([row["recovery"] for row in certified])
    travel = np.array([row["J"] for row in certified])
    edge_counts = np.array([row["positive_slope_edges"] for row in certified])
    edge_pass = np.array([row["edge_psd_pass_fraction"] for row in certified])
    total_edges = int(edge_counts.sum())
    passed_edges = float(np.dot(edge_counts, edge_pass))
    out = {
        "n": len(rows),
        "certified_n": len(certified),
        "spearman": spearman_rho(hhi, recovery),
        "kendall": kendall_tau_b(hhi, recovery),
        "pairwise": [
            pairwise_summary(hhi, recovery, travel, j_ue, threshold)
            for threshold in (0.0, 0.02, 0.05)
        ],
        "recovery_range": [float(recovery.min()), float(recovery.max())],
        "recovery_sd": float(np.std(recovery, ddof=1)),
        "profile_psd_pass_n": int(sum(row["profile_psd_certified"] for row in certified)),
        "profile_psd_pass_rate": float(np.mean(
            [row["profile_psd_certified"] for row in certified]
        )),
        "edge_psd_pass_rate": float(passed_edges / total_edges),
        "max_exact_violation_ratio": float(max(
            row["max_exact_violation_ratio"] for row in certified
        )),
    }
    coarse = [row for row in certified if row["coarse_bound_available"]]
    out["coarse_bound_available"] = bool(coarse)
    if coarse:
        coarse_edge_counts = np.array([row["positive_slope_edges"] for row in coarse])
        coarse_edge_pass = np.array([row["coarse_bound_pass_fraction"] for row in coarse])
        out.update(
            profile_coarse_pass_n=int(sum(
                row["profile_coarse_certified"] for row in coarse
            )),
            profile_coarse_pass_rate=float(np.mean([
                row["profile_coarse_certified"] for row in coarse
            ])),
            edge_coarse_pass_rate=float(
                np.dot(coarse_edge_counts, coarse_edge_pass) / coarse_edge_counts.sum()
            ),
            max_coarse_violation_ratio=float(max(
                row["max_coarse_violation_ratio"] for row in coarse
            )),
        )
    else:
        out.update(
            profile_coarse_pass_n=None,
            profile_coarse_pass_rate=None,
            edge_coarse_pass_rate=None,
            max_coarse_violation_ratio=None,
        )
    return out


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


__all__ = [
    "CG_TOL_REL",
    "DIRICHLET_A",
    "K",
    "LAM_TIMES_M",
    "SHARE_SEED",
    "VI_TOL",
    "SE",
    "build_sioux_falls",
    "network_oracle",
    "profile_record",
    "share_design",
    "solve_benchmarks",
    "solve_profile",
    "summarize_group",
    "write_json",
]
