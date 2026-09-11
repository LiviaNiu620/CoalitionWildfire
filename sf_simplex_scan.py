"""Systematic four-fleet ownership-simplex scan on Sioux Falls.

The design is locked in ``plan/experiment-protocol.md``. The scan keeps the
number of fleet decision centers fixed at four, draws 200 ownership vectors
from five predeclared symmetric Dirichlet distributions, and uses a positive
prior slope. Every retained profile is checked against the full network with
player-specific shortest paths.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

import sf_equilibrium as SE
from sioux_falls_loader import DATA_DIR, build_sioux_falls, parse_net


ALPHA = 0.9
K = 4
LAM_TIMES_M = 6.78
KAPPA = 0.5
DIRICHLET_A = (0.2, 0.5, 1.0, 2.0, 5.0)
SAMPLES_PER_A = 40
SHARE_SEED = 20260731
BOOTSTRAP_SEED = 20260732
BOOTSTRAP_REPS = 10_000
VI_TOL = 1e-6
CG_TOL_REL = 1e-6
RAW_PATH = Path("sf_simplex_scan_results.json")
SUMMARY_PATH = Path("sf_simplex_scan_summary.json")


def chunk_raw_path(index):
    return Path(f"sf_simplex_scan_chunk_{index:02d}.json")


def chunk_summary_path(index):
    return Path(f"sf_simplex_scan_chunk_{index:02d}_summary.json")


def average_ranks(values):
    """One-based average ranks with exact-tie handling."""
    values = np.asarray(values)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + 1 + stop)
        start = stop
    return ranks


def spearman_rho(x, y):
    rx = average_ranks(x)
    ry = average_ranks(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def kendall_tau_b(x, y):
    """Kendall tau-b from all unordered pairs, including tie correction."""
    i, j = np.triu_indices(len(x), 1)
    sx = np.sign(np.asarray(x)[i] - np.asarray(x)[j])
    sy = np.sign(np.asarray(y)[i] - np.asarray(y)[j])
    concordant = np.sum(sx * sy > 0)
    discordant = np.sum(sx * sy < 0)
    ties_x_only = np.sum((sx == 0) & (sy != 0))
    ties_y_only = np.sum((sy == 0) & (sx != 0))
    denominator = np.sqrt(
        (concordant + discordant + ties_x_only)
        * (concordant + discordant + ties_y_only)
    )
    if denominator == 0:
        return np.nan
    return float((concordant - discordant) / denominator)


def share_design(samples_per_a: int = SAMPLES_PER_A) -> list[dict]:
    rng = np.random.default_rng(SHARE_SEED)
    design = []
    for a in DIRICHLET_A:
        draws = rng.dirichlet(np.full(K, a), size=samples_per_a)
        for j, shares in enumerate(draws):
            design.append({
                "id": f"a{a:g}-{j:03d}",
                "dirichlet_a": a,
                "shares": shares.tolist(),
                "HHI": float(np.dot(shares, shares)),
            })
    return design


def network_oracle(sf):
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


def solve_on_shared_routes(sf, routes, shares, lam_beta, theta_bar, oracle,
                           warm=None, max_rounds=25, verbose=False):
    """Solve one profile while expanding the shared route set as needed."""
    graph, edge_idx = oracle
    game = SE.SFGame(sf, ALPHA, shares, lam_beta, routes, mode="eq",
                     theta_bar=theta_bar)
    if warm is not None and warm[1].shape == (K, game.P):
        current_warm = warm
    else:
        current_warm = None

    for rnd in range(max_rounds):
        sol = game.solve(tol=VI_TOL, warm=current_warm, verbose=False)
        weights_h, _ = game.edge_weights(sol["fH"], sol["fF"])
        scale = float(np.max(game.A.T @ weights_h))
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
        if not violations:
            sol.update(
                certified=bool(sol["gap"] <= VI_TOL),
                cg_rounds=rnd + 1,
                n_paths=routes.total(),
                max_reduced_cost=float(worst),
                route_cost_scale=scale,
            )
            return sol, game

        old_keys = routes.keys()
        added = sum(routes.add(oi, path) for oi, path, _ in violations)
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
        current_warm = remap_warm(old_keys, routes.keys(), sol, K)

    sol.update(
        certified=False,
        cg_rounds=max_rounds,
        n_paths=routes.total(),
        max_reduced_cost=float(worst),
        route_cost_scale=scale,
    )
    return sol, game


def solve_benchmark(sf, lam_beta, mode):
    sol, game = SE.solve_with_column_generation(
        sf, ALPHA, [1.0], lam_beta, mode=mode, tol=1e-7, verbose=False
    )
    if not sol["certified"]:
        raise RuntimeError(f"{mode} benchmark failed certification")
    return sol, game


def profile_record(spec, sol, game, j_ue, j_so):
    cp = SE.bpr_deriv(sol["x"], game.t0, game.cap)
    own_edge = sol["fF"] @ game.A.T
    theta = game.perceived_slope(cp, own_edge)
    pos = cp > 1e-12
    effective = theta[:, pos] / cp[None, pos]
    beta = game.betas
    return {
        **spec,
        "Sigma_coverage": float(np.sum(1.0 / beta)),
        "beta_min": float(beta.min()),
        "beta_max": float(beta.max()),
        "beta_effective_min": float(effective.min()),
        "beta_effective_max": float(effective.max()),
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


def pairwise_summary(hhi, recovery, travel_time, j_ue, min_delta=0.0):
    i, j = np.triu_indices(len(hhi), 1)
    dh = hhi[i] - hhi[j]
    dr = recovery[i] - recovery[j]
    eligible = (np.abs(dh) > 1e-12) & (np.abs(dh) >= min_delta)
    inverted = eligible & (dh * dr < 0)
    differences = np.abs(travel_time[i] - travel_time[j])[inverted]
    eligible_count = int(eligible.sum())
    return {
        "min_abs_delta_HHI": min_delta,
        "eligible_pairs": eligible_count,
        "inverted_pairs": int(inverted.sum()),
        "inversion_rate": (
            float(inverted.sum() / eligible_count) if eligible_count else None
        ),
        "median_abs_delta_J": float(np.median(differences)) if differences.size else None,
        "max_abs_delta_J": float(np.max(differences)) if differences.size else None,
        "median_abs_delta_J_pct_UE": (
            float(100 * np.median(differences) / j_ue) if differences.size else None
        ),
        "max_abs_delta_J_pct_UE": (
            float(100 * np.max(differences) / j_ue) if differences.size else None
        ),
        "share_inversions_over_0.10pct_UE": (
            float(np.mean(differences >= 0.001 * j_ue)) if differences.size else None
        ),
        "share_inversions_over_0.25pct_UE": (
            float(np.mean(differences >= 0.0025 * j_ue)) if differences.size else None
        ),
        "share_inversions_over_0.50pct_UE": (
            float(np.mean(differences >= 0.005 * j_ue)) if differences.size else None
        ),
    }


def bootstrap_rank(hhi, recovery, reps=BOOTSTRAP_REPS):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n = len(hhi)
    values = np.empty((reps, 3))
    for b in range(reps):
        idx = rng.integers(0, n, n)
        hb = hhi[idx]
        rb = recovery[idx]
        values[b, 0] = spearman_rho(hb, rb)
        tau = kendall_tau_b(hb, rb)
        values[b, 1] = tau
        ii, jj = np.triu_indices(n, 1)
        dh = hb[ii] - hb[jj]
        dr = rb[ii] - rb[jj]
        valid = np.abs(dh) > 1e-12
        values[b, 2] = np.mean((dh[valid] * dr[valid]) < 0)
    quantiles = np.quantile(values, (0.025, 0.975), axis=0)
    return {
        "replicates": reps,
        "seed": BOOTSTRAP_SEED,
        "spearman_95pct": quantiles[:, 0].tolist(),
        "kendall_95pct": quantiles[:, 1].tolist(),
        "inversion_rate_95pct": quantiles[:, 2].tolist(),
        "interpretation": "design-set resampling, not market-population uncertainty",
    }


def analyze(raw, bootstrap_reps=BOOTSTRAP_REPS):
    rows = [row for row in raw["profiles"] if row["certified"]]
    hhi = np.array([row["HHI"] for row in rows])
    recovery = np.array([row["recovery"] for row in rows])
    travel_time = np.array([row["J"] for row in rows])
    spearman = spearman_rho(hhi, recovery)
    kendall = kendall_tau_b(hhi, recovery)
    by_a = {}
    for a in DIRICHLET_A:
        group = [row for row in rows if row["dirichlet_a"] == a]
        gh = np.array([row["HHI"] for row in group])
        gr = np.array([row["recovery"] for row in group])
        gj = np.array([row["J"] for row in group])
        by_a[str(a)] = {
            "n": len(group),
            "HHI_range": [float(gh.min()), float(gh.max())],
            "recovery_range": [float(gr.min()), float(gr.max())],
            "spearman": spearman_rho(gh, gr) if len(group) > 1 else None,
            "kendall": kendall_tau_b(gh, gr) if len(group) > 1 else None,
            "pairwise": pairwise_summary(gh, gr, gj, raw["J_UE"]),
        }
    summary = {
        "design_n": len(raw["profiles"]),
        "certified_n": len(rows),
        "failed_n": len(raw["profiles"]) - len(rows),
        "spearman": spearman,
        "kendall": kendall,
        "pairwise": [
            pairwise_summary(hhi, recovery, travel_time, raw["J_UE"], threshold)
            for threshold in (0.0, 0.02, 0.05)
        ],
        "bootstrap": bootstrap_rank(hhi, recovery, bootstrap_reps),
        "by_dirichlet_concentration": by_a,
    }
    return summary


def merge_chunks(count, bootstrap_reps):
    chunks = [json.loads(chunk_raw_path(index).read_text()) for index in range(count)]
    first = chunks[0]
    profiles = [row for chunk in chunks for row in chunk["profiles"]]
    ids = [row["id"] for row in profiles]
    if len(ids) != len(set(ids)):
        raise RuntimeError("duplicate profile IDs across scan chunks")
    expected = SAMPLES_PER_A * len(DIRICHLET_A)
    if len(profiles) != expected:
        raise RuntimeError(f"expected {expected} profiles, found {len(profiles)}")
    for chunk in chunks[1:]:
        for key in ("J_UE", "J_SO", "recoverable_gap"):
            if not np.isclose(chunk[key], first[key], rtol=0, atol=1e-6):
                raise RuntimeError(f"inconsistent {key} across chunks")
    protocol = dict(first["protocol"])
    protocol.pop("chunk_index", None)
    protocol.pop("chunks", None)
    protocol["parallel_chunks"] = count
    raw = {
        "protocol": protocol,
        "J_UE": first["J_UE"],
        "J_SO": first["J_SO"],
        "recoverable_gap": first["recoverable_gap"],
        "profiles": sorted(profiles, key=lambda row: row["id"]),
    }
    RAW_PATH.write_text(json.dumps(raw, indent=2) + "\n")
    summary = analyze(raw, bootstrap_reps)
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    if summary["failed_n"]:
        raise SystemExit(f"{summary['failed_n']} profiles failed certification")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-per-a", type=int, default=SAMPLES_PER_A)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--chunks", type=int, default=1)
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--merge-chunks", action="store_true")
    args = parser.parse_args()

    if args.merge_chunks:
        merge_chunks(args.chunks, args.bootstrap_reps)
        return
    if not 0 <= args.chunk_index < args.chunks:
        parser.error("--chunk-index must be between zero and --chunks minus one")

    raw_path = RAW_PATH if args.chunks == 1 else chunk_raw_path(args.chunk_index)
    summary_path = (
        SUMMARY_PATH if args.chunks == 1 else chunk_summary_path(args.chunk_index)
    )

    design = share_design(args.samples_per_a)
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)
    ue, ue_game = solve_benchmark(sf, lam_beta, "ue")
    so, _ = solve_benchmark(sf, lam_beta, "so")
    theta_bar = KAPPA * SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle(sf)
    routes = SE.RouteSet(sf["od_list"])

    existing = {}
    if args.resume:
        resume_paths = [path for path in (raw_path, RAW_PATH) if path.exists()]
        for path in resume_paths:
            prior = json.loads(path.read_text())
            existing.update({row["id"]: row for row in prior.get("profiles", [])})

    # Build a shared route basis from predeclared anchors before the scan.
    extremes = sorted(design, key=lambda row: row["HHI"])
    anchor_specs = [
        [0.25] * 4,
        [0.7, 0.1, 0.1, 0.1],
        extremes[0]["shares"],
        extremes[-1]["shares"],
    ]
    warm = None
    print("Building the shared certified route basis", flush=True)
    for shares in anchor_specs:
        sol, _ = solve_on_shared_routes(
            sf, routes, shares, lam_beta, theta_bar, oracle, warm=warm,
            verbose=args.verbose,
        )
        if not sol["certified"]:
            raise RuntimeError("anchor route-basis solve failed certification")
        warm = (sol["fH"], sol["fF"])
    print(f"Shared basis contains {routes.total()} paths", flush=True)

    profiles = []
    warm = None
    raw = {
        "protocol": {
            "alpha": ALPHA,
            "K": K,
            "lambda_times_M": LAM_TIMES_M,
            "kappa": KAPPA,
            "dirichlet_a": list(DIRICHLET_A),
            "samples_per_a": args.samples_per_a,
            "share_seed": SHARE_SEED,
            "VI_tolerance": VI_TOL,
            "omitted_route_relative_tolerance": CG_TOL_REL,
            "chunks": args.chunks,
            "chunk_index": args.chunk_index,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "recoverable_gap": float(ue["J"] - so["J"]),
        "profiles": profiles,
    }
    full_design = sorted(design, key=lambda row: row["HHI"])
    ordered_design = [
        row for index, row in enumerate(full_design)
        if index % args.chunks == args.chunk_index
    ]
    design_ids = {row["id"] for row in ordered_design}
    existing = {key: value for key, value in existing.items() if key in design_ids}
    for index, spec in enumerate(ordered_design, 1):
        if spec["id"] in existing:
            profiles.append(existing[spec["id"]])
            continue
        sol, game = solve_on_shared_routes(
            sf, routes, spec["shares"], lam_beta, theta_bar, oracle,
            warm=warm, verbose=args.verbose,
        )
        record = profile_record(spec, sol, game, ue["J"], so["J"])
        profiles.append(record)
        warm = (sol["fH"], sol["fF"])
        raw["profiles"] = profiles
        raw_path.write_text(json.dumps(raw, indent=2) + "\n")
        print(
            f"[{index:03d}/{len(ordered_design)}] {spec['id']:<10} "
            f"HHI={spec['HHI']:.3f} recovery={record['recovery']:.3f} "
            f"gap={record['VI_gap']:.1e} cert={record['certified']} "
            f"paths={record['paths']}",
            flush=True,
        )

    profiles.sort(key=lambda row: row["id"])
    raw["profiles"] = profiles
    raw_path.write_text(json.dumps(raw, indent=2) + "\n")
    summary = analyze(raw, args.bootstrap_reps)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    if summary["failed_n"]:
        raise SystemExit(f"{summary['failed_n']} profiles failed certification")


if __name__ == "__main__":
    main()
