"""Targeted initialization audit for the positive-prior Sioux Falls scan."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import sf_equilibrium as SE
from sf_simplex_scan import (
    ALPHA, CG_TOL_REL, KAPPA, LAM_TIMES_M, VI_TOL,
)
from sf_simplex_scan import network_oracle, solve_benchmark, solve_on_shared_routes
from sioux_falls_loader import build_sioux_falls


RAW = Path("sf_simplex_scan_results.json")
OUTPUT = Path("sf_simplex_multistart_results.json")


def widest_inversion(rows, threshold=0.05):
    hhi = np.array([row["HHI"] for row in rows])
    recovery = np.array([row["recovery"] for row in rows])
    travel_time = np.array([row["J"] for row in rows])
    i, j = np.triu_indices(len(rows), 1)
    eligible = (np.abs(hhi[i] - hhi[j]) >= threshold) & (
        (hhi[i] - hhi[j]) * (recovery[i] - recovery[j]) < 0
    )
    candidates = np.flatnonzero(eligible)
    chosen = candidates[
        np.argmax(np.abs(travel_time[i[candidates]] - travel_time[j[candidates]]))
    ]
    return rows[i[chosen]], rows[j[chosen]]


def structured_initial(game, use_last):
    f_h = np.zeros(game.P)
    f_f = np.zeros((game.K, game.P))
    for oi in range(game.n_od):
        columns = np.flatnonzero(game.masks[oi])
        selected = columns[-1] if use_last else columns[0]
        f_h[selected] = game.dem_hdv[oi]
        f_f[:, selected] = game.dem_firm[oi]
    return f_h, f_f


def main():
    raw = json.loads(RAW.read_text())
    first, second = widest_inversion(raw["profiles"])
    cases = [
        {"id": "four-equal-anchor", "shares": [0.25] * 4},
        {"id": "skewed-anchor", "shares": [0.7, 0.1, 0.1, 0.1]},
        {"id": first["id"], "shares": first["shares"]},
        {"id": second["id"], "shares": second["shares"]},
    ]

    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)
    ue, ue_game = solve_benchmark(sf, lam_beta, "ue")
    theta_bar = KAPPA * SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    oracle = network_oracle(sf)
    routes = SE.RouteSet(sf["od_list"])

    design = sorted(raw["profiles"], key=lambda row: row["HHI"])
    anchors = [
        [0.25] * 4,
        [0.7, 0.1, 0.1, 0.1],
        design[0]["shares"],
        design[-1]["shares"],
    ]
    for shares in anchors:
        sol, _ = solve_on_shared_routes(
            sf, routes, shares, lam_beta, theta_bar, oracle
        )
        if not sol["certified"]:
            raise RuntimeError("failed to rebuild the certified route basis")

    graph, edge_idx = oracle
    output_cases = []
    for case in cases:
        baseline, baseline_game = solve_on_shared_routes(
            sf, routes, case["shares"], lam_beta, theta_bar, oracle
        )
        starts = []
        for label in ("random-1", "random-2", "random-3", "first-route", "last-route"):
            game = SE.SFGame(
                sf, ALPHA, case["shares"], lam_beta, routes, mode="eq",
                theta_bar=theta_bar,
            )
            if label.startswith("random"):
                initial = game.initial(int(label.split("-")[1]))
            else:
                initial = structured_initial(game, label == "last-route")
            sol = game.solve(tol=VI_TOL, warm=initial)
            weight_h, _ = game.edge_weights(sol["fH"], sol["fF"])
            scale = float(np.max(game.A.T @ weight_h))
            violations, worst = SE.certificate(
                game, sol, graph, edge_idx, CG_TOL_REL * scale
            )
            starts.append({
                "start": label,
                "VI_gap": float(sol["gap"]),
                "omitted_route_slack": float(worst),
                "route_cost_scale": scale,
                "certified": bool(sol["gap"] <= VI_TOL and not violations),
                "relative_max_edge_flow_deviation": float(
                    np.max(np.abs(sol["x"] - baseline["x"]))
                    / max(np.max(baseline["x"]), 1e-12)
                ),
                "relative_TSTT_deviation": float(
                    abs(sol["J"] - baseline["J"]) / baseline["J"]
                ),
            })
            print(
                f"{case['id']:<18} {label:<12} gap={sol['gap']:.2e} "
                f"cert={starts[-1]['certified']} "
                f"dx={starts[-1]['relative_max_edge_flow_deviation']:.2e}",
                flush=True,
            )
        output_cases.append({
            **case,
            "baseline_J": float(baseline["J"]),
            "baseline_gap": float(baseline["gap"]),
            "starts": starts,
            "all_certified": all(item["certified"] for item in starts),
            "max_relative_edge_flow_deviation": max(
                item["relative_max_edge_flow_deviation"] for item in starts
            ),
            "max_relative_TSTT_deviation": max(
                item["relative_TSTT_deviation"] for item in starts
            ),
        })

    output = {
        "selection": (
            "equal and skewed anchors plus the result-selected largest TSTT "
            "inversion with abs(Delta HHI) >= 0.05"
        ),
        "random_seeds": [1, 2, 3],
        "structured_starts": ["first-route", "last-route"],
        "paths": routes.total(),
        "cases": output_cases,
        "all_certified": all(case["all_certified"] for case in output_cases),
        "max_relative_edge_flow_deviation": max(
            case["max_relative_edge_flow_deviation"] for case in output_cases
        ),
        "max_relative_TSTT_deviation": max(
            case["max_relative_TSTT_deviation"] for case in output_cases
        ),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2), flush=True)


if __name__ == "__main__":
    main()
