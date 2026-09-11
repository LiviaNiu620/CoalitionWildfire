"""Cross-topology validation of coalition regimes on small BPR networks."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import networkx as nx
import numpy as np

import sf_equilibrium as SE
from coalition_private_costs import company_accounting_costs
from sf_revision_common import LAM_TIMES_M, solve_profile


OUT = Path("cross_network_regime_results.json")
SUMMARY = Path("cross_network_regime_summary.json")
ALPHAS = (0.3, 0.5, 0.7, 0.9)
GAMMAS = (0.25, 0.5, 0.75, 1.0)
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
SHARES = np.full(4, 0.25)


def partitions(n=4):
    blocks = []
    def rec(i):
        if i == n:
            yield tuple(tuple(b) for b in blocks)
            return
        for j in range(len(blocks)):
            blocks[j].append(i); yield from rec(i + 1); blocks[j].pop()
        blocks.append([i]); yield from rec(i + 1); blocks.pop()
    return sorted(rec(0), key=lambda p: (len(p), p))


PARTITIONS = partitions()


def build_network(name):
    if name == "braess":
        edges = [
            (0, 1, 1900.0, 0.8), (1, 3, 100000.0, 9.0),
            (0, 2, 100000.0, 9.0), (2, 3, 1900.0, 0.8),
            (1, 2, 100000.0, 0.5),
        ]
        ods = [(0, 3, 5000.0)]
    elif name == "grid":
        edges = [
            (0, 1, 2600.0, 2.0), (1, 2, 1700.0, 1.0),
            (3, 4, 1800.0, 1.0), (4, 5, 2800.0, 2.0),
            (6, 7, 2400.0, 2.2), (7, 8, 1900.0, 1.1),
            (0, 3, 1800.0, 1.0), (3, 6, 2600.0, 2.0),
            (1, 4, 2500.0, 2.0), (4, 7, 1750.0, 1.0),
            (2, 5, 1800.0, 1.0), (5, 8, 2600.0, 2.0),
            (1, 3, 1500.0, 0.8), (4, 6, 1600.0, 0.9),
            (2, 4, 1650.0, 0.9), (5, 7, 1550.0, 0.8),
        ]
        ods = [(0, 8, 4200.0), (1, 8, 1500.0), (0, 7, 1300.0)]
    else:
        raise ValueError(name)
    graph = nx.DiGraph(); edge_idx = {}
    for e, (u, v, cap, t0) in enumerate(edges):
        graph.add_edge(u, v, weight=t0, eidx=e); edge_idx[(u, v)] = e
    od_list = []
    for origin, dest, demand in ods:
        node_paths = list(nx.all_simple_paths(graph, origin, dest, cutoff=8))
        edge_paths = [[edge_idx[(p[i], p[i + 1])] for i in range(len(p) - 1)] for p in node_paths]
        od_list.append({"o": origin, "d": dest, "demand": demand, "paths": edge_paths})
    return {
        "edges": len(edges),
        "t0": np.array([e[3] for e in edges]),
        "cap": np.array([e[2] for e in edges]),
        "od_list": od_list,
    }, (graph, edge_idx)


def solve_benchmark(sf, oracle, mode):
    total = sum(o["demand"] for o in sf["od_list"])
    routes = SE.RouteSet(sf["od_list"])
    game = SE.SFGame(sf, 0.5, [1.0], LAM_TIMES_M / (0.9 * total), routes, mode=mode)
    sol = game.solve(tol=1e-7)
    violations, worst = SE.certificate(game, sol, *oracle, 1e-6)
    if violations or sol["gap"] > 1e-6:
        raise RuntimeError(f"{mode} benchmark failed")
    sol.update(certified=True, max_reduced_cost=float(worst), route_cost_scale=max(float(np.max(game.A.T @ game.edge_weights(sol['fH'], sol['fF'])[0])), 1e-12))
    return sol, game


def main():
    profiles = []
    benchmarks = {}
    for network in ("braess", "grid"):
        sf, oracle = build_network(network)
        ue, ue_game = solve_benchmark(sf, oracle, "ue")
        so, _ = solve_benchmark(sf, oracle, "so")
        benchmarks[network] = {"J_UE": ue["J"], "J_SO": so["J"], "nodes": len(oracle[0]), "edges": sf["edges"], "ods": len(sf["od_list"])}
        slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
        total = sum(o["demand"] for o in sf["od_list"])
        lam_beta = LAM_TIMES_M / (0.9 * total)
        routes = SE.RouteSet(sf["od_list"])
        for alpha in ALPHAS:
            for gamma in GAMMAS:
                for pid, partition in enumerate(PARTITIONS, 1):
                    sol, game = solve_profile(
                        sf, routes, oracle, alpha, SHARES, lam_beta,
                        coalitions=partition, objective_types=TYPES,
                        coordination_gamma=gamma, management_slope=slope,
                        max_rounds=8,
                    )
                    profiles.append({
                        "network": network, "alpha": alpha, "gamma": gamma,
                        "partition_id": f"P{pid:02d}",
                        "partition": [list(b) for b in partition],
                        "n_coalitions": len(partition),
                        "J": sol["J"],
                        "recovery": (ue["J"] - sol["J"]) / (ue["J"] - so["J"]),
                        "VI_gap": sol["gap"],
                        "omitted_route_slack": sol["max_reduced_cost"],
                        "route_cost_scale": sol["route_cost_scale"],
                        "certified": bool(sol["certified"]),
                        **company_accounting_costs(sol, game),
                    })
                    print(network, alpha, gamma, pid, len(partition), sol["gap"], sol["certified"], flush=True)
    raw = {"protocol": {"alphas": ALPHAS, "gammas": GAMMAS, "objective_types": TYPES.tolist(), "partitions": len(PARTITIONS)}, "benchmarks": benchmarks, "profiles": profiles}
    OUT.write_text(json.dumps(raw, indent=2) + "\n")
    states = []
    for network in benchmarks:
        for alpha in ALPHAS:
            for gamma in GAMMAS:
                rows = [r for r in profiles if r["network"] == network and r["alpha"] == alpha and r["gamma"] == gamma]
                best = min(rows, key=lambda r: r["J"])
                singleton = next(r for r in rows if r["n_coalitions"] == 4)
                grand = next(r for r in rows if r["n_coalitions"] == 1)
                states.append({
                    "network": network, "alpha": alpha, "gamma": gamma,
                    "best_partition_id": best["partition_id"], "best_K": best["n_coalitions"],
                    "best_recovery": best["recovery"],
                    "gain_vs_singleton_pp": 100 * (best["recovery"] - singleton["recovery"]),
                    "grand_regret_pp": 100 * (best["recovery"] - grand["recovery"]),
                })
    summary = {"states": states, "certified_rows": sum(r["certified"] for r in profiles), "rows": len(profiles), "max_VI_gap": max(r["VI_gap"] for r in profiles), "max_relative_omitted_route_slack": max(r["omitted_route_slack"] / r["route_cost_scale"] for r in profiles)}
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
