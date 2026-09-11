"""Smoke test for wildfire-state capacity changes with ordinary HDVs.

This is a planning/implementation check, not manuscript evidence. It keeps
commercial AV order OD obligations fixed and changes only the road state.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sf_equilibrium import SFGame, RouteSet, bpr_cost


OUT = Path("wildfire_snapshot_smoke_results.json")
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
SHARES = np.full(4, 0.25)
PARTITIONS = {
    "singleton": ((0,), (1,), (2,), (3,)),
    "grand": ((0, 1, 2, 3),),
}


def solve_case(name: str, capacities: list[float], paths: list[list[int]]) -> dict:
    sf = {
        "edges": 3,
        "t0": np.array([1.0, 1.1, 1.2]),
        "cap": np.array(capacities, dtype=float),
        "od_list": [{"demand": 100.0, "paths": paths}],
    }
    rows = []
    for partition_name, partition in PARTITIONS.items():
        game = SFGame(
            sf,
            alpha=0.6,
            shares=SHARES,
            lam_beta=0.01,
            routes=RouteSet(sf["od_list"]),
            mode="eq",
            coalitions=partition,
            objective_types=TYPES,
            coordination_gamma=1.0,
            management_slope=[0.01, 0.01, 0.01],
        )
        sol = game.solve(seed=0, tol=1e-8, max_iter=12000, verbose=False)
        x = np.asarray(sol["x"], dtype=float)
        c = bpr_cost(x, game.t0, game.cap)
        rows.append(
            {
                "partition": partition_name,
                "n_coalitions": len(partition),
                "J": float(np.dot(x, c)),
                "x": x.tolist(),
                "cost": c.tolist(),
                "VI_gap": float(sol["gap"]),
                "iterations": int(sol["iters"]),
                "certified_on_fixed_routes": bool(sol["gap"] <= 1e-8),
            }
        )
    return {
        "scenario": name,
        "capacities": capacities,
        "open_paths": paths,
        "profiles": rows,
    }


def main() -> None:
    scenarios = [
        solve_case("no_hazard", [100.0, 100.0, 100.0], [[0], [1], [2]]),
        solve_case("capacity_degradation", [100.0, 45.0, 100.0], [[0], [1], [2]]),
        solve_case("single_edge_closure", [100.0, 100.0, 100.0], [[0], [1]]),
    ]
    result = {
        "status": "smoke_only",
        "planning_data": True,
        "model": "ordinary HDV + commercial AV orders + temporary coalition",
        "note": "Fixed-route smoke test; not wildfire evidence and not a full Sioux Falls certification.",
        "scenarios": scenarios,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
