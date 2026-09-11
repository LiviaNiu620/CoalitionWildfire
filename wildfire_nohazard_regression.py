"""Regression of the hazard adapter against the committed full Sioux Falls scan."""
from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np

from sf_coalition_partition_scan import PARTITIONS, PARTITION_IDS
from sf_revision_common import LAM_TIMES_M, SE, network_oracle, solve_benchmarks, solve_profile
from wildfire_scenarios import apply_snapshot, build_snapshot, hazard_oracle


BASELINE = Path("sf_coalition_partition_scan_results.json")
OUT = Path("wildfire_nohazard_regression_results.json")
SHARES = np.full(4, 0.25)
TYPES = np.array([0.5, 0.5, 1.5, 1.5])
ALPHAS = (0.5, 0.9)
GAMMAS = (0.5, 1.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--od-limit", type=int, default=20)
    parser.add_argument("--max-profiles", type=int, default=4)
    args = parser.parse_args()
    od_limit = None if args.od_limit <= 0 else args.od_limit
    base, scenarios = build_snapshot(od_limit=od_limit, demand_scale=1.0)
    no_hazard = next(s for s in scenarios if s["scenario_id"] == "no_hazard")
    sf = apply_snapshot(base, no_hazard)
    total = sum(item["demand"] for item in sf["od_list"])
    lam_beta = LAM_TIMES_M / (0.9 * total)
    ue, ue_game, so = solve_benchmarks(sf)
    slope = SE.bpr_deriv(ue["x"], ue_game.t0, ue_game.cap)
    routes = SE.RouteSet(sf["od_list"])
    oracle = hazard_oracle(no_hazard)
    baseline = json.loads(BASELINE.read_text()) if od_limit is None else None
    lookup = ({(r["alpha"], r["gamma"], r["partition_id"]): r for r in baseline["profiles"]}
              if baseline is not None else {})
    rows = []
    target = min(args.max_profiles, len(ALPHAS) * len(GAMMAS) * len(PARTITIONS))
    for alpha in ALPHAS:
        for gamma in GAMMAS:
            warm = {}
            for partition in PARTITIONS:
                pid = PARTITION_IDS[partition]
                sol, game = solve_profile(
                    sf, routes, oracle, alpha, SHARES, lam_beta,
                    coalitions=partition, objective_types=TYPES,
                    coordination_gamma=gamma, management_slope=slope,
                    warm=warm.get(pid), max_rounds=25,
                )
                row = {
                    "alpha": alpha,
                    "gamma": gamma,
                    "partition_id": pid,
                    "new_J": float(sol["J"]),
                    "VI_gap": float(sol["gap"]),
                    "certified": bool(sol["certified"]),
                }
                if baseline is not None:
                    old = lookup[(alpha, gamma, pid)]
                    row.update(
                        baseline_J=float(old["J"]),
                        abs_delta_J=float(sol["J"] - old["J"]),
                        relative_delta_J=float(
                            (sol["J"] - old["J"]) / max(abs(old["J"]), 1e-12)
                        ),
                    )
                rows.append(row)
                warm[pid] = (sol["fH"], sol["fF"])
                # Checkpoint after every profile so long runs remain inspectable.
                OUT.write_text(json.dumps({"status": "running", "profiles": rows}, indent=2) + "\n")
                if len(rows) >= target:
                    break
            if len(rows) >= target:
                break
        if len(rows) >= target:
            break
    result = {
        "status": "regression",
        "source_baseline": str(BASELINE),
        "network": "Sioux Falls full OD" if od_limit is None else f"Sioux Falls OD limit {od_limit}",
        "scenario_id": "no_hazard",
        "profiles": rows,
        "baseline_comparable": baseline is not None,
        "max_abs_delta_J": (max(abs(r["abs_delta_J"]) for r in rows)
                            if baseline is not None else None),
        "max_relative_delta_J": (max(abs(r["relative_delta_J"]) for r in rows)
                                 if baseline is not None else None),
        "max_VI_gap": max(r["VI_gap"] for r in rows),
        "certified_rows": sum(r["certified"] for r in rows),
        "passed": bool(
            all(r["certified"] for r in rows)
            and (baseline is None or max(abs(r["relative_delta_J"]) for r in rows) <= 1e-8)
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
