"""Design-regret, near-tie, greedy-merge, and governance-cost analysis."""
from __future__ import annotations

import itertools
import json
from math import comb
from pathlib import Path


RAW = Path("sf_coalition_partition_scan_results.json")
OUT = Path("sf_coalition_design_analysis.json")
EPS_REL = 1e-8


def canon(partition):
    return tuple(sorted(tuple(sorted(block)) for block in partition))


def coarsenings(partition):
    blocks = [tuple(block) for block in partition]
    for i, j in itertools.combinations(range(len(blocks)), 2):
        merged = tuple(sorted(blocks[i] + blocks[j]))
        yield canon([blocks[k] for k in range(len(blocks)) if k not in (i, j)] + [merged])


def coordination_burden(row):
    n = row["n_companies"]
    pairs = sum(comb(len(block), 2) for block in row["partition"])
    return row["gamma"] ** 2 * pairs / comb(n, 2)


def lower_envelope(rows, gap):
    # Loss = residual congestion gap + lambda * normalized coordination burden.
    intercept = {r["partition_id"]: (r["J"] - min(x["J"] for x in rows)) / gap for r in rows}
    slope = {r["partition_id"]: coordination_burden(r) for r in rows}
    cuts = {0.0}
    for a, b in itertools.combinations(rows, 2):
        ds = slope[a["partition_id"]] - slope[b["partition_id"]]
        if abs(ds) < 1e-14:
            continue
        lam = (intercept[b["partition_id"]] - intercept[a["partition_id"]]) / ds
        if lam > 0:
            cuts.add(float(lam))
    cuts = sorted(cuts)
    probes = []
    for left, right in zip(cuts, cuts[1:] + [None]):
        probe = left + 1.0 if right is None else (left + right) / 2
        values = {r["partition_id"]: intercept[r["partition_id"]] + probe * slope[r["partition_id"]] for r in rows}
        best = min(values.values())
        ids = sorted(pid for pid, value in values.items() if abs(value - best) <= 1e-10)
        if probes and probes[-1]["optimal_partition_ids"] == ids:
            probes[-1]["lambda_max"] = right
        else:
            probes.append({"lambda_min": left, "lambda_max": right, "optimal_partition_ids": ids})
    return probes


def main():
    raw = json.loads(RAW.read_text())
    j_ue, j_so = raw["J_UE"], raw["J_SO"]
    gap = j_ue - j_so
    states = []
    policy_regret = {name: [] for name in ("greedy_merge", "two_step_merge", "always_grand", "always_singleton")}
    for alpha in sorted({r["alpha"] for r in raw["profiles"]}):
        for gamma in sorted({r["gamma"] for r in raw["profiles"]}):
            rows = [r for r in raw["profiles"] if r["alpha"] == alpha and r["gamma"] == gamma]
            by_partition = {canon(r["partition"]): r for r in rows}
            best_j = min(r["J"] for r in rows)
            eps = EPS_REL * j_ue
            optimal = sorted(r["partition_id"] for r in rows if r["J"] - best_j <= eps)
            ranked = sorted(rows, key=lambda r: r["J"])

            current = canon([(i,) for i in range(4)])
            greedy_path = [by_partition[current]["partition_id"]]
            while len(current) > 1:
                candidates = [by_partition[p] for p in coarsenings(current)]
                choice = min(candidates, key=lambda r: r["J"])
                if choice["J"] >= by_partition[current]["J"] - eps:
                    break
                current = canon(choice["partition"])
                greedy_path.append(choice["partition_id"])
            greedy = by_partition[current]
            current2 = canon([(i,) for i in range(4)])
            lookahead_path = [by_partition[current2]["partition_id"]]
            while len(current2) > 1:
                candidates = []
                for first in coarsenings(current2):
                    candidates.append((by_partition[first]["J"], first))
                    for second in coarsenings(first):
                        candidates.append((by_partition[second]["J"], first))
                future_j, next_partition = min(candidates, key=lambda item: item[0])
                if future_j >= by_partition[current2]["J"] - eps:
                    break
                current2 = next_partition
                lookahead_path.append(by_partition[current2]["partition_id"])
            lookahead = by_partition[current2]
            grand = next(r for r in rows if r["n_coalitions"] == 1)
            singleton = next(r for r in rows if r["n_coalitions"] == 4)

            def regret(row):
                return {
                    "delta_J": row["J"] - best_j,
                    "delta_J_pct_UE": 100 * (row["J"] - best_j) / j_ue,
                    "recovery_regret_pp": 100 * (ranked[0]["recovery"] - row["recovery"]),
                }

            policy_regret["greedy_merge"].append(regret(greedy))
            policy_regret["two_step_merge"].append(regret(lookahead))
            policy_regret["always_grand"].append(regret(grand))
            policy_regret["always_singleton"].append(regret(singleton))
            states.append({
                "alpha": alpha,
                "gamma": gamma,
                "global_optimal_partition_ids": optimal,
                "global_best_J": best_j,
                "best_minus_second_J": ranked[1]["J"] - ranked[0]["J"],
                "best_minus_second_pct_UE": 100 * (ranked[1]["J"] - ranked[0]["J"]) / j_ue,
                "greedy_path": greedy_path,
                "greedy_partition_id": greedy["partition_id"],
                "greedy_regret": regret(greedy),
                "two_step_path": lookahead_path,
                "two_step_partition_id": lookahead["partition_id"],
                "two_step_regret": regret(lookahead),
                "always_grand_regret": regret(grand),
                "always_singleton_regret": regret(singleton),
                "coordination_cost_regimes": lower_envelope(rows, gap),
            })
    summary = {}
    for name, values in policy_regret.items():
        summary[name] = {
            "mean_recovery_regret_pp": sum(v["recovery_regret_pp"] for v in values) / len(values),
            "max_recovery_regret_pp": max(v["recovery_regret_pp"] for v in values),
            "mean_delta_J_pct_UE": sum(v["delta_J_pct_UE"] for v in values) / len(values),
            "max_delta_J_pct_UE": max(v["delta_J_pct_UE"] for v in values),
        }
    output = {
        "protocol": {
            "source": str(RAW),
            "near_tie_relative_TSTT_tolerance": EPS_REL,
            "greedy_rule": "start from singletons; choose the best improving pairwise merge; stop when no merge improves TSTT",
            "coordination_burden": "gamma^2 times within-coalition company-pair share",
            "penalized_loss": "TSTT gap from state best divided by UE-SO gap plus lambda times coordination burden",
        },
        "states": states,
        "policy_regret_summary": summary,
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
