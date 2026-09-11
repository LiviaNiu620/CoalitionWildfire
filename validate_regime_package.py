"""Fail-closed checks for regime, regret, and cross-network results."""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def close(actual, expected, tol=1e-9):
    if not math.isclose(actual, expected, rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"expected {expected}, found {actual}")

def finite(v):
    if isinstance(v, dict):
        for x in v.values(): finite(x)
    elif isinstance(v, list):
        for x in v: finite(x)
    elif isinstance(v, float) and not math.isfinite(v): raise AssertionError("nonfinite value")

def main():
    design = json.loads((ROOT / "sf_coalition_design_analysis.json").read_text())
    cross = json.loads((ROOT / "cross_network_regime_results.json").read_text())
    summary = json.loads((ROOT / "cross_network_regime_summary.json").read_text())
    if len(design["states"]) != 4: raise AssertionError("expected four Sioux Falls design states")
    if design["policy_regret_summary"]["greedy_merge"]["max_recovery_regret_pp"] > 1e-8:
        raise AssertionError("greedy merge failed to recover a global optimum")
    sf_tie = next(s for s in design["states"] if s["alpha"] == 0.9 and s["gamma"] == 0.5)
    if sf_tie["global_optimal_partition_ids"] != ["P04", "P05"]:
        raise AssertionError("Sioux Falls near-tie set changed")
    close(sf_tie["best_minus_second_J"], 5.578622221946716e-7)
    close(design["policy_regret_summary"]["always_grand"]["mean_recovery_regret_pp"], 4.907916380977595)
    close(design["policy_regret_summary"]["always_grand"]["max_recovery_regret_pp"], 18.694975105831336)
    close(design["policy_regret_summary"]["always_singleton"]["mean_recovery_regret_pp"], 13.547596572936506)
    sf_cost = next(s for s in design["states"] if s["alpha"] == 0.5 and s["gamma"] == 1.0)
    thresholds = [r["lambda_min"] for r in sf_cost["coordination_cost_regimes"]]
    for actual, expected in zip(thresholds, (0.0, 0.02934229348648187, 0.26776179961604457, 0.7434622137106147)):
        close(actual, expected)
    if summary["rows"] != 2 * 4 * 4 * 15: raise AssertionError("cross-network grid incomplete")
    if summary["certified_rows"] != summary["rows"]: raise AssertionError("uncertified cross-network profile")
    if summary["max_VI_gap"] > 1.0001e-6: raise AssertionError("VI tolerance exceeded")
    if summary["max_relative_omitted_route_slack"] > 1.0001e-6: raise AssertionError("oracle tolerance exceeded")
    if summary["policy_regret_summary"]["greedy_merge"]["states"] != 32:
        raise AssertionError("cross-network policy grid incomplete")
    if summary["policy_regret_summary"]["two_step_merge"]["max_recovery_regret_pp"] > 1e-8:
        raise AssertionError("two-step merge failed cross-network exactness check")
    cross_policy = summary["policy_regret_summary"]
    if cross_policy["greedy_merge"]["exact_state_count"] != 29:
        raise AssertionError("one-step exact-state count changed")
    if cross_policy["two_step_merge"]["exact_state_count"] != 32:
        raise AssertionError("two-step exact-state count changed")
    close(cross_policy["greedy_merge"]["mean_recovery_regret_pp"], 0.3315742189269947)
    close(cross_policy["greedy_merge"]["max_recovery_regret_pp"], 8.339364454550314)
    close(cross_policy["always_grand"]["max_recovery_regret_pp"], 18.205112281909464)
    close(cross_policy["always_singleton"]["max_recovery_regret_pp"], 42.28882706589836)
    states = summary["states"]
    for network, partial_count in (("braess", 2), ("grid", 3)):
        network_states = [s for s in states if s["network"] == network]
        if len(network_states) != 16: raise AssertionError(f"{network} state grid incomplete")
        if any(s["best_K"] != 1 for s in network_states if s["alpha"] <= 0.7):
            raise AssertionError(f"{network} low-penetration regime changed")
        if sum(s["best_K"] > 1 for s in network_states) != partial_count:
            raise AssertionError(f"{network} partial-regime count changed")
    finite(design); finite(cross); finite(summary)
    print("regime and cross-network package validation passed")

if __name__ == "__main__": main()
