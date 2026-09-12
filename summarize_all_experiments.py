"""Create a single auditable summary of completed baseline and wildfire results."""
from __future__ import annotations

import json
import statistics
from pathlib import Path


OUT = Path("all_experiments_summary.json")


def full_partition_summary(path, label):
    data = json.loads(Path(path).read_text())
    rows = data["profiles"]
    states = []
    for key in sorted(set((r["scenario_id"], r["alpha"], r["gamma"]) for r in rows)):
        group = [r for r in rows if (r["scenario_id"], r["alpha"], r["gamma"]) == key]
        best = min(group, key=lambda r: r["J"])
        crit = min(group, key=lambda r: r["critical_vc_max"])
        states.append({
            "source": label,
            "scenario_id": key[0],
            "alpha": key[1],
            "gamma": key[2],
            "best_partition_id": best["partition_id"],
            "best_K": best["n_coalitions"],
            "best_J": best["J"],
            "critical_best_partition_id": crit["partition_id"],
            "critical_best_vc": crit["critical_vc_max"],
            "critical_best_J": crit["J"],
            "J_range": max(r["J"] for r in group) - min(r["J"] for r in group),
            "critical_vc_range": max(r["critical_vc_max"] for r in group) - min(r["critical_vc_max"] for r in group),
            "certified_rows": sum(r["certified"] for r in group),
        })
    return {"source": label, "rows": len(rows), "certified_rows": sum(r["certified"] for r in rows), "states": states}


def main():
    mild = full_partition_summary("wildfire_full_partition_audit_results.json", "mild_and_critical")
    multiple = full_partition_summary("wildfire_full_partition_multiple_results.json", "multiple_corridor")
    info = json.loads(Path("wildfire_information_delay_results.json").read_text())["profiles"]
    information = []
    for case in sorted({r["information_case"] for r in info}):
        group = [r for r in info if r["information_case"] == case]
        information.append({
            "information_case": case,
            "rows": len(group),
            "certified_rows": sum(r["certified"] for r in group),
            "mean_physical_J": statistics.mean(r["J_physical"] for r in group),
            "mean_true_risk_cost": statistics.mean(r["true_risk_cost"] for r in group),
            "mean_critical_vc": statistics.mean(r["critical_vc_max"] for r in group),
        })
    prototype = {}
    for fn in ("prototype/exp_pred_n6.json", "prototype/exp_pred_n8.json"):
        rows = json.loads(Path(fn).read_text())
        prototype[fn] = {
            "rows": len(rows),
            "gamma_05_sign_accuracy_mean": statistics.mean(r["sign_acc"] for r in rows if r.get("gamma") == 0.5),
            "gamma_05_pred2_k2_exact": sum(r["pred2_k2"][0] == 0 for r in rows if r.get("gamma") == 0.5),
            "gamma_05_pred2_k2_mean_solves": statistics.mean(r["pred2_k2"][1] for r in rows if r.get("gamma") == 0.5),
            "gamma_1_sign_accuracy_mean": (statistics.mean(r["sign_acc"] for r in rows if r.get("gamma") == 1.0)
                                           if any(r.get("gamma") == 1.0 for r in rows) else None),
        }
    result = {
        "status": "completed_experiments_summary",
        "wildfire_full_partition": mild,
        "wildfire_multiple_corridor": multiple,
        "information_delay": information,
        "merge_score_prototype": prototype,
        "interpretation": {
            "full_partition": "Hazard topology and AV share change the TSTT-best coalition structure; critical-v/c constrained choices can differ.",
            "information": "Stale public risk information lowers perceived cost but increases true risk cost and critical v/c.",
            "merge_score": "Fixed-support prediction is strong at gamma=0.5 and unreliable near gamma=1; use fallback certification near degeneracy.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
