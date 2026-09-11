"""Fail-closed validation for the complete n=4 coalition partition scan."""
from __future__ import annotations
import json, math
from pathlib import Path
ROOT = Path(__file__).resolve().parent
def finite(v):
    if isinstance(v, dict):
        for x in v.values(): finite(x)
    elif isinstance(v, list):
        for x in v: finite(x)
    elif isinstance(v, float) and not math.isfinite(v): raise AssertionError("nonfinite result")
def main():
    raw=json.loads((ROOT/"sf_coalition_partition_scan_results.json").read_text())
    summary=json.loads((ROOT/"sf_coalition_partition_scan_summary.json").read_text())
    rows=raw["profiles"]
    if len(rows)!=60: raise AssertionError(f"expected 60 rows, found {len(rows)}")
    if summary["partition_count"]!=15: raise AssertionError("expected 15 set partitions")
    if not all(r["certified"] for r in rows): raise AssertionError("uncertified partition profile")
    if max(r["VI_gap"] for r in rows)>1.0001e-6: raise AssertionError("VI tolerance exceeded")
    if max(r["omitted_route_slack"]/r["route_cost_scale"] for r in rows)>1.0001e-6:
        raise AssertionError("omitted-route tolerance exceeded")
    for alpha in (0.5,0.9):
        for gamma in (0.5,1.0):
            group=[r for r in rows if r["alpha"]==alpha and r["gamma"]==gamma]
            if len(group)!=15: raise AssertionError("partition grid incomplete")
            if {r["partition_id"] for r in group}!={f"P{i:02d}" for i in range(1,16)}:
                raise AssertionError("partition IDs incomplete")
    finite(raw); finite(summary)
    print("coalition partition scan package validation passed")
if __name__ == "__main__": main()
