"""Fail-closed validation for identity-preserving/pooling decomposition."""
from __future__ import annotations
import json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def finite(v):
    if isinstance(v,dict):
        for x in v.values(): finite(x)
    elif isinstance(v,list):
        for x in v: finite(x)
    elif isinstance(v,float) and not math.isfinite(v): raise AssertionError("nonfinite result")
def main():
    raw=json.loads((ROOT/"sf_coalition_pooling_results.json").read_text())
    summary=json.loads((ROOT/"sf_coalition_pooling_summary.json").read_text())
    rows=raw["profiles"]
    if len(rows)!=4: raise AssertionError("expected 4 pooling rows")
    if not all(r["certified"] for r in rows): raise AssertionError("uncertified pooling profile")
    if max(r["VI_gap"] for r in rows)>1.0001e-6: raise AssertionError("VI tolerance exceeded")
    if max(r["omitted_route_slack"]/r["route_cost_scale"] for r in rows)>1.0001e-6:
        raise AssertionError("omitted-route tolerance exceeded")
    finite(raw); finite(summary)
    print("coalition pooling package validation passed")
if __name__ == "__main__": main()
