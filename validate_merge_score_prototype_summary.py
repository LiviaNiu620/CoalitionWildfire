"""Validator for archived merge-score prototype evidence."""
from __future__ import annotations

import json
from pathlib import Path


d = json.loads((Path(__file__).resolve().parent / "merge_score_prototype_summary.json").read_text())
assert d["status"] == "historical_prototype_summary"
assert d["rerun_status"] == "not_rerun"
assert "gamma_1" in d["evidence"]
assert "n6_gamma_0_5" in d["evidence"]
print("merge-score prototype summary validation passed")
