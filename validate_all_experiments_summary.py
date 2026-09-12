"""Validator for the completed-experiment summary."""
from __future__ import annotations

import json
from pathlib import Path


d = json.loads((Path(__file__).resolve().parent / "all_experiments_summary.json").read_text())
assert d["status"] == "completed_experiments_summary"
assert d["wildfire_full_partition"]["certified_rows"] == 120
assert d["wildfire_multiple_corridor"]["certified_rows"] == 60
assert all(x["certified_rows"] == 60 for x in d["information_delay"])
assert d["merge_score_prototype"]
print("all experiments summary validation passed")
