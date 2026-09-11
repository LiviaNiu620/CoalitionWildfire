# CoalitionDesign — Coalition Design in Mixed-Autonomy Networks

Research repository for *Coalition Design under Heterogeneous AV Objectives: Regimes and Low-Regret Selection in Mixed-Autonomy Networks* (Luyao Niu, Ruolin Li).

## Layout

| Path | Contents |
|---|---|
| `paper_response_interoperability.tex`, `refs.bib`, `response_paper/` | Manuscript source (main file, chapters, figures and figure generators) |
| `docs/Coalition_1_compiled_2026-09-11.pdf` | Compiled manuscript snapshot (2026-09-11; note: bibliography not resolved in this build) |
| `sioux_falls_data/`, `sioux_falls_loader.py` | Sioux Falls network inputs |
| `sf_*.py` / `sf_*_results.json` / `sf_*_summary.json` | Sioux Falls solver, scans, and certified outputs |
| `cross_network_*` | Braess / directed-grid regime audits and policy analysis |
| `proof_*.py`, `verify_*.py`, `regime_*.py` | Symbolic checks and proof-support scripts |
| `validate_*.py` | Fail-closed validators for result packages |
| `plan/` | Research plan and experiment protocol (wildfire-disruption extension) |
| `docs/claude_review/` | Review of the manuscript and detailed revision plan (2026-09-11) |
| `docs/notes/` | Short summaries of review findings and revision plan |
| `prototype/` | Small-network prototype supporting the revision plan: frozen vs self-consistent slopes, adjoint merge scores, predicted two-step policy, n=6/8 scalability |
| `RESTORE_MANIFEST.md` | Provenance of the restored paper package |

## Status (2026-09-11)

See `docs/claude_review/02_revision_plan_with_prototype.md` for the current revision priorities
(frozen-slope robustness, γ=1 curvature issue, adjoint-score policies, scalability, larger networks).
