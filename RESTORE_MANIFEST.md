# Restored paper package

This directory is a filtered restore from commit `9ff72248071413fea4b9a6c895d1bd140f88121e` of `lauraNiu/Mixed_Autonomy_Coalition` (the local checkout's configured remote for the repository).

Included material:

- `paper_response_interoperability.tex`, `refs.bib`, and the existing `response_paper/` chapter and figure tree;
- Sioux Falls network inputs under `sioux_falls_data/`;
- the `sf_*` experiment, solver, scan, recomputation, table, and validation scripts and their committed JSON outputs;
- Braess/directed-grid cross-network scripts and committed summaries;
- coalition-objective, response/regime proof, symbolic-check, and validation support code needed by the manuscript's reproducibility appendix;
- the response-paper figure generators and committed PDF/PNG/SVG/JSON figure assets;
- the paper response letter and small table/data manifests.

Intentionally excluded are unrelated or superseded paper families, meeting slides, proposal/technote artifacts, historical planning notes, archived figures, compiled PDFs outside `response_paper/figures/`, and temporary page-render images.

The restore is faithful to the selected paths as they existed at the commit. In particular, that commit deleted `response_paper/chapters/04_ownership_interoperability.tex` while `paper_response_interoperability.tex` still contains an `\\input{response_paper/chapters/04_ownership_interoperability}` line. The missing chapter was not reconstructed from an earlier commit, so the source tree preserves the commit's exact state rather than silently mixing revisions.
