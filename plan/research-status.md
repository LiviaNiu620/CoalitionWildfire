# Research Status

Date: 2026-09-11

## Completed

- Restored and filtered the paper/code/data package from commit `9ff72248071413fea4b9a6c895d1bd140f88121e`.
- Ran all existing numerical package validators successfully.
- Ran the symbolic canonical-regime check successfully.
- Added and ran `validate_potential_structure.py` on a two-edge, one-OD, two-company coalition instance.
- Potential certificate: maximum finite-difference Jacobian asymmetry `1.11e-10`; maximum potential-gradient error `3.58e-9`; minimum symmetric-part eigenvalue `7.50e-7`; passed.
- Created the wildfire research plan, experiment protocol, and method-experiment traceability matrix.
- Ran the Phase 1 implementation smoke test with ordinary HDVs, fixed commercial AV OD obligations, singleton/grand coalitions, capacity degradation, and a single-edge closure. All six fixed-route profiles converged with VI gaps below `1.2e-9`.
- Ran the reduced-OD Sioux Falls wildfire pilot: 240/240 profiles certified across four hazard snapshots, two AV penetrations, two governance intensities, and all 15 n=4 partitions. The pilot showed hazard-dependent best-K changes, but remains planning evidence because it uses 20 OD pairs.
- Rebuilt V0/V2 on the same spatial scenario catalog: 240/240 profiles per calibration, with no best-K or best-partition changes across the 16 pilot states and essentially identical grand regret. This is a valid null robustness result; the earlier apparent ranking shifts mixed different scenario catalogs and are not used.
- Replaced the incidence-ranked hazard catalog with a spatial-footprint catalog: footprint nodes, interior closure candidates, boundary degraded edges, and separate critical edges are recorded. The V0/V2 comparison is now controlled for scenario geometry.
- Added `public_edge_penalty` to the solver as a fixed common perceived-cost term for future risk/information experiments and extended the potential certificate to a nonzero penalty.
- Completed the full n=4 hazard audit: 120/120 full-OD profiles across mild capacity degradation and critical-edge closure, all 15 partitions, two AV penetrations, and two governance intensities are certified. The formal epsilon selector now runs on all 15 partitions per state.
- Completed the perceived-state/true-state information-delay pilot: 180/180 certified reduced-OD profiles across perfect, half-strength stale, and no-risk-information cases. Stale information lowers perceived penalties but increases true risk-weighted flow and critical v/c.
- Completed the third multiple-corridor full partition audit: 60/60 full-OD profiles certified; the combined wildfire hazard audit now contains 180/180 full-OD profiles. Added `all_experiments_summary.json` to aggregate certified baseline, wildfire, information-delay, and merge-score prototype evidence.
- Merge-score prototype evidence is archived but not rerun: Anaconda SciPy binary incompatibility prevents a clean rerun, so prototype numbers remain historical evidence and are not treated as new results.
- Added `all_experiments_summary.json` and its validator to aggregate baseline, full-OD wildfire, information-delay, slope-robustness, and historical prototype evidence. The current combined evidence supports state-dependent coalition ranking, not universal grand/partial optimality.
- Completed the full 15-partition mild/critical audit and the third multiple-corridor full audit: 180/180 full-OD profiles across three hazard classes, two AV penetrations, two governance intensities, and all n=4 partitions are certified. Full-partition epsilon selection is available through the selector's `--raw` input.
- Ran and passed a reduced-OD no-hazard adapter regression harness. A full 528-OD 60-profile regression was attempted and interrupted after roughly six minutes because the current route-based projection solver produced no checkpoint output; it generated no result file and is recorded as a scalability bottleneck rather than a numerical pass.
- Completed the full-OD endpoint hazard audit: 32/32 certified profiles for mild capacity degradation and critical-edge closure across representative grand/partial/singleton partitions, two AV penetrations, and two governance intensities.
- Added full-OD endpoint Pareto and epsilon-constraint diagnostics. Critical-edge limits can select a different endpoint from TSTT minimization, but these are four-partition pilot diagnostics and not final policy thresholds.
- Corrected the manuscript's potential/formation wording and added an explicit information-structure assumption, BPR fixed-reference limitation, and gamma=1 tangent-curvature caveat.
- Added the optional solver switch `management_slope_mode="state"` for state-dependent BPR-slope diagnostics while preserving the default `"reference"` baseline behavior. The state-mode toy solve converged with zero reported VI gap and accounting-cost extraction succeeded.
- Completed the full n=4 hazard audit and information-delay pilot: 180/180 full-OD hazard profiles and 180/180 reduced-OD information profiles are certified.
- Re-ran the prototype symbolic check successfully: frozen-slope Jacobian asymmetry is exactly zero, while true state-dependent atomic internalization is asymmetric. Archived merge-score numeric results are summarized in `merge_score_prototype_summary.json`; a clean numeric rerun is blocked by the local SciPy binary mismatch and is not claimed as new evidence.

## Current gate

The original model is baseline-complete and the static wildfire snapshot layer is experimentally complete for three synthetic hazard classes. Final wildfire claims remain conditional on synthetic scenario construction; LA calibration, AV-specific route access, score-guided policy, and larger-n scalability remain pending.

## Do not claim yet

- real-LA wildfire coalition ranking;
- household evacuation or life-safety improvement;
- universal hazard-dependent grand/partial transition;
- computational savings from merge-score selection in the main solver.

## Next executable task

Implement and validate a main-solver merge-score policy, then add AV-specific feasible-edge sets and n=6 scalability before calibrating synthetic hazard states with public Los Angeles data.
