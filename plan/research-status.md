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
- Ran and passed a reduced-OD no-hazard adapter regression harness. A full 528-OD 60-profile regression was attempted and interrupted after roughly six minutes because the current route-based projection solver produced no checkpoint output; it generated no result file and is recorded as a scalability bottleneck rather than a numerical pass.
- Completed the full-OD endpoint hazard audit: 32/32 certified profiles for mild capacity degradation and critical-edge closure across representative grand/partial/singleton partitions, two AV penetrations, and two governance intensities.
- Added full-OD endpoint Pareto and epsilon-constraint diagnostics. Critical-edge limits can select a different endpoint from TSTT minimization, but these are four-partition pilot diagnostics and not final policy thresholds.
- Corrected the manuscript's potential/formation wording and added an explicit information-structure assumption, BPR fixed-reference limitation, and gamma=1 tangent-curvature caveat.
- Added the optional solver switch `management_slope_mode="state"` for state-dependent BPR-slope diagnostics while preserving the default `"reference"` baseline behavior. The state-mode toy solve converged with zero reported VI gap and accounting-cost extraction succeeded.

## Current gate

The original model is baseline-complete but wildfire-incomplete. No final wildfire claim is supported yet. A static Sioux Falls hazard adapter and reduced-OD pilot now exist; the full-OD hazard adapter still needs incremental checkpointing, no-hazard equivalence, and runtime optimization before formal results are generated.

## Do not claim yet

- wildfire-specific coalition ranking;
- hazard-dependent grand/partial transition;
- emergency-corridor benefit;
- robustness to hazard-state or regulator information error;
- computational savings from merge-score selection.

## Next executable task

Optimize and checkpoint the full-OD Sioux Falls static hazard-snapshot adapter, then run no-hazard equivalence, single-closure, capacity-degradation, and critical-edge scenarios before scaling to full partition audits.
