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
- Corrected the manuscript's potential/formation wording and added an explicit information-structure assumption, BPR fixed-reference limitation, and gamma=1 tangent-curvature caveat.
- Added the optional solver switch `management_slope_mode="state"` for state-dependent BPR-slope diagnostics while preserving the default `"reference"` baseline behavior. The state-mode toy solve converged with zero reported VI gap and accounting-cost extraction succeeded.

## Current gate

The original model is baseline-complete but wildfire-incomplete. No wildfire claim is supported yet. A minimal capacity/closure smoke adapter now exists, but the Sioux Falls hazard-state adapter is still pending. It must alter edge availability, capacity, free-flow time, risk attributes, and critical-edge metadata while preserving the existing no-hazard output exactly.

## Do not claim yet

- wildfire-specific coalition ranking;
- hazard-dependent grand/partial transition;
- emergency-corridor benefit;
- robustness to hazard-state or regulator information error;
- computational savings from merge-score selection.

## Next executable task

Implement the Sioux Falls static hazard-snapshot adapter and run no-hazard regression, single-closure, capacity-degradation, and critical-edge scenarios before scaling to full partition audits.
