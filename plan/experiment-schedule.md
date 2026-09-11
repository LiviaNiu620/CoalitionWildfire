# Experiment Schedule and Gates

## Stage A: Integrity baseline (completed)

Purpose: freeze the original coalition mechanism before application changes.

Tasks: run existing validators; verify potential gradient and Jacobian symmetry; preserve baseline JSON; document the fixed-reference BPR interpretation.

Outputs: `potential_structure_check.json`, existing certified result packages, and validator logs.

Gate A: passed.

## Stage B: Sioux Falls hazard adapter

Purpose: provide the minimum executable wildfire application.

Tasks:

1. Define `hazard_scenarios.json` schema.
2. Implement edge closure, capacity/free-flow multipliers, risk attributes, public-information metadata, and critical-edge metadata.
3. Ensure the shortest-path oracle excludes closed edges.
4. Add no-hazard regression against the original solver.
5. Add unit tests for demand conservation and route validity.

Data: current Sioux Falls TNTP files plus synthetic scenario definitions.

Gate B: no-hazard results match baseline within `1e-8` relative TSTT; every hazard profile passes VI and omitted-route tolerance `1e-6`; no route contains a closed edge.

## Stage C: Pilot partition audit

Purpose: test whether hazard state can change coalition ranking before launching a large grid.

Design: n=4, all 15 partitions, ordinary HDVs, fixed commercial AV OD, two AV loads, two governance intensities, and four hazard scenarios (none, capacity loss, one critical closure, multiple degradation). Run fixed-reference and state-slope diagnostic modes separately.

Data: Sioux Falls plus Stage B scenarios.

Gate C: all profiles certified; at least one nontrivial hazard comparison has measurable effect above solver tolerance; report null results if no ranking changes.

## Stage D: Main hazard-state regime scan

Purpose: estimate the state map for singleton, partial, and grand coordination.

Design: hazard severity x disruption topology x AV load x order-corridor exposure x objective composition x gamma. Use complete n=4 audits for selected states and Latin-hypercube or stratified sampling for the broader state space.

Data: calibrated synthetic order portfolios, scenario definitions, and public road-state ranges.

Gate D: raw/summary/validator chain complete; best partitions reported as tolerance-based sets; effect sizes reported in absolute time, percent TSTT, and recovery points where defined.

## Stage E: Robustness

Purpose: identify calibration artifacts and information sensitivity.

Tasks: reference-curvature V0/V1/V2 comparisons; public information delay; noisy regulator estimates; optional bounded-rational HDV sensitivity; order-retention fraction sensitivity.

Data: Stage D profiles plus generated error draws with recorded seeds.

Gate E: report cross-model regret, ranking agreement, transition movement, and slope-mismatch diagnostics. High-load partial claims survive only if supported across the declared robustness range.

## Stage F: Score-guided selection

Purpose: connect the response theorem to computational policy.

Tasks: implement adjoint/fixed-support merge scores; compare signs and rankings with certified finite merges; use fallback solves when gamma is near one, curvature is singular, the score is near zero, or support changes.

Data: Stage D/E current-state equilibria and candidate merge outcomes.

Gate F: report sign accuracy, Spearman ranking, regret, solve counts, runtime, and fallback rate.

## Stage G: Scale and empirical grounding

Purpose: test computational value and limited real-world validity.

Tasks: n=6 exhaustive, n=8 selected-state, n>=10 heuristic; construct an LA network case; calibrate capacity/closure scenarios with public data; seek operator order data or use labeled calibrated synthetic portfolios.

Gate G: distinguish observed, calibrated, and synthetic fields in every dataset; report runtime and coverage; make no empirical company-behavior claim without operator data.
