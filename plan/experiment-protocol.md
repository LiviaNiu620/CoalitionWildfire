# Experiment Protocol

## Research objective

Evaluate whether a regulator-mandated temporary operational coalition among commercial AV operators changes congestion and emergency-critical corridor loading in wildfire-disrupted networks when accepted AV order OD obligations are preserved and ordinary HDVs remain informed but uncontrolled.

## Phase 0 baseline experiments

Dataset: committed Sioux Falls network, Braess network, directed grid, and committed machine-readable result files.

Baselines: singleton partition, all declared n=4 partitions, grand partition, current greedy merge, current two-step merge, and the existing reference-curvature regimes.

Metrics: TSTT, recovery, VI gap, omitted-route slack, partition regret, HHI, and certified profile count.

Acceptance: existing validators pass; no committed result changes; potential residual and conservation checks pass.

## Phase 1 static hazard snapshots

Dataset: the Sioux Falls network plus spatial-footprint scenario files generated from explicit edge-state rules. Each scenario records `scenario_id`, footprint nodes, open/closed edges, degraded edges, capacity multiplier, free-flow multiplier, risk score, public-information delay, and a critical-edge set that is disjoint from degraded edges whenever the network permits.

Factors:

- hazard severity: none, mild, moderate, severe;
- disruption pattern: local closure, corridor closure, multiple closures, asymmetric north/south degradation;
- AV order load: low, medium, high;
- order exposure: low, mixed, high overlap with critical edges;
- AV penetration: 0.1, 0.3, 0.5, 0.7, 0.9;
- governance intensity: 0, 0.25, 0.5, 0.75, 1;
- partition: all 15 n=4 partitions for audit states.

Metrics: network TSTT, commercial AV order delay/detour, critical-edge overload, risk-weighted road use, governance burden, best partition set, best K, and regret of fixed grand/singleton policies.

Acceptance: no-hazard outputs match Phase 0; every hazard profile has finite costs, valid routes, VI gap <= 1e-6, and relative omitted-route slack <= 1e-6.

## Phase 2 reference-curvature robustness

Treat the same-hazard-state HDV-only derivative as the primary wildfire reference. Compare it with the legacy base-network derivative and with `bbar = kappa_b * c_prime(x_ref)` for `kappa_b in {0.25, 0.5, 1, 2, 4}` and `x_ref` equal to hazard-state HDV-only and no-coordination reference flows. Re-run complete n=4 audits for selected hazard/order states.

Report: best K, best partition, grand regret, partial-optimum frequency, transition boundary, and sensitivity intervals. Do not update `bbar` inside an iteration unless a new fixed-point model is explicitly specified.

## Phase 3 information and HDV robustness

Use the optional fixed `public_edge_penalty` interface for common public risk/closure penalties. Perturb public road information with zero, one-window, and multi-window delay by separating perceived state from realized evaluation state. Compare deterministic Wardrop with a bounded-rational or logit route-choice sensitivity. Keep the potential model as the primary benchmark and label alternatives as behavioral robustness.

Report: partition changes, TSTT change, critical-edge overload, AV order delay, and regret relative to the full-information policy.

## Phase 4 regulator information error

Generate true states from `tau`, order OD, and hazard scenarios. Let the regulator select a partition using noisy `tau_hat`, `d_hat`, or hazard state `s_hat`; evaluate the selected partition under the true state.

Report: selection regret, false-grand rate, false-partial rate, critical-edge overload, and AV service loss versus error magnitude.

## Phase 5 score-guided policy

For each current partition, estimate the fixed-support merge-direction score for feasible pairwise merges. Rank candidates by predicted improvement, then certify the selected finite merge with the continuation solver.

Report: score sign accuracy, rank correlation with realized finite TSTT change, regret, saved equilibrium solves, support-change rate, and false-positive/false-negative counts.

## Phase 6 scalability

Use n=6 exhaustive partitions where feasible, n=8 selected-state audits, and n>=10 score-guided or beam-search policies. Report runtime, memory, equilibrium solve count, partition coverage, and regret against the best available audit.

## Reproducibility contract

Every experiment must write a raw JSON or CSV profile file, a summary file, a validator, a figure/table generator, and a protocol record containing scenario definitions, seeds, tolerances, software versions, and source data identifiers. Planning or synthetic values must never enter manuscript results prose before being replaced by labeled real outputs.
