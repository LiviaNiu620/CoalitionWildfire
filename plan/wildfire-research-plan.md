# Research Plan: Coalition Design for Commercial AV Traffic in Wildfire-Disrupted Networks

## Scope lock

The paper studies temporary regulator-mandated operational coalitions among commercial AV operators. AV companies continue serving accepted pickup-dropoff orders; member OD obligations and order ownership are preserved. The regulator selects a feasible coalition structure and coordination intensity but does not prescribe vehicle-level routes or reassign orders. HDVs remain ordinary human-driven traffic: publicly informed and route-responsive, but not directly controlled. Wildfire is an exogenous network and information state, not an evacuee-demand or rescue-fleet model.

Out of scope for this paper: household evacuation demand, rescue dispatch, shelter assignment, vulnerable-household pickup, legal firm merger, voluntary merge-and-split stability, fire-spread PDEs, and LLM emergency coordination.

## Model readiness decision

The existing coalition continuation model is mathematically usable as a baseline but is not yet wildfire-ready. The current code supports static networks, fixed OD demand, ordinary HDVs, identity-preserving partitions, symmetric PSD governance, and certified BPR continuation profiles. It does not yet support hazard-dependent edge sets/capacities, public hazard information, emergency-critical corridor metrics, hazard-specific reference-curvature calibration, or rolling updates.

The current model is therefore **baseline-complete, application-incomplete**. Experiments can begin now only as a baseline and solver-validation phase. Wildfire claims must wait until the hazard-state extension and its tests pass.

## Research phases

### Phase 0: Freeze and verify the baseline

Purpose: establish a reproducible reference before changing the model.

Verify the potential/VI equivalence, route-flow feasibility, column-generation certificate, partition scan, composition scan, degree scan, and existing cross-network policy outputs. Record exact package versions, tolerances, seeds, and commit hashes.

Exit criteria: all existing validators pass; no baseline JSON is overwritten; a potential residual test is added; current results are labeled `normal-network baseline`.

### Phase 1: Correct and freeze the static mathematical model

Purpose: make the model internally consistent before adding wildfire inputs.

Changes: state the convex potential for the fixed-reference governance family; distinguish general behavioral VI from the potential special case; state ordinary informed HDV behavior; define regulator-selected temporary operational coalitions; expose reference curvature as a parameter; retain identity-preserving OD constraints.

Exit criteria: symbolic Jacobian/potential residual is zero; numerical gradient checks pass; a written information structure and assumption ledger are approved; no claim says the BPR experiment is the exact nonlinear own-travel-time atomic game.

### Phase 2: Add static wildfire-state snapshots

Purpose: introduce wildfire relevance without building a full dynamic fire simulator.

Add scenario inputs for open/closed edges, capacity reduction, free-flow-time changes, risk penalties, public road information, and emergency-critical edge sets. Keep orders and ordinary HDV demand fixed within each snapshot.

Exit criteria: the solver reproduces the no-hazard case exactly; closed-edge feasibility is enforced; hazard-state costs are finite and monotone where required; all scenario profiles pass VI and omitted-route certificates.

### Phase 3: Test state-contingent coalition regimes

Purpose: determine when singleton, partial, or grand operational coordination is best.

Cross hazard severity, capacity loss, closure topology, AV order load, order exposure to critical corridors, objective composition, governance intensity, and governance burden. Use full partition audits for n=4 as an audit baseline.

Exit criteria: each claim maps to a raw profile file, summary file, validator, table, and figure; best partitions are reported as sets under a declared tolerance.

### Phase 4: Reference-curvature and behavioral robustness

Purpose: test whether the high-AV-load regime is an artifact of calibration or HDV behavior.

Compare reference slopes from hazard-state HDV-only equilibrium, no-coordination equilibrium, and scaled versions of each. Add delayed/publicly incomplete road information and bounded-rational or logit HDV sensitivity as secondary models.

Exit criteria: report transition stability, partition regret, and effect-size ranges over reference calibrations and HDV information cases.

### Phase 5: Response-aware policy selection

Purpose: connect the analytical merge-direction score to the actual algorithm.

Implement local score estimation on a fixed support, rank candidate merges by score, then validate selected merges with the certified continuation solver. Report score sign accuracy, ranking correlation, finite-regret, support-change failures, and equilibrium solves saved.

Exit criteria: no score-based claim is made without a finite-endpoint audit; support changes are explicitly detected and reported.

### Phase 6: Scale and limited empirical grounding

Purpose: test computational value and real-world plausibility.

Run n=6 exhaustive audits where feasible, n=8 reduced audits, and n>=10 heuristic searches. Validate network speeds/closures against public LA data. Treat AV order portfolios and objective types as calibrated scenarios unless confidential data become available.

Exit criteria: report runtime, memory, partition coverage, and regret; distinguish empirical calibration from synthetic scenario generation.

## Main claims allowed after completion

- Wildfire state changes the value and ranking of temporary AV operational coalitions.
- Partial coordination can be preferred to singleton or grand coordination under heterogeneous order exposure, hazard topology, or governance burden.
- Coalition structure is more informative than HHI for emergency network performance.
- The response-aware policy is useful only conditionally on fixed-support diagnostics and has measurable regret outside that condition.

Claims not allowed without additional work: improved real evacuation clearance, lives saved, universal optimality, voluntary firm stability, or generalization to all LA operators.
