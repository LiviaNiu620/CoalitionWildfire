# Introduction Architecture

## Paragraph 1: Wildfire as a network disruption

Open with wildfire-induced road closures, capacity loss, risk, and rapidly changing public road information. Do not present a new fire-spread model and do not introduce AVs as rescue vehicles.

## Paragraph 2: Commercial AV orders remain a traffic flow

Explain that accepted commercial pickup-dropoff orders may remain active within a short emergency decision epoch. These orders continue to consume scarce road capacity. Their OD and company ownership are preserved.

## Paragraph 3: Programmability is fragmented

Ordinary HDVs observe public road information and reroute but cannot be directly controlled. AVs are centrally dispatchable within firms, yet authority is fragmented across operators. This creates an institutional rather than vehicle-level control problem.

## Paragraph 4: Regulator-mandated operational coalitions

The regulator cannot reassign orders or prescribe individual routes. It can mandate temporary participation in an operational coalition that changes joint routing authority and cross-operator internalization. State clearly that this is not a legal merger or voluntary coalition-formation equilibrium.

## Paragraph 5: Why grand coordination is not automatic

Grand coordination may internalize more cross-company congestion but can be poorly matched to heterogeneous order portfolios, objective types, route access, hazard exposure, HDV adjustment, and governance burden. The scientific question is whom to coordinate and how strongly, not whether coordination is always beneficial.

## Paragraph 6: Mechanism and gap

Introduce the chain: governance changes an identity-preserving AV response; the network maps it into edge flow; ordinary HDVs reroute; the residual activated response affects congestion. Explain that ownership concentration and coalition size omit this mechanism.

## Paragraph 7: Research questions

1. When does a temporary AV coalition improve network and critical-corridor performance under wildfire disruption?
2. Which operators should be grouped when objectives and order-corridor exposure differ?
3. Can a response-aware local policy select a low-regret partition without enumerating the full partition lattice?

## Paragraph 8: Contributions

1. A hazard-conditioned, identity-preserving coalition continuation model with ordinary informed but uncontrolled HDVs.
2. A corrected fixed-reference convex-potential characterization and support-conditional activated-response mechanism.
3. Hazard-state, reference-curvature, information-error, and order-exposure experiments that identify singleton, partial, and grand regimes without treating HHI as causal.
4. A merge-score-guided policy evaluated by sign accuracy, finite regret, support-change fallback, equilibrium solves, and runtime.

## Scope sentence

The paper does not model household evacuation demand, rescue dispatch, shelter assignment, fire-spread dynamics, voluntary coalition stability, pricing, new-order acceptance, or inter-epoch empty-vehicle repositioning.
