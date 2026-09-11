# Coalition Design for Commercial AV Traffic in Wildfire-Disrupted Networks

## Project summary

Wildfire creates a temporary traffic state with road closures, capacity loss, risk, and incomplete public information. Commercial autonomous-vehicle (AV) orders do not necessarily disappear during that state. Each operator retains accepted pickup-dropoff obligations, while ordinary human-driven vehicles (HDVs) remain responsive users that can observe public road information but cannot be directly routed by the regulator. AV routing is programmable at the company level, yet routing authority is fragmented across operators.

This project studies whether a regulator can reduce the congestion externality of commercial AV orders by mandating temporary operational coalitions. The regulator selects which operators must share routing authority and the intensity of cross-operator coordination. It does not reassign passenger orders, change pickup-dropoff OD pairs, or prescribe vehicle-level routes. Coalition members continue to route their own accepted orders subject to identity-preserving OD constraints. The central object is therefore coalition design, not wildfire spread, rescue dispatch, or household evacuation modeling.

## Research question

When wildfire-induced network disruption changes road accessibility and public information, which regulator-mandated temporary coalition structure among commercial AV operators minimizes network congestion, emergency-critical corridor overload, commercial order delay, and coordination burden?

## Research hypotheses

H1. Coalition structure changes network performance because joint routing authority changes the aggregate AV response even when order ownership and OD obligations are fixed.

H2. The preferred coalition count is state-dependent. Grand coordination is not uniformly optimal because order exposure, objective heterogeneity, hazard topology, HDV rerouting, and governance burden can favor partial coordination.

H3. Company mass or HHI is not a sufficient emergency coordination state. Holding mass concentration fixed, order exposure to critical corridors and member objective composition can change the preferred partition.

H4. A fixed reference-curvature calibration can alter high-AV-load rankings. State-dependent slope and reference-state robustness are therefore required before interpreting a grand-to-partial transition.

H5. A response-aware merge score can reduce equilibrium solves, but it is a fixed-support screening tool rather than a global optimality guarantee.

## Institutional setting

The regulator is assumed to have temporary emergency coordination authority. This authority permits mandatory participation in an operational coordination protocol, data exchange, and joint dispatch governance. It does not imply legal ownership consolidation or unrestricted command over private vehicles. The authority is an explicit institutional scenario; its jurisdictional validity must be checked separately against California and local emergency-management rules.

The study uses short decision epochs. Within an epoch, accepted commercial orders and ordinary HDV demand are fixed. Across epochs, new orders, cancellations, and state updates are future extensions. No AV is assigned a rescue or household-evacuation task.

## Mathematical model

For hazard state \(\omega\) and decision epoch \(t\), the available network is \(G_t^\omega=(V,E_t^\omega)\). Closed edges are removed. Open but degraded edges have state-dependent capacity \(\kappa_{e,t}^\omega\), free-flow time \(t_{e,t}^{0,\omega}\), and optional risk attribute \(r_{e,t}^\omega\). Physical cost is

\[
c_{e,t}^\omega(x_e)=t_{e,t}^{0,\omega}\left[1+B\left(x_e/\kappa_{e,t}^\omega\right)^4\right].
\]

Total edge flow is ordinary HDV flow plus commercial AV order flow and, if used, exogenous background flow:

\[
x_{e,t}=u_{e,t}^\omega+h_{e,t}+\sum_i f_{e,t}^i.
\]

Company \(i\)'s accepted order demand is \(d_{w,t}^i\), with fixed pickup-dropoff OD. A coalition \(C\) uses the identity-preserving feasible set

\[
\mathcal Y_{C,t}^{\mathrm{id}}=\{y_t^C\geq0:\Gamma_{C,t}y_t^C=d_t^C\}.
\]

The coalition objective retains the existing structure:

\[
\Phi_{C,t}(y^C;z_{-C})
=\sum_e\int_0^{F_{C,e}}c_{e,t}^\omega(x_{-C,e}+s)\,ds
+\Psi_{C,t}(y^C),
\]

with a symmetric PSD quadratic governance term. Under fixed-reference curvature, the continuation equilibrium is the minimizer of

\[
\mathcal P_{\Pi,t}^\omega
=\sum_e\int_0^{x_{e,t}}c_{e,t}^\omega(u)\,du
+\sum_{C\in\Pi}\Psi_{C,t}(y_t^C).
\]

HDVs use public experienced road costs and optional public risk penalties. They are informed but not directly controllable. The regulator selects \(\Pi_t\) and coordination intensity \(\gamma_t\), then evaluates the continuation equilibrium using

\[
L_t^\omega(\Pi_t)=J_{\mathrm{TSTT}}
+\lambda_AJ_{\mathrm{order}}
+\lambda_RJ_{\mathrm{risk}}
+\lambda_CJ_{\mathrm{critical}}
+\lambda_GC_{\mathrm{gov}}.
\]

## Work packages

### WP1: Model and solver integrity

Correct the potential statement, separate general VI claims from the fixed-reference potential special case, document the information structure, and retain the original solver as a regression baseline. Add the state-dependent slope diagnostic without changing the default reference-slope results.

### WP2: Hazard-state network adapter

Implement static hazard snapshots with edge closure, capacity degradation, free-flow-time changes, risk values, public-information state, and critical-edge metadata. The no-hazard adapter must reproduce the existing solver output.

### WP3: Coalition regime experiments

Run complete n=4 partition audits across hazard severity, closure topology, AV order load, order-corridor exposure, governance intensity, and objective composition. Report best partition sets, best coalition count, TSTT, order delay, critical-edge overload, risk-weighted use, and governance cost.

### WP4: Robustness and policy selection

Compare reference-state calibrations and state-dependent slopes; perturb public road information and regulator estimates; implement fixed-support merge-score ranking with certified finite-endpoint fallback.

### WP5: Scale and empirical grounding

Run n=6 exhaustive and n=8 selected-state audits, then n>=10 heuristic searches. Calibrate hazard scenarios and traffic states with public Los Angeles data. Treat AV order portfolios and objective types as controlled or calibrated scenarios unless operator data are obtained.

## Data plan

Network: OpenStreetMap, Caltrans GIS, LADOT data, and Caltrans PeMS.

Wildfire: CAL FIRE incident/perimeter records, NIFC archives, NASA FIRMS, NOAA wind/weather, Caltrans QuickMap closures, and local emergency GIS.

Traffic: PeMS, LADOT sensors, and, subject to license, StreetLight, Replica, or INRIX.

Commercial AV: ideal data are company-level accepted-order OD, timestamps, delays, cancellations, service zones, and fleet exposure. These are likely proprietary. Without them, use calibrated synthetic OD portfolios and state the limitation explicitly.

Institutional: California and local emergency coordination authorities, AV permit conditions, public-private data-sharing agreements, and emergency traffic-management plans.

## Evaluation and deliverables

Deliverables are a validated hazard adapter, raw and summary profile files, fail-closed validators, hazard regime tables and figures, reference-slope robustness results, information-error results, score-policy runtime/regret results, and a reproducible model/experiment ledger.

The project will not claim improved evacuation clearance, lives saved, household rescue, or universal coalition optimality because household evacuation demand and rescue operations are outside scope.

## Risks and mitigation

If emergency coalition authority is not legally available, interpret the intervention as a proposed temporary coordination protocol and report the institutional assumption. If commercial AV order data are unavailable, separate calibrated scenario analysis from empirical validation. If state-dependent BPR response is non-monotone, retain the fixed-reference potential as the primary benchmark, report multi-start diagnostics, and avoid uniqueness claims.
