# Wildfire Coalition Model Specification

## Scope lock

This is a coalition-design paper, not a household-evacuation or rescue paper. HDVs remain ordinary human-driven vehicles. Wildfire enters through an exogenous emergency network state and public information state. Commercial AV operators keep accepted order ownership and pickup-dropoff OD obligations.

## Objects and authority

The regulator selects a temporary operational partition \(\Pi\) and coordination intensity \(\gamma\). It does not select individual vehicle routes, reassign orders, or impose legal ownership changes. A coalition must have joint routing/dispatch authority or shared routing governance; otherwise a partition is behaviorally inert.

## State

Use a scenario or snapshot state

\[
s_t^\omega=(E_t^\omega,t_{t}^{0,\omega},\kappa_t^\omega,r_t^\omega,\mathcal I_t^\omega,d_t^A,d_t^H,E_{\mathrm{critical}}^\omega).
\]

The first implementation uses static snapshots. A rolling horizon is a later extension.

## Flow and demand

\[
x_{e,t}=u_{e,t}^\omega+h_{e,t}+\sum_i f_{e,t}^i.
\]

\(h\) is ordinary HDV traffic. \(d^H\) is ordinary HDV demand. \(d^i_{w,t}\) is accepted commercial AV order demand for company \(i\) in a short decision epoch.

## Coalition feasibility

\[
\mathcal Y_C^{\mathrm{id}}=\{y^C\geq0:\Gamma_Cy^C=d^C\}.
\]

Do not use operational pooling in the main treatment. Pooling changes service authority and belongs in an extension.

## Coalition continuation objective

\[
\Phi_C(y^C;z_{-C},s)
=\sum_e\int_0^{F_e^C}c_e^s(x_{-C,e}+u)\,du
+\frac12(y^C)^\top H_C(\gamma)y^C+q_C^\top y^C.
\]

Under symmetric PSD fixed-reference governance, the VI is the gradient condition of

\[
\mathcal P_\Pi^s(y,h)
=\sum_e\int_0^{x_e}c_e^s(u)\,du
+\sum_{C\in\Pi}\Psi_C(y^C).
\]

The exact activated-response theorem remains affine and support-conditional. BPR experiments are nonlinear continuation experiments; fixed-reference and state-dependent slope results must be reported separately.

## HDV information

Baseline HDVs observe public road travel-time and closure information and respond through deterministic Wardrop conditions. They do not observe private company objectives or private order portfolios and are not directly routed. Information delay and bounded rationality are robustness treatments.

## Regulator objective

\[
L^s(\Pi)=J_{\mathrm{TSTT}}^s(\Pi)+\lambda_AJ_{\mathrm{order}}^s(\Pi)+\lambda_RJ_{\mathrm{risk}}^s(\Pi)+\lambda_CJ_{\mathrm{critical}}^s(\Pi)+\lambda_GC_{\mathrm{gov}}(\Pi,\gamma).
\]

\(J_{\mathrm{critical}}\) is an emergency-sensitive corridor penalty or capacity-protection metric, not an explicit evacuee-flow objective.

## Assumption ledger

1. Wildfire state is exogenous to the routing continuation solve.
2. Accepted orders are fixed within a short decision epoch.
3. Company OD obligations and order ownership are preserved.
4. Temporary operational coalition authority is available as an institutional scenario.
5. HDVs are ordinary informed but uncontrolled users.
6. Physical costs are common across vehicle classes.
7. Fixed-reference governance is symmetric PSD in the primary potential benchmark.
8. Reference curvature is a declared calibration parameter, not an unobserved implementation detail.
9. The primary partition is regulator-selected, not a voluntary stable formation outcome.
10. Real wildfire and operator data are used for calibration/validation only when available; otherwise scenarios are labeled synthetic or calibrated.
