# Revision plan + prototype findings (2026-09-11)

Prototype code: prototype/ (numpy/scipy; Braess + 3x3 grid; n=4/6/8).

## Key findings
- Frozen slope b̄=c'(x^UE) vs self-consistent (true atomic, damped diagonalization b̄←c'(x*)): in 5/12 toy BPR states optimal partition differs; in every case frozen picks partial (K=2) while true picks grand; frozen choice regret under true behavior 2.5–12.5 recovery pp. => Must run V0/V1/V2/V3/V4 robustness on Sioux Falls before claiming grand→partial at high α.
- Adjoint merge scores: one adjoint solve gives all pairwise m_ij; M_{A,B}=Σ_{i∈A,j∈B} m_ij; verified vs finite difference (4–5 digits).
- Scalability (solve counts): n=6 exhaustive 203 / one-step 36 / two-step 114 / predicted two-step (k=2) 13 (exact at γ=0.5); n=8 exhaustive 4140 / one-step 85 / two-step 533 / predicted two-step 19–20 (exact at γ=0.5).
- γ=1 issue: Ω_{C,e}(1)=b̄√τ√τᵀ is rank one → tangent-curvature condition (19) fails for multi-member coalitions → Thm 1/P_C undefined at γ=1; first-order/fixed-support predictions unreliable there (sign acc 0.2–0.6). Fix: treat zero-curvature coalition directions as Wardrop-like and add to screening subspace L=[L_H, L_0].
- Homogeneous τ=1 (affine grid n=6, 18 states): grand always optimal; heterogeneous τ produces partial optima → novelty: objective heterogeneity (need some τ>1 in bottleneck: block contributes m/(1+γ(m−1)) ≥1 when τ=1) drives partial dominance, unlike collusion literature (Hayrapetyan 2006, Bhaskar–Fleischer–Huang 2010, Huang 2013).

## Plan order
Wk1–2: potential-game statement, γ=1 theory fix, diagonalization, Q2 V0/V1/V4 on SF. Wk3–4: adjoint scores + predicted policy, n=6/8. Wk5–6: G-a (homogeneous benchmark), G-b (globalized merge identity), optional G-c (externality-alignment bound). Wk6–9: multi-dim objectives (τ, q, SVO, ODD), Anaheim, dense states & random instances.
