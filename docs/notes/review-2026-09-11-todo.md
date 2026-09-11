# Coalition Design paper — review findings & TODO (2026-09-11)

## Key findings
- PDF citations all render as "?" and References is empty — bib not compiled.
- Notation clashes: η (BPR exponent vs merge path), B (0.15 vs slope matrix), N, H, K overloaded.
- Current model is an exact potential game: F_Π = ∇[Σ_e ∫c_e + Σ_C Ψ_C] (sympy-verified symmetric Jacobian with frozen b̄). Text claims "not assumed to possess a common potential" — must fix. True atomic splittable with BPR is NOT symmetric; frozen b̄_e = c'_e(x^UE) may drive the grand→partial transition at α=0.9 → needs robustness check.
- Information: implicit complete information (players observe x^{-C}; b̄ uses global UE slopes; within-coalition full sharing; deterministic Wardrop HDVs; designer knows all τ, OD, HDV support). Add explicit Information-structure assumption (steady-state interpretation) + designer-misspecification regret experiment.
- Theory is local (fixed support, affine, 2-link bottleneck); merge score M_Π not used by the one-/two-step policies (re-solve equilibria instead).
- Experiments: full audit only n=4; SF regime claim from 4 states; Braess/3x3 grid toy; effect sizes small; always-grand exact in 27/32 vs one-step 29/32.

## TODO
Phase A (1–2 wk): fix bib/notation/redundancy/abstract ≤250 words; information assumption; potential-game proposition; cite Hayrapetyan 2006, Huang 2013, Cominetti 2009, Battifarano & Qian TS 2023, Toso et al. 2024.
Phase B (3–5 wk): score-guided merge policy (sign accuracy, regret, # solves saved); true-atomic-slope robustness; n=6/8 exhaustive + n≥10 heuristics; dense (α,γ,demand) phase map; one large network (Anaheim/Chicago Sketch).
Phase C (2–4 wk): richer objectives (q_C, SO-seeking operator); IR/core/Shapley + AV/HDV distribution; designer τ-misspecification; attempt a global theorem.

## Venue
Primary: Transportation Science (or TR-B) after Phase A+B. TAC only if rewritten as theory paper (global PoA/monotonicity/convergence results). T-ITS as fallback (compress, emphasize scale/application).
