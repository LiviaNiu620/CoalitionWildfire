# Prototype experiment index (2026-09-11)

All runs are small-network prototypes (numpy/scipy; sympy for symbolic checks). Every equilibrium is
KKT/VI-certified (residual ~1e-12). Regret units differ by script (see last column).

| # | Question | Script | Output | Key result |
|---|---|---|---|---|
| 1 | Bottleneck derivative (Thm 2); potential-game structure | `symbolic_checks.py` (sympy) | `symbolic_checks.log` | residual 0; frozen-slope Jacobian symmetric (exact potential game); true atomic Jacobian asymmetric: 9/5 (f1-f2) x^2 |
| 2 | Adjoint merge score vs finite difference | `test1.py` | `test1.log` | agree to 4-5 significant digits (Braess, grid; alpha=0.5, 0.9) |
| 3 | n=4 audit, affine Braess+grid, 40 states; one-step / two-step / score / surrogate / grand; score sign accuracy | `exp_n4.py` | `exp_n4.json` | grand mostly optimal; screened regime (Braess alpha=0.3) makes first-order scores 0 -> score policies stop too early (regret = relative TSTT) |
| 4 | Frozen b̄=c'(x_UE) vs self-consistent b̄=c'(x*) (true atomic), BPR, 12 states | `exp_bpr.py` (+ `bpr.py`) | `exp_bpr.json` | 5/12 states differ; always frozen->partial, true->grand |
| 4b | Regret of frozen-model choice under true behaviour, in recovery points | `recov.py` (+ `so.py`) | `recov.log` | 2.5-12.5 recovery pp |
| 5 | n=6 random heterogeneous instances (24 states), exhaustive 203 partitions vs policies | `exp_n6.py` | `exp_n6.json` | one-step 36 solves (22/24 exact), two-step 114 (23/24), surrogate+refine 8-18 (23/24); regret = fraction of (J_singleton - J*) |
| 6 | Homogeneous tau=1 vs heterogeneous tau (n=6, 18 states each) | `exp_homog.py` | `homog.log` | homogeneous: grand always optimal; heterogeneous: partial optima appear |
| 7 | Structure of partial optima (contiguity in tau order) | `contig.py` | `contig.log` | 2 of 3 partial optima contiguous in tau |
| 8 | n=8 exhaustive (4140 partitions), 6 states | `exp_n8.py` | `exp_n8.json`, `n8.log` | one-step 85 solves (4/6), two-step 533 (6/6), surrogate+refine 10 (4/6) |
| 9 | First-order-score screened two-step, n=8 gamma=1 hard states | `exp_n8b.py` (+ `policies2.py`) | `n8b.log` | fails at gamma=1 (regret 0.003-0.65) |
| 10 | Fixed-support-prediction two-step (Thm-1 style), n=6 and n=8 | `exp_pred.py` (+ `policies2.py`) | `exp_pred_n6.json`, `pred6.log`, `pred8.log`, `pred8g05.log`, `exp_pred_n8.json` | gamma=0.5: exact with 13 (n=6) / 19-20 (n=8) solves; gamma=1: sign accuracy 0.2-0.6 |
| 11 | Tangent-curvature condition (19) for a merged 2-member coalition | `gamma1_curvature.py` | `gamma1_curvature.log` | min eigenvalue of Z^T H_C Z ~ 0 (route-space degeneracy; rank-one Omega at gamma=1) |

Notes
- `pred8.log` is a partial run (gamma=1, interrupted after 1 state); `exp_pred_n8.json` holds the gamma=0.5 run.
- Networks: `nets.py` (affine Braess with the paper's parameters; 3x3 directed grid, 3 ODs); `bpr.py` (BPR versions).
