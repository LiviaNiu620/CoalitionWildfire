# Coalition-design prototype (small networks only)

numpy/scipy prototype supporting the revision plan. Not a replacement for the Sioux Falls solver.

- core.py      : affine potential-game equilibrium as convex QP (active-set with KKT certification);
                 adjoint merge scores `pair_scores` (one linear solve -> all pairwise m_ij); partition utilities.
- nets.py      : affine Braess (paper's parameters) and a 3x3 directed grid (3 ODs).
- bpr.py       : BPR networks; 'frozen' (bbar = c'(x_UE), paper's model) vs 'self' (damped diagonalization,
                 bbar -> c'(x*), i.e. true state-dependent atomic internalization); reports TRUE-model VI gap.
- so.py        : system-optimum benchmark (BPR).
- policies2.py : fixed-support prediction (Thm-1 style) and predicted two-step policy.
- test1.py     : adjoint vs finite-difference check.
- exp_n4.py    : n=4 audit; one-step / two-step / score policies.
- exp_bpr.py, recov.py : frozen vs self-consistent slope robustness (Q2); regret in recovery points.
- exp_n6.py, exp_n8.py, exp_n8b.py, exp_pred.py : scalability (solve counts) for n=6/8.
- exp_homog.py, contig.py : homogeneous vs heterogeneous tau; structure of partial optima.
- *.json / *.log : raw outputs used in the revision document.

Regret units: exp_n4 = relative TSTT; exp_n6/n8/pred = fraction of (J_singleton - J*) gap; recov.py = recovery points.
