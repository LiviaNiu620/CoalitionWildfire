"""Check the tangent-curvature condition (19): min eigenvalue of Z^T H_C Z for a merged 2-member coalition."""
import numpy as np
from scipy.linalg import null_space
from core import E_of_partition
from nets import grid
for g in [0.5, 0.9, 1.0]:
    inst = grid(0.9, [0.5, 1.5], g)
    H = inst.hessian(E_of_partition([[0, 1]], 2))
    R = inst.R
    Hc = H[:2*R, :2*R] - inst.H_road[:2*R, :2*R]      # management curvature H_C only
    Z = null_space(inst.A[:2*inst.W, :2*R])
    ev = np.linalg.eigvalsh(Z.T @ Hc @ Z)
    print('gamma %.1f  min eig(Z^T H_C Z) = %.3e   max = %.3e' % (g, ev.min(), ev.max()))
