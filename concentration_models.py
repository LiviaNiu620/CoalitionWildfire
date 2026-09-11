"""
Two-link (and 3-link) minimal models for "concentration of AV control".

Goal: test whether reorganizing AV concentration (group mass m = alpha/k) can
break the selfish-HDV "filling" plateau, and whether system cost J(m) is U-shaped.

We compare THREE mechanisms that differ in how HDV filling is (im)perfect:
  MECH 1: both links congestible (no perfect escape road)
  MECH 2: elastic / switching-cost HDV demand (filling not instantaneous)
  MECH 3: weaving bottleneck coupling (interaction term A<->B)  [added as variant]

For each mechanism we compute the mixed atomic-nonatomic equilibrium by solving
the coupled best-response conditions as a variational inequality via projected
fixed-point iteration with careful interior handling.

Player structure (symmetric):
  HDV: nonatomic mass (1-alpha), selfish Wardrop.
  AV : mass alpha in k equal atomic groups of mass m; each group internalizes its
       OWN marginal congestion with intensity tau in [0,1]:
         perceived marginal cost on link L = c_L(x_L) + tau * m * c_L'(x_L)
  Everyone routes on the same 2 links (A,B). x_L = total flow on L.
"""

import numpy as np


# ============================================================
#  Mechanism 1: both links congestible
#     c_A(x) = a1*x
#     c_B(x) = b0 + b1*x
# ============================================================
class Mech1:
    name = "both-congestible"
    a1 = 1.0
    b0 = 0.3
    b1 = 0.5

    @classmethod
    def cA(cls, x):  return cls.a1 * x
    @classmethod
    def cAp(cls, x): return cls.a1
    @classmethod
    def cB(cls, x):  return cls.b0 + cls.b1 * x
    @classmethod
    def cBp(cls, x): return cls.b1


# ============================================================
#  Mechanism 3 (variant): weaving coupling
#     Using link B cost depends also on x_A (conflict at merge)
#     We keep it simple: c_A = a1*x_A ; c_B = b0 + b1*x_B + w*x_A
#  (HDV filling A raises B's cost -> filling has a self-limiting externality)
# ============================================================
class Mech3:
    name = "weaving-coupled"
    a1 = 1.0
    b0 = 0.3
    b1 = 0.5
    w  = 0.4   # weaving coupling

    @classmethod
    def cA(cls, xA, xB):  return cls.a1 * xA
    @classmethod
    def cB(cls, xA, xB):  return cls.b0 + cls.b1 * xB + cls.w * xA


TOTAL = 1.0


# ------------------------------------------------------------
#  Generic lower-level equilibrium for MECH1-style (separable) costs
#  Variables: fraction on A for HDV (hA in [0,1-alpha]) and per-group aA in [0,m].
#  Symmetric groups.  Solve by projected fixed point on marginal-cost gaps.
# ------------------------------------------------------------
def equilibrium_separable(M, alpha, k, tau, iters=20000, lam=0.02):
    hdv = 1.0 - alpha
    m = alpha / k if k > 0 else 0.0
    hA = 0.5 * hdv
    aA = 0.5 * m
    for _ in range(iters):
        xA = hA + k * aA
        xB = TOTAL - xA
        # HDV selfish gradient: cA - cB  (>0 => move off A)
        gH = M.cA(xA) - M.cB(xB)
        # AV group internalized gradient
        gA = (M.cA(xA) + tau * m * M.cAp(xA)) - (M.cB(xB) + tau * m * M.cBp(xB))
        # gradient descent on fractions (move away from higher-cost link)
        hA = np.clip(hA - lam * gH, 0.0, hdv)
        aA = np.clip(aA - lam * gA, 0.0, m)
    xA = hA + k * aA
    xB = TOTAL - xA
    J = xA * M.cA(xA) + xB * M.cB(xB)
    return dict(m=m, tau=tau, hA=hA, aA=aA, xA=xA, xB=xB, J=J,
                cA=M.cA(xA), cB=M.cB(xB))


def refs_separable(M, alpha):
    # all selfish Wardrop: cA(xA)=cB(1-xA)
    # a1*xA = b0 + b1*(1-xA) -> xA*(a1+b1) = b0+b1 -> xA
    xA = (M.b0 + M.b1) / (M.a1 + M.b1)
    xA = min(max(xA, 0.0), 1.0)
    xB = 1 - xA
    J_wardrop = xA * M.cA(xA) + xB * M.cB(xB)
    # system optimal: min xA*cA + (1-xA)*cB
    # d/dxA [a1 xA^2 + (b0+b1(1-xA))(1-xA)] = 0
    # = a1*... let's just grid
    grid = np.linspace(0, 1, 100001)
    Js = grid * M.cA(grid) + (1 - grid) * M.cB(1 - grid)
    i = np.argmin(Js)
    return dict(wardrop_J=J_wardrop, wardrop_xA=xA,
                opt_J=Js[i], opt_xA=grid[i])


# ------------------------------------------------------------
#  Endogenous tau via obedience (for separable mechanisms)
#  A sacrificed AV (pushed toward the high-private-cost link) can deviate.
#  Temptation = private cost gap it would gain by deviating.
#  Enforceable wedge = tau*m*c' (the internalization the group sustains).
#  Obedience holds iff enforceable >= temptation. Platform picks largest feasible
#  tau that ALSO minimizes J (we scan tau on a grid under the obedience filter).
# ------------------------------------------------------------
def endogenous_tau_and_eq(M, alpha, k, tau_grid=None):
    if tau_grid is None:
        tau_grid = np.linspace(0.0, 1.0, 51)
    m = alpha / k if k > 0 else 0.0
    best = None
    tau_bar = 0.0  # largest obedience-feasible tau
    for tau in tau_grid:
        eq = equilibrium_separable(M, alpha, k, tau)
        xA, xB = eq["xA"], eq["xB"]
        # which link is the group being pushed toward on the margin?
        # temptation for a vehicle asked to take the less-private-attractive link:
        # private gap between the two links
        priv_gap = abs(M.cA(xA) - M.cB(xB))
        enforceable = tau * m * max(M.cAp(xA), M.cBp(xB))
        obedient = enforceable + 1e-9 >= priv_gap
        if obedient:
            tau_bar = max(tau_bar, tau)
            if best is None or eq["J"] < best["J"]:
                best = eq
    if best is None:
        # nothing obedient beyond tau=0 -> group acts selfish (tau=0)
        best = equilibrium_separable(M, alpha, k, 0.0)
        tau_bar = 0.0
    best = dict(best)
    best["tau_bar"] = tau_bar
    return best


if __name__ == "__main__":
    for M in [Mech1]:
        print("=" * 60)
        print("MECHANISM:", M.name)
        alpha = 0.4
        r = refs_separable(M, alpha)
        print(f"alpha={alpha}")
        print(f"  Wardrop(all selfish): xA={r['wardrop_xA']:.4f}  J={r['wardrop_J']:.5f}")
        print(f"  System optimal      : xA={r['opt_xA']:.4f}  J={r['opt_J']:.5f}")
        print()
        print("  --- exogenous tau=1, sweep concentration k ---")
        print("   k      m      xA      xB       J       (vs Wardrop)")
        for k in [1, 2, 3, 4, 6, 8, 16, 32, 128, 1024]:
            eq = equilibrium_separable(M, alpha, k, tau=1.0)
            print(f"  {k:5d} {eq['m']:.4f} {eq['xA']:.4f} {eq['xB']:.4f}  {eq['J']:.5f}   {eq['J']-r['wardrop_J']:+.5f}")
        print()
        print("  --- endogenous tau(m) via obedience, sweep k ---")
        print("   k      m    tau_bar   xA       J       (vs Wardrop)")
        for k in [1, 2, 3, 4, 6, 8, 16, 32, 128, 1024]:
            eq = endogenous_tau_and_eq(M, alpha, k)
            print(f"  {k:5d} {eq['m']:.4f}  {eq['tau_bar']:.3f}   {eq['xA']:.4f}  {eq['J']:.5f}   {eq['J']-r['wardrop_J']:+.5f}")
