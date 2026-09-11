"""
Two-link minimal model for "concentration of AV control" study.

Network: single OD, two parallel links.
    Link A: c_A(x) = a_coef * x        (congestible)
    Link B: c_B(x) = b_const           (constant "escape" road)

Populations:
    HDV: mass (1-alpha), nonatomic, selfish  -> Wardrop
    AV : mass alpha, split into k equal groups of mass m = alpha/k.
         Each group is an ATOMIC player that internalizes its own marginal
         congestion with intensity tau in [0,1] (transfer-free market power / dividend).

Perceived marginal cost of link A for an AV group of mass m:
    mc_A^group = c_A(x_A) + tau * m * c_A'(x_A)
where x_A is TOTAL flow on A (HDV + all AV groups).

We solve the mixed atomic-nonatomic equilibrium at the lower level, then
sweep concentration (k) to get system cost J(m).

Two versions of tau:
    (1) exogenous tau = 1  (upper bound; isolates market-power badness)
    (2) endogenous tau(m) from an obedience constraint (the real deal; tests Hard Point 3)
"""

import numpy as np

# ---------------- primitives ----------------
A_COEF = 1.0      # c_A(x) = A_COEF * x
B_CONST = 0.6     # c_B(x) = B_CONST
TOTAL_DEMAND = 1.0


def cA(x):
    return A_COEF * x

def cA_prime(x):
    return A_COEF

def cB(x):
    return B_CONST


# ---------------- lower-level equilibrium ----------------
# Unknowns: hA (HDV mass on A), and for a SYMMETRIC configuration all groups
# behave identically, so let aA = mass on A from EACH group (0..m). Total AV on A = k*aA.
#
# HDV Wardrop: nonatomic, selfish.
#   if hA in (0, 1-alpha): c_A(xA) == c_B  (indifferent)
#   corner cases handled.
# AV group best response (symmetric atomic):
#   group compares perceived mc on A vs B, given total flow xA.
#   perceived mc_A = c_A(xA) + tau*m*c_A'(xA);  mc_B = c_B
#
# We find fixed point by solving for the common congestion level.

def lower_equilibrium(alpha, k, tau):
    """Return dict with flows and total system cost J for symmetric config."""
    m = alpha / k if k > 0 else 0.0

    # We search over xA (total flow on link A) as the single scalar that both
    # populations respond to, then check consistency.
    # HDV chooses A iff c_A(xA) < c_B ; indifferent iff equal.
    # AV group chooses A-fraction so that perceived mc balances (interior) or corner.

    # Strategy: iterate on (hA, aA_each).
    # Use fixed-point iteration.
    hdv_mass = 1.0 - alpha

    hA = 0.5 * hdv_mass
    aA_each = 0.5 * m
    for _ in range(2000):
        xA = hA + k * aA_each

        # --- HDV best response (nonatomic, selfish) ---
        # wants A if cA(xA) < cB, B if >, indifferent if ==
        cA_val = cA(xA)
        if cA_val < B_CONST - 1e-12:
            hA_new = hdv_mass
        elif cA_val > B_CONST + 1e-12:
            hA_new = 0.0
        else:
            hA_new = hA  # indifferent, keep (pins via AV side)

        # --- AV group best response (atomic, internalizing) ---
        # perceived mc on A = cA(xA) + tau*m*cA'(xA); compare to cB
        # NOTE: group takes others' flow as given; use xA but its own marginal term uses m
        perc_mcA = cA(xA) + tau * m * cA_prime(xA)
        if perc_mcA < B_CONST - 1e-12:
            aA_new = m
        elif perc_mcA > B_CONST + 1e-12:
            aA_new = 0.0
        else:
            aA_new = aA_each

        # damping
        lam = 0.1
        hA = (1 - lam) * hA + lam * hA_new
        aA_each = (1 - lam) * aA_each + lam * aA_new

        # clip
        hA = min(max(hA, 0.0), hdv_mass)
        aA_each = min(max(aA_each, 0.0), m)

    xA = hA + k * aA_each
    xB = TOTAL_DEMAND - xA
    # system total travel cost
    J = xA * cA(xA) + xB * cB(xB)
    return dict(m=m, hA=hA, aA_each=aA_each, xA=xA, xB=xB, J=J,
                cA=cA(xA), cB=B_CONST, tau=tau)


# ---------------- reference costs ----------------
def cost_all_selfish_baseline(alpha):
    """Everyone (HDV + AV) selfish nonatomic -> Wardrop. tau irrelevant, k->inf equiv."""
    # pure Wardrop: cA(xA)=cB => xA = B_CONST/A_COEF, capped at demand
    xA = min(B_CONST / A_COEF, TOTAL_DEMAND)
    xB = TOTAL_DEMAND - xA
    J = xA * cA(xA) + xB * cB(xB)
    return dict(xA=xA, xB=xB, J=J)


def cost_system_optimal():
    """Full social optimum (single planner over everyone)."""
    # min over xA of xA*cA(xA) + (1-xA)*cB
    # = A*xA^2 + B_CONST*(1-xA); derivative 2A*xA - B_CONST = 0 -> xA = B/(2A)
    xA = min(max(B_CONST / (2 * A_COEF), 0.0), TOTAL_DEMAND)
    xB = TOTAL_DEMAND - xA
    J = xA * cA(xA) + xB * cB(xB)
    return dict(xA=xA, xB=xB, J=J)


# ---------------- endogenous tau via obedience ----------------
def endogenous_tau(alpha, k):
    """
    Obedience-constrained internalization.

    The platform wants tau as high as possible (=1) to internalize, but a
    'sacrificed' AV (recommended to B while A looks cheaper privately) can deviate.

    Private temptation to deviate from B back to A:  Delta = c_B - c_A(xA)  (if >0)
    The benefit of staying that flows back to *itself* scales with the marginal
    congestion relief it personally enjoys, ~ (its share) of tau*m*c_A' effect.

    Simplified obedience: the group can enforce internalization intensity tau only
    if the per-vehicle temptation does not exceed the per-vehicle share of the
    coordinated benefit. Because a single vehicle's contribution to the group is
    diluted by 1/(m * scale), the enforceable tau shrinks with m.

    We model enforceable tau_bar(m) = min(1, kappa / m) style ceiling derived below,
    then the platform uses tau = min(1, tau_bar(m)).
    """
    m = alpha / k if k > 0 else 0.0
    if m <= 1e-12:
        return 0.0  # infinitesimal group: cannot internalize (nonatomic)

    # Solve self-consistently: given tau, get equilibrium, compute temptation,
    # check obedience, bisect tau to the largest feasible value.
    def obedience_slack(tau):
        eq = lower_equilibrium(alpha, k, tau)
        xA = eq["xA"]
        # temptation for a vehicle told to go B: private gain of switching to A
        temptation = max(B_CONST - cA(xA), 0.0)
        # coordinated marginal benefit internalized per unit, that the vehicle
        # 'owes' the group by staying: tau * m * cA'(xA).
        # A single autonomous vehicle only feels a fraction of this returned to itself.
        # Model returned-to-self benefit as: tau * m * cA'(xA)  (the internalization
        # wedge the group is trying to sustain). Obedience holds if this wedge is at
        # least the temptation (staying is individually rational given the wedge).
        enforceable = tau * m * cA_prime(xA)
        # slack >= 0 means obedience satisfied
        return enforceable - temptation, eq

    # if even tau=1 satisfies obedience, return 1
    slack1, _ = obedience_slack(1.0)
    if slack1 >= 0:
        return 1.0
    # if tau=0 (no wedge, no temptation? temptation exists but enforceable=0)
    slack0, _ = obedience_slack(0.0)
    if slack0 < 0 and slack1 < 0:
        # even full internalization can't hold obedience once temptation>enforceable
        # find largest tau with slack>=0; if none, tau_bar = 0
        # bisection between 0 and 1
        lo, hi = 0.0, 1.0
        feasible = False
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            s, _ = obedience_slack(mid)
            if s >= 0:
                feasible = True
                hi = mid  # want largest feasible? slack increasing in tau, so search
                lo = mid
                break
            else:
                lo = mid
        # slack is monotone increasing in tau (enforceable ~ tau), so root is threshold
        # largest feasible tau is 1 if slack1>=0 else 0
        return 0.0 if slack1 < 0 else 1.0

    # general bisection: slack increasing in tau -> find min tau feasible, use up to 1
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        s, _ = obedience_slack(mid)
        if s >= 0:
            hi = mid
        else:
            lo = mid
    tau_threshold = hi
    # feasible region is [tau_threshold, 1]; platform picks the tau that minimizes J.
    return tau_threshold  # placeholder; refined in sweep


if __name__ == "__main__":
    alpha = 0.4
    print("=== references ===")
    print("all-selfish baseline:", cost_all_selfish_baseline(alpha))
    print("system optimal      :", cost_system_optimal())
    print()
    print("=== exogenous tau=1, sweep k ===")
    for k in [1, 2, 4, 8, 16, 64, 256]:
        eq = lower_equilibrium(alpha, k, tau=1.0)
        print(f"k={k:4d} m={eq['m']:.4f}  xA={eq['xA']:.4f} xB={eq['xB']:.4f}  J={eq['J']:.5f}")
