"""Large-sample regression checks for the corrected two-subgroup technical note.

Complements test_two_subgroup.py: that file pins exact numbers, this one runs
high-volume randomized sweeps (1e4 closed-form, 8e4 mergers) and re-verifies a
second irreducibility witness located in the *concave* branch of y(m)=beta(m)m.

Both witnesses are kept on purpose, because they point in opposite directions:
  concave branch  m=1.5 > 2/lambda=1 : equalizing helps, rho(1.5,1.5) > rho(1.8,1.2)
  convex branch   m=0.5 < 2/lambda=1 : equalizing hurts, rho(0.5,0.5) < rho(0.8,0.2)
Under any constant beta0 both pairs give identical rho, so neither ordering --
let alone their reversal -- is reproducible by scale-independent coordination.
"""

import numpy as np

from two_subgroup import solve_two_subgroup


def _solve(masses, lam, demand, beta_func=None):
    """Wrapper supplying alpha implied by the profile (Sum m_k = alpha N)."""
    masses = np.asarray(masses, dtype=float)
    return solve_two_subgroup(
        masses, lam=lam, alpha=float(masses.sum()) / demand, N=demand,
        beta_func=beta_func,
    )


def check_closed_form():
    """x1 from the bisection must equal the closed form (H + n xUE)/(1+n)."""
    rng = np.random.default_rng(11)
    max_error = 0.0
    for _ in range(10_000):
        k = int(rng.integers(1, 8))
        demand = float(rng.uniform(1.01, 5.0))
        av_mass = float(rng.uniform(0.01, demand))
        masses = rng.dirichlet(np.ones(k)) * av_mass
        lam = float(rng.uniform(0.1, 5.0))
        result = _solve(masses, lam, demand)
        if result["regime"] == "S":
            continue                      # closed form applies in regime (C)
        x_closed = (result["H"] + result["n_int"] * result["xUE"]) \
            / (1 + result["n_int"])
        max_error = max(max_error, abs(result["x1"] - x_closed))
    assert max_error < 1e-10, max_error
    return max_error


def check_mergers():
    """Merging two operators never lowers rho (Prop.: merging always helps)."""
    rng = np.random.default_rng(12)
    worst_change = 0.0
    for _ in range(80_000):
        k = int(rng.integers(2, 8))
        demand = float(rng.uniform(1.01, 6.0))
        av_mass = float(rng.uniform(0.01, demand))
        masses = list(rng.dirichlet(np.ones(k)) * av_mass)
        lam = float(rng.uniform(0.1, 6.0))
        i, j = rng.choice(k, 2, replace=False)
        merged = [m for idx, m in enumerate(masses) if idx not in (i, j)]
        merged.append(masses[i] + masses[j])
        before = _solve(masses, lam, demand)["rho"]
        after = _solve(merged, lam, demand)["rho"]
        worst_change = min(worst_change, after - before)
        assert after + 1e-10 >= before, (masses, lam, demand, before, after)
    return worst_change


def check_irreducibility_witness(structure_a, structure_b, demand, expect):
    """Endogenous beta separates A from B; every constant beta0 does not.

    expect: 'a>b' or 'a<b', the endogenous ordering this branch should show.
    """
    endogenous_a = _solve(structure_a, 2.0, demand)["rho"]
    endogenous_b = _solve(structure_b, 2.0, demand)["rho"]
    if expect == "a>b":
        assert endogenous_a > endogenous_b + 1e-6, (endogenous_a, endogenous_b)
    else:
        assert endogenous_a < endogenous_b - 1e-6, (endogenous_a, endogenous_b)

    max_exogenous_gap = 0.0
    for beta0 in np.linspace(0.001, 1.0, 1000):
        beta_func = lambda m, b=beta0: np.full_like(np.asarray(m, dtype=float), b)
        exogenous_a = _solve(structure_a, 2.0, demand, beta_func)["rho"]
        exogenous_b = _solve(structure_b, 2.0, demand, beta_func)["rho"]
        max_exogenous_gap = max(max_exogenous_gap, abs(exogenous_a - exogenous_b))
    assert max_exogenous_gap < 1e-10, max_exogenous_gap
    return endogenous_a, endogenous_b, max_exogenous_gap


if __name__ == "__main__":
    print("closed-form max error   :", check_closed_form())
    print("worst merger rho change :", check_mergers())
    # concave branch of y(m)=beta(m)m (each m=1.5 > 2/lambda=1): equalizing helps
    print("witness concave (m>2/l) :",
          check_irreducibility_witness([1.5, 1.5], [1.8, 1.2], 3.0, "a>b"))
    # convex branch (each m=0.5 < 2/lambda=1): equalizing hurts -- opposite sign
    print("witness convex  (m<2/l) :",
          check_irreducibility_witness([0.5, 0.5], [0.8, 0.2], 1.5, "a<b"))
