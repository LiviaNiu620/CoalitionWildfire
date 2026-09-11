"""
Two Sioux Falls robustness experiments that test claims the main text
currently makes only conditionally.

(A) Edge-specific fidelity.  The main model gives every operator a single
    scalar beta(m_k), justified by "coverage-balanced sampling" (A6).  The
    honest objection is that an operator's data are concentrated where its
    vehicles actually go, so fidelity should be edge-specific.  We implement

        beta_{k,e} = 1 - exp(-lambda * q_e * m_k),

    where q_e is the share of network travel that occurs on edge e under the
    free-flow assignment, normalised to mean one.  This is exogenous (it does
    not depend on the equilibrium f^k_e), so it does not create the nonlinear
    game that a beta_e(f^k_e) specification would; it only breaks the
    single-scalar reduction.  The question is whether the Sigma ordering
    failure reported in Table 2 gets worse, which is what the theory predicts:
    with edge-specific fidelity there is no scalar Sigma at all.

(B) Affine link costs.  Proposition 1(i) proves edge-flow uniqueness for
    affine costs but not for BPR.  All headline numbers use BPR, so the
    qualitative conclusions rest on a case where uniqueness is not certified.
    We rerun the whole table with c_e(x) = t0_e * (1 + 0.15 * x / kappa_e),
    i.e. the BPR form with exponent one, where Prop. 1(i) applies exactly.
    If the ordering survives, the conclusions do not depend on the
    uncertified case.
"""
from __future__ import annotations

import json
import numpy as np

import sf_equilibrium as SE
from sioux_falls_loader import build_sioux_falls

ALPHA = 0.9
LAM_TIMES_M = 6.78
CONFIGS = [
    ("8 equal", [1 / 8] * 8),
    ("4 equal", [1 / 4] * 4),
    ("2 equal", [1 / 2] * 2),
    ("(0.7,0.1,0.1,0.1)", [0.7, 0.1, 0.1, 0.1]),
    ("monopoly", [1.0]),
]


# ----------------------------------------------------------------------
# (A) edge-specific fidelity
# ----------------------------------------------------------------------
def edge_exposure(sf):
    """
    q_e = share of vehicle-distance on edge e under the free-flow assignment,
    normalised to mean one over used edges.  Exogenous by construction.
    """
    E = sf["edges"]
    x = np.zeros(E)
    for od in sf["od_list"]:
        x[od["paths"][0]] += od["demand"]
    q = np.zeros(E)
    used = x > 0
    q[used] = x[used] / x[used].mean()
    q[~used] = q[used].min() if used.any() else 1.0
    return q


class EdgeBetaGame(SE.SFGame):
    """SFGame with beta_{k,e} instead of beta_k."""

    def set_edge_beta(self, q, lam_beta):
        self.q = q
        q = np.asarray(q)
        exposure = np.outer(self.masses, q) if q.ndim == 1 else self.masses[:, None] * q
        self.betas_e = 1.0 - np.exp(-lam_beta * exposure)  # (K,E)

    def operator(self, fH, fF):
        x = self.edge_flow(fH, fF)
        ce = SE.bpr_cost(x, self.t0, self.cap)
        cp = SE.bpr_deriv(x, self.t0, self.cap)
        pc = self.A.T @ ce
        gH = pc.copy()
        xk = fF @ self.A.T                       # (K,E)
        gF = pc[None, :] + ((self.betas_e * xk * cp[None, :]) @ self.A)
        return gH, gF, x, ce, cp

    def edge_weights(self, fH, fF):
        x = self.edge_flow(fH, fF)
        ce = SE.bpr_cost(x, self.t0, self.cap)
        cp = SE.bpr_deriv(x, self.t0, self.cap)
        xk = fF @ self.A.T
        return ce, ce[None, :] + self.betas_e * xk * cp[None, :]


# ----------------------------------------------------------------------
# (B) affine costs
# ----------------------------------------------------------------------
def patch_affine(on):
    """Switch the module-level BPR exponent between 4 (default) and 1."""
    SE.BPR_ETA = 4.0 if not on else 1.0


# ----------------------------------------------------------------------
def solve_variant(sf, shares, lam_beta, mode, edge_beta=None, tol=1e-7):
    """Column generation with an optional edge-specific-beta game class."""
    import os
    from sioux_falls_loader import parse_net, DATA_DIR
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    G, edge_idx = SE.build_graph(links)
    routes = SE.RouteSet(sf["od_list"])
    cls = EdgeBetaGame if edge_beta is not None else SE.SFGame
    game = cls(sf, ALPHA, shares, lam_beta, routes, mode=mode)
    if edge_beta is not None:
        game.set_edge_beta(edge_beta, lam_beta)
    warm = None
    for rnd in range(25):
        sol = game.solve(seed=0, tol=tol, warm=warm)
        wH, wF = game.edge_weights(sol["fH"], sol["fF"])
        scale = float(np.max(game.A.T @ wH))
        viol, worst = SE.certificate(game, sol, G, edge_idx, 1e-6 * scale)
        if not viol:
            sol.update(certified=True, n_paths=routes.total(),
                       max_reduced_cost=worst, cg_rounds=rnd + 1)
            return sol, game
        old = routes.keys()
        if sum(routes.add(oi, ep) for oi, ep, _ in viol) == 0:
            sol.update(certified=False, n_paths=routes.total(),
                       max_reduced_cost=worst, cg_rounds=rnd + 1)
            return sol, game
        pos = {k: i for i, k in enumerate(old)}
        game.rebuild()
        if edge_beta is not None:
            game.set_edge_beta(edge_beta, lam_beta)
        fH = np.zeros(game.P)
        fF = np.zeros((game.K, game.P))
        for j, k in enumerate(routes.keys()):
            i = pos.get(k)
            if i is not None:
                fH[j] = sol["fH"][i]
                fF[:, j] = sol["fF"][:, i]
        warm = (fH, fF)
    sol.update(certified=False, n_paths=routes.total(),
               max_reduced_cost=worst, cg_rounds=25)
    return sol, game


def operator_specific_exposure(sf, K, specialization):
    """Deterministic geographic specialization around the common exposure q_e."""
    import os
    from sioux_falls_loader import parse_net, DATA_DIR
    common = edge_exposure(sf)
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    origins = np.array([u for u, _, _, _ in links])
    zones = (origins - 1) % K
    q = np.empty((K, sf["edges"]))
    for k in range(K):
        raw = common * (1.0 + specialization * (zones == k))
        q[k] = raw * common.sum() / raw.sum()
    return q


def run_specialization_sweep():
    """Test whether the two-operator interior optimum survives geographic specialization."""
    patch_affine(False)
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(o["demand"] for o in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)
    ue, _ = solve_variant(sf, [1.0], lam_beta, "ue", None)
    so, _ = solve_variant(sf, [1.0], lam_beta, "so", None)
    gap = ue["J"] - so["J"]
    structures = [("8 equal", [1 / 8] * 8), ("2 equal", [0.5, 0.5]), ("monopoly", [1.0])]
    levels = [0.0, 0.25, 0.50, 0.75]
    rows = {}
    print("\nOperator-specific geographic exposure")
    print(f"{'specialization':>15}{'8 equal':>12}{'2 equal':>12}{'monopoly':>12}")
    for level in levels:
        vals = {}
        for name, shares in structures:
            q = operator_specific_exposure(sf, len(shares), level)
            sol, _ = solve_variant(sf, shares, lam_beta, "eq", q)
            vals[name] = (ue["J"] - sol["J"]) / gap
        rows[str(level)] = vals
        print(f"{level:>15.2f}{100*vals['8 equal']:>11.1f}%"
              f"{100*vals['2 equal']:>11.1f}%{100*vals['monopoly']:>11.1f}%")
    return rows


def run_experiment(name, edge_beta_on, affine_on):
    patch_affine(affine_on)
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(o["demand"] for o in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)
    q = edge_exposure(sf) if edge_beta_on else None

    print("\n" + "=" * 78)
    print(f"{name}   (edge-specific beta = {edge_beta_on}, affine costs = {affine_on})")
    if edge_beta_on:
        print(f"  exposure q_e: min={q.min():.3f} median={np.median(q):.3f} "
              f"max={q.max():.3f}")
    print("=" * 78)

    ue, _ = solve_variant(sf, [1.0], lam_beta, "ue", None)
    so, _ = solve_variant(sf, [1.0], lam_beta, "so", None)
    gap = ue["J"] - so["J"]
    print(f"  J_UE={ue['J']:.0f} (cert={ue['certified']})  "
          f"J_SO={so['J']:.0f} (cert={so['certified']})  gap={gap:.0f} "
          f"({100*gap/ue['J']:.2f}%)")

    rows = []
    print(f"  {'config':<20}{'Sigma':>9}{'HHI':>8}{'J':>13}{'rho':>9}"
          f"{'gap':>10}{'cert':>7}")
    for nm, sh in CONFIGS:
        sol, game = solve_variant(sf, sh, lam_beta, "eq", q)
        if edge_beta_on:
            # descriptive analogue: average over edges, weighted by exposure
            b_bar = (game.betas_e * q[None, :]).sum(axis=1) / q.sum()
            Sigma = float(np.sum(1.0 / b_bar))
        else:
            Sigma = float(np.sum(1.0 / game.betas))
        hhi = float(np.sum(np.array(sh) ** 2))
        rho = (ue["J"] - sol["J"]) / gap
        rows.append(dict(name=nm, Sigma=Sigma, HHI=hhi, J=sol["J"], rho=rho,
                         gap=sol["gap"], certified=bool(sol["certified"]),
                         n_paths=sol["n_paths"]))
        print(f"  {nm:<20}{Sigma:>9.2f}{hhi:>8.3f}{sol['J']:>13.0f}"
              f"{100*rho:>8.1f}%{sol['gap']:>10.1e}{str(sol['certified']):>7}")
    patch_affine(False)
    return dict(name=name, J_UE=ue["J"], J_SO=so["J"], gap=gap, rows=rows)


if __name__ == "__main__":
    out = {}
    out["edge_beta"] = run_experiment("(A) edge-specific fidelity, BPR eta=4",
                                      True, False)
    out["affine"] = run_experiment("(B) affine link costs (Prop 1(i) applies)",
                                   False, True)
    out["operator_specific_exposure"] = run_specialization_sweep()
    with open("sf_robustness_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nwritten: sf_robustness_results.json")
