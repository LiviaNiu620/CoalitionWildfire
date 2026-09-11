"""
Sweep the prior-slope multiplier kappa on Sioux Falls.

Why this exists
---------------
The first version of the numerical section ran every experiment with
theta_bar = 0, i.e. beta0 = 0.  Under the Bayesian reading of Sec. 3.4 that is
not an admissible parameter value: beta0 = theta_bar_e / theta_e and the model
requires theta_bar_e > 0 (a zero-mean prior on a strictly positive slope is
exactly the inconsistency the reviewer identified).  Running everything at
beta0 = 0 therefore evaluated the model outside its own assumptions, and in
particular left the over-internalisation regime (R+) untested.

Parametrisation
---------------
The prior mean slope on edge e is

    theta_bar_e = kappa * c'_e(x^UE_e),

so at user-equilibrium flows beta0_e = kappa exactly.  kappa < 1 means the
prior understates congestibility (regime R-), kappa > 1 means it overstates it
(regime R+).  The perceived slope is

    theta_hat^k_e = theta_bar_e + beta(m_k) * (c'_e - theta_bar_e),

which is the shrinkage formula of eq. (posterior), and beta^eff_{k,e} =
theta_hat^k_e / c'_e.
"""
import json
import numpy as np

from sioux_falls_loader import build_sioux_falls, parse_net, DATA_DIR
import sf_equilibrium as SE
import os

ALPHA = 0.9
LAM_TIMES_M = 6.78
KAPPAS = [0.0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50]
CONFIGS = [
    ("8 equal", [1 / 8] * 8),
    ("4 equal", [1 / 4] * 4),
    ("2 equal", [1 / 2] * 2),
    ("(0.7,0.1,0.1,0.1)", [0.7, 0.1, 0.1, 0.1]),
    ("monopoly", [1.0]),
]


def solve(sf, shares, lam_beta, mode, theta_bar=None, tol=1e-7):
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    G, edge_idx = SE.build_graph(links)
    routes = SE.RouteSet(sf["od_list"])
    game = SE.SFGame(sf, ALPHA, shares, lam_beta, routes, mode=mode,
                     theta_bar=theta_bar)
    warm = None
    for rnd in range(25):
        sol = game.solve(seed=0, tol=tol, warm=warm)
        wH, _ = game.edge_weights(sol["fH"], sol["fF"])
        scale = float(np.max(game.A.T @ wH))
        viol, worst = SE.certificate(game, sol, G, edge_idx, 1e-6 * scale)
        if not viol:
            sol.update(certified=True, n_paths=routes.total(),
                       max_reduced_cost=worst)
            return sol, game
        old = routes.keys()
        if sum(routes.add(oi, ep) for oi, ep, _ in viol) == 0:
            break
        pos = {k: i for i, k in enumerate(old)}
        game.rebuild()
        fH = np.zeros(game.P); fF = np.zeros((game.K, game.P))
        for j, k in enumerate(routes.keys()):
            i = pos.get(k)
            if i is not None:
                fH[j] = sol["fH"][i]; fF[:, j] = sol["fF"][:, i]
        warm = (fH, fF)
    sol.update(certified=False, n_paths=routes.total(), max_reduced_cost=worst)
    return sol, game


def main():
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(o["demand"] for o in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total)

    ue, gue = solve(sf, [1.0], lam_beta, "ue")
    so, _ = solve(sf, [1.0], lam_beta, "so")
    gap = ue["J"] - so["J"]
    cp_ue = SE.bpr_deriv(ue["x"], gue.t0, gue.cap)     # c'_e at user equilibrium
    print(f"J_UE={ue['J']:.0f}  J_SO={so['J']:.0f}  gap={gap:.0f} "
          f"({100*gap/ue['J']:.2f}%)")
    print(f"c'_e(x^UE): min={cp_ue[cp_ue>0].min():.3e} "
          f"median={np.median(cp_ue):.3e} max={cp_ue.max():.3e}")

    out = {"J_UE": ue["J"], "J_SO": so["J"], "gap": gap, "sweeps": {}}
    for kappa in KAPPAS:
        tb = None if kappa == 0.0 else kappa * cp_ue
        regime = "beta0=0 (inadmissible)" if kappa == 0.0 else \
                 ("R- under-internalisation" if kappa < 1 else
                  "R+ OVER-internalisation")
        print(f"\n--- kappa = {kappa}   [{regime}] ---")
        print(f"  {'config':<20}{'J':>13}{'rho':>9}{'beta_eff rng':>20}"
              f"{'gap':>10}{'cert':>7}")
        rows = []
        for nm, sh in CONFIGS:
            sol, game = solve(sf, sh, lam_beta, "eq", theta_bar=tb)
            cp = SE.bpr_deriv(sol["x"], game.t0, game.cap)
            xk = sol["fF"] @ game.A.T
            th = game.perceived_slope(cp, xk)
            pos = cp > 1e-12
            beff = th[:, pos] / cp[None, pos]
            rho = (ue["J"] - sol["J"]) / gap
            rows.append(dict(name=nm, J=sol["J"], rho=rho,
                             beff_min=float(beff.min()),
                             beff_max=float(beff.max()),
                             gap=sol["gap"], certified=bool(sol["certified"])))
            print(f"  {nm:<20}{sol['J']:>13.0f}{100*rho:>8.1f}%"
                  f"   [{beff.min():.2f},{beff.max():.2f}]".ljust(20)
                  + f"{sol['gap']:>10.1e}{str(sol['certified']):>7}")
        out["sweeps"][str(kappa)] = rows

    print("\n" + "=" * 78)
    print("rho by structure and prior multiplier kappa")
    print(f"{'kappa':>8}" + "".join(f"{n.split()[0]:>11}" for n, _ in CONFIGS))
    for kappa in KAPPAS:
        r = out["sweeps"][str(kappa)]
        print(f"{kappa:>8.2f}" + "".join(f"{100*x['rho']:>10.1f}%" for x in r))
    print("=" * 78)
    with open("sf_kappa_sweep.json", "w") as f:
        json.dump(out, f, indent=2)
    print("written: sf_kappa_sweep.json")


if __name__ == "__main__":
    main()
