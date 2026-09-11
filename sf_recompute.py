"""
Recompute the Sioux Falls table with certified numerics.

Produces, for each ownership structure:
    Sigma = sum_k 1/beta(m_k),  HHI,  J,  rho,  VI gap,  reduced-cost slack,
    multi-start edge-flow dispersion,  min eigenvalue of the symmetric Jacobian.

Benchmarks J^UE and J^SO are computed with the same column-generation loop, so
all three quantities in rho refer to the same (certified) route set.

Calibration note
----------------
The paper reports lambda_beta * alpha * N ~ 6.78, which reproduces the Sigma
column (8 equal -> 14.0, 4 equal -> 4.9, 2 equal -> 2.1, skewed -> 7.1,
monopoly -> 1.0).  The committed sf_report.py used lambda_beta = 3e-4, i.e.
lambda_beta * alpha * N = 97.4, under which every operator has beta ~ 1 and
Sigma degenerates to K.  We use the value consistent with the reported Sigma.
"""
import json
import numpy as np

from sioux_falls_loader import build_sioux_falls
from sf_equilibrium import (SFGame, RouteSet, solve_with_column_generation,
                            bpr_cost, bpr_deriv, build_graph, certificate)
from sioux_falls_loader import parse_net, DATA_DIR
import os

ALPHA = 0.9
LAM_TIMES_M = 6.78           # lambda_beta * alpha * N, as reported in the paper
CONFIGS = [
    ("8 equal",          [1 / 8] * 8),
    ("4 equal",          [1 / 4] * 4),
    ("2 equal",          [1 / 2] * 2),
    ("(0.7,0.1,0.1,0.1)", [0.7, 0.1, 0.1, 0.1]),
    ("monopoly",         [1.0]),
]
N_SEEDS = 3
MAIN_TOL = 1e-7      # relative complementarity residual for the reported solution
MULTI_TOL = 1e-6     # looser: multi-start only needs to show the same edge flow


def min_sym_eigenvalue(game, sol):
    """Minimum eigenvalue of the normalised symmetric edge Jacobian, over
    positive-flow edges.  Reported as a diagnostic only: it is evaluated at the
    computed profile and therefore certifies nothing globally (see Prop. 1)."""
    x = sol["x"]
    pos = x > 1e-9
    cp = bpr_deriv(x, game.t0, game.cap)
    worst = np.inf
    K = game.K
    for e in np.where(pos)[0]:
        tau = np.array([game.betas[k] * (game.A[e] @ game.fF_k(sol, k)) * 3.0 / x[e]
                        for k in range(K)])
        one = np.ones(K + 1)
        M = np.outer(one, one)
        M[:K, :K] += np.diag(game.betas)
        for k in range(K):
            ek = np.zeros(K + 1)
            ek[k] = 1.0
            M += 0.5 * tau[k] * (np.outer(ek, one) + np.outer(one, ek))
        worst = min(worst, float(np.linalg.eigvalsh(M)[0]))
    return worst


def _fF_k(self, sol, k):
    return sol["fF"][k]


SFGame.fF_k = _fF_k


def main():
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total_dem = sum(o["demand"] for o in sf["od_list"])
    lam_beta = LAM_TIMES_M / (ALPHA * total_dem)

    print("=" * 78)
    print("Sioux Falls, certified recomputation")
    print(f"  links={sf['edges']}  OD pairs={len(sf['od_list'])}  "
          f"total demand={total_dem:.0f}")
    print(f"  alpha={ALPHA}  AV mass={ALPHA*total_dem:.0f}  "
          f"lambda_beta={lam_beta:.6e}  lambda_beta*AV mass={LAM_TIMES_M}")
    print(f"  BPR: b={0.15}, eta={4}")
    print("=" * 78)

    print("\n--- benchmarks ---")
    cache_path = "sf_benchmarks.json"
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            bm = json.load(f)
        J_UE, J_SO = bm["J_UE"], bm["J_SO"]
        ue_paths, so_paths = bm["ue_paths"], bm["so_paths"]
        print(f"  (cached)  J^UE = {J_UE:.1f}   J^SO = {J_SO:.1f}")
    else:
        sol_ue, _ = solve_with_column_generation(sf, ALPHA, [1.0], lam_beta,
                                                 mode="ue", verbose=True)
        J_UE, ue_paths = sol_ue["J"], sol_ue["n_paths"]
        print(f"  J^UE = {J_UE:.1f}   certified={sol_ue['certified']}  "
              f"paths={ue_paths}  gap={sol_ue['gap']:.2e}")
        sol_so, _ = solve_with_column_generation(sf, ALPHA, [1.0], lam_beta,
                                                 mode="so", verbose=True)
        J_SO, so_paths = sol_so["J"], sol_so["n_paths"]
        print(f"  J^SO = {J_SO:.1f}   certified={sol_so['certified']}  "
              f"paths={so_paths}  gap={sol_so['gap']:.2e}")
        with open(cache_path, "w") as f:
            json.dump(dict(J_UE=J_UE, J_SO=J_SO, ue_paths=ue_paths,
                           so_paths=so_paths), f, indent=2)
    gap_abs = J_UE - J_SO
    print(f"  recoverable gap = {gap_abs:.1f} "
          f"({100*gap_abs/J_UE:.2f}% of J^UE)", flush=True)

    rows = []
    cfg_cache = "sf_config_cache.json"
    done = json.load(open(cfg_cache)) if os.path.exists(cfg_cache) else {}
    for name, shares in CONFIGS:
        if name in done:
            rows.append(done[name])
            print(f"\n--- {name} (cached) rho={100*done[name]['rho']:.1f}% ---")
            continue
        print(f"\n--- {name} ---")
        sol, game = solve_with_column_generation(sf, ALPHA, shares, lam_beta,
                                                 mode="eq", verbose=True,
                                                 tol=MAIN_TOL)
        betas = game.betas
        Sigma = float(np.sum(1.0 / betas))
        hhi = float(np.sum(np.array(shares) ** 2))
        rho = (J_UE - sol["J"]) / gap_abs

        # multi-start on THIS configuration (not only the symmetric one)
        devs, dJs = [], []
        for s in range(1, N_SEEDS + 1):
            alt = game.solve(seed=s, tol=MULTI_TOL)
            devs.append(float(np.max(np.abs(alt["x"] - sol["x"]))
                              / max(np.max(sol["x"]), 1e-9)))
            dJs.append(abs(alt["J"] - sol["J"]) / sol["J"])
        mineig = min_sym_eigenvalue(game, sol)

        rows.append(dict(name=name, shares=shares, Sigma=Sigma, HHI=hhi,
                         J=sol["J"], rho=rho, gap=sol["gap"],
                         reduced=sol["max_reduced_cost"],
                         certified=bool(sol["certified"]),
                         n_paths=sol["n_paths"], cg_rounds=sol["cg_rounds"],
                         beta_min=float(betas.min()), beta_max=float(betas.max()),
                         multistart_dx=max(devs), multistart_dJ=max(dJs),
                         min_eig=mineig))
        done[name] = rows[-1]
        with open(cfg_cache, "w") as f:
            json.dump(done, f, indent=2)
        print(f"  Sigma={Sigma:.3f} HHI={hhi:.3f} J={sol['J']:.1f} "
              f"rho={100*rho:.2f}%  gap={sol['gap']:.2e} "
              f"reduced={sol['max_reduced_cost']:+.2e} "
              f"multistart dx={max(devs):.2e} dJ={max(dJs):.2e} "
              f"min_eig={mineig:+.4f}")

    print("\n" + "=" * 78)
    print(f"{'config':<20}{'Sigma':>8}{'HHI':>8}{'J':>14}{'rho':>9}"
          f"{'VI gap':>10}{'red.cost':>11}{'multi dx':>10}{'min eig':>9}")
    for r in rows:
        print(f"{r['name']:<20}{r['Sigma']:>8.2f}{r['HHI']:>8.3f}{r['J']:>14.0f}"
              f"{100*r['rho']:>8.1f}%{r['gap']:>10.1e}{r['reduced']:>+11.1e}"
              f"{r['multistart_dx']:>10.1e}{r['min_eig']:>+9.3f}")
    print("=" * 78)

    out = dict(alpha=ALPHA, total_demand=total_dem, lam_beta=lam_beta,
               lam_times_M=LAM_TIMES_M, J_UE=J_UE, J_SO=J_SO,
               gap=gap_abs, ue_paths=ue_paths,
               so_paths=so_paths, rows=rows)
    with open("sf_certified_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("written: sf_certified_results.json")


if __name__ == "__main__":
    main()
