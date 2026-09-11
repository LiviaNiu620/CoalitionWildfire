"""
Emit the two Sioux Falls LaTeX tables from sf_certified_results.json.

Table 1 (tab:sf)      : the substantive table -- Sigma, HHI, J, rho.
Table 2 (tab:sf-diag) : numerical certificates -- VI residual, reduced-cost
                        slack, paths after column generation, multi-start
                        dispersion, minimum eigenvalue diagnostic.
"""
import json

with open("sf_certified_results.json") as f:
    R = json.load(f)

rows = R["rows"]
JUE, JSO, gap = R["J_UE"], R["J_SO"], R["gap"]


def fnum(x, d=0):
    return f"{x:,.{d}f}".replace(",", "{,}")


print(r"% ---------- Table: substantive ----------")
print(r"\begin{table}[t]\centering")
print(r"\begin{tabular}{lccrr}")
print(r"\toprule")
print(r"Size distribution & $\Sigma=\sum_k 1/\beta(m_k)$ & $\HHI$ & $J$ & $\rho$ \\")
print(r"\midrule")
for r in rows:
    print(f"{r['name']} & {r['Sigma']:.1f} & {r['HHI']:.3f} & "
          f"${fnum(r['J'])}$ & ${100*r['rho']:.1f}\\%$ \\\\")
print(r"\bottomrule")
print(r"\end{tabular}")
print(r"\caption{Sioux Falls ($\alpha=0.9$, $\lambda_\beta\alpha N=6.78$). "
      rf"Benchmarks $J^{{UE}}={fnum(JUE)}$ and $J^{{SO}}={fnum(JSO)}$ "
      rf"give a recoverable gap of ${fnum(gap)}$ "
      rf"(${100*gap/JUE:.2f}\%$ of $J^{{UE}}$); all three are computed on "
      r"column-generated route sets with the certificates of "
      r"Table~\ref{tab:sf-diag}.}")
print(r"\label{tab:sf}")
print(r"\end{table}")

print()
print(r"% ---------- Table: certificates ----------")
print(r"\begin{table}[t]\centering\small")
print(r"\begin{tabular}{lccccc}")
print(r"\toprule")
print(r"Size distribution & VI residual & reduced-cost slack & paths "
      r"& multi-start $\max|\Delta x|/\max x$ & $\lambda_{\min}$ \\")
print(r"\midrule")
for r in rows:
    print(f"{r['name']} & ${r['gap']:.1e}$ & ${r['reduced']:+.1e}$ & "
          f"{r['n_paths']} & ${r['multistart_dx']:.1e}$ & "
          f"${r['min_eig']:+.3f}$ \\\\")
print(r"\bottomrule")
print(r"\end{tabular}")
print(r"\caption{Numerical certificates. The VI residual is \eqref{eq:vigap}, "
      r"the normalized complementarity gap over \emph{all} players. The "
      r"reduced-cost slack is the largest amount by which an omitted route "
      r"undercuts a player's held routes, computed exactly by one Dijkstra "
      r"run per (player, origin) under that player's own edge weights; a "
      r"nonpositive value certifies that no omitted route is profitable. "
      r"Multi-start dispersion is the largest relative edge-flow deviation "
      r"across three Dirichlet-random initializations. "
      r"$\lambda_{\min}$ is the minimum eigenvalue of the normalized "
      r"symmetric edge Jacobian at the computed profile; by (C2)--(C3) of "
      r"Appendix~\ref{app:wellposed} it certifies nothing globally and is "
      r"reported as a diagnostic only.}")
print(r"\label{tab:sf-diag}")
print(r"\end{table}")

print()
print("% ---------- plain summary ----------")
for r in rows:
    print(f"% {r['name']:<20} Sigma={r['Sigma']:7.2f} HHI={r['HHI']:.3f} "
          f"J={r['J']:12.1f} rho={100*r['rho']:6.2f}%  "
          f"cert={r['certified']} cg={r['cg_rounds']} "
          f"beta=[{r['beta_min']:.3f},{r['beta_max']:.3f}]")
