#!/usr/bin/env python3
"""Recompute and visualize the Sioux Falls two-fleet/monopoly displacement.

The figure uses the certified edge-specific-fidelity model.  Positive edge cost
differences mean that the two-fleet equilibrium lowers the edge's TSTT
contribution relative to the monopoly; negative values are offsetting losses.
The accompanying scatter relates those cost changes to the change in the
aggregate perceived own-flow response term.
"""

from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import sf_equilibrium as SE  # noqa: E402
import sf_robustness as SR  # noqa: E402
from sioux_falls_loader import build_sioux_falls, parse_net, DATA_DIR  # noqa: E402


def solve(shares, sf, q, lam_beta):
    sol, game = SR.solve_variant(sf, shares, lam_beta, "eq", q, tol=1e-7)
    if not sol.get("certified", False) or sol["gap"] >= 1.1e-7:
        raise RuntimeError(f"uncertified solution for {shares}: {sol}")
    x = sol["x"]
    edge_cost = x * SE.bpr_cost(x, game.t0, game.cap)
    xk = sol["fF"] @ game.A.T
    response = np.sum(game.betas_e * xk, axis=0)
    return sol, game, edge_cost, response


def main():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.2,
        "axes.titlesize": 10.2,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.2,
        "legend.fontsize": 8.1,
        "legend.frameon": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.16,
    })

    SR.patch_affine(False)
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(o["demand"] for o in sf["od_list"])
    lam_beta = SR.LAM_TIMES_M / (SR.ALPHA * total)
    q = SR.edge_exposure(sf)

    two, game2, cost2, response2 = solve([0.5, 0.5], sf, q, lam_beta)
    mono, gamem, costm, responsem = solve([1.0], sf, q, lam_beta)

    savings = costm - cost2
    response_change = responsem - response2
    positive = float(savings[savings > 0].sum())
    negative = float(-savings[savings < 0].sum())
    net = float(savings.sum())
    corr = float(np.corrcoef(savings, response_change)[0, 1])

    links = parse_net(str(Path(DATA_DIR) / "SiouxFalls_net.tntp"))
    labels = np.array([f"{u}→{v}" for u, v, _, _ in links])
    order = np.argsort(savings)
    colors = np.where(savings[order] >= 0, "#0072B2", "#D55E00")

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.05),
                             gridspec_kw={"width_ratios": [1.18, 1.0]})
    ax = axes[0]
    ax.bar(np.arange(len(order)), savings[order], color=colors, width=0.86,
           edgecolor="none")
    ax.axhline(0, color="#4B5563", linewidth=0.8)
    ax.set_xlabel("directed edges, sorted by contribution change")
    ax.set_ylabel("monopoly minus two-fleet edge TSTT")
    ax.set_title("Gains on some edges are partly offset elsewhere")
    ax.set_xticks([])
    ax.text(0.03, 0.96,
            f"gross gains = {positive:,.0f}\n"
            f"offsetting losses = {negative:,.0f}\n"
            f"net = {net:,.0f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.1,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                  "edgecolor": "#D1D5DB", "alpha": 0.92})

    ax = axes[1]
    scatter = ax.scatter(response_change, savings, c=q, cmap="viridis", s=25,
                         alpha=0.88, edgecolor="white", linewidth=0.35)
    coef = np.polyfit(response_change, savings, 1)
    xx = np.linspace(response_change.min(), response_change.max(), 100)
    ax.plot(xx, np.polyval(coef, xx), color="#4B5563", linestyle="--",
            linewidth=1.0)
    ax.axhline(0, color="#9CA3AF", linewidth=0.7)
    ax.axvline(0, color="#9CA3AF", linewidth=0.7)
    ax.set_xlabel(r"monopoly minus two-fleet response term")
    ax.set_ylabel("monopoly minus two-fleet edge TSTT")
    ax.set_title("Spatial response and cost displacement")
    ax.text(0.05, 0.95, f"edge correlation = {corr:.2f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.2)
    colorbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.03)
    colorbar.set_label(r"free-flow exposure $q_e$", fontsize=8.2)
    colorbar.ax.tick_params(labelsize=7.2)

    fig.tight_layout(w_pad=2.0)
    fig.savefig(OUT / "fig_sf_displacement.pdf")
    fig.savefig(OUT / "fig_sf_displacement.png", dpi=300)
    plt.close(fig)

    payload = {
        "two_fleet_J": float(two["J"]),
        "monopoly_J": float(mono["J"]),
        "gross_savings": positive,
        "offsetting_losses": negative,
        "net_savings": net,
        "edge_correlation": corr,
        "two_fleet_gap": float(two["gap"]),
        "monopoly_gap": float(mono["gap"]),
        "two_fleet_paths": int(two["n_paths"]),
        "monopoly_paths": int(mono["n_paths"]),
    }
    with open(OUT / "fig_sf_displacement_data.json", "w") as handle:
        json.dump(payload, handle, indent=2)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
