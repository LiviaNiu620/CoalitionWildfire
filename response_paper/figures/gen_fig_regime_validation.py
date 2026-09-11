"""Regime theorem, cross-network validation, and policy regret figure."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
COLORS = ["#3B6FB6", "#0B8A83", "#E07A3F", "#7A5AA6"]

def main():
    cross = json.loads((ROOT / "cross_network_regime_summary.json").read_text())
    design = json.loads((ROOT / "sf_coalition_design_analysis.json").read_text())
    fig = plt.figure(figsize=(10.5, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 0.95])
    ax = fig.add_subplot(gs[0, :2]); sigma = np.linspace(0, 3, 400); xue = 1.0
    for ratio, color, label in [(0.2, COLORS[0], r"$H/x^{UE}=0.2$"), (0.6, COLORS[2], r"$H/x^{UE}=0.6$")]:
        h = ratio * xue; x = (h + sigma * xue) / (1 + sigma); j = x ** 2 - xue * x; j -= j.min()
        ax.plot(sigma, j, lw=2.2, color=color, label=label)
        star = 1 - 2 * ratio
        if star > 0:
            ax.axvline(star, color=color, ls="--", lw=1); ax.text(star + .03, .01, r"$\Sigma^*$", color=color, fontsize=9)
    ax.set_xlabel(r"effective coalition response $\Sigma$"); ax.set_ylabel("TSTT above within-curve minimum")
    ax.set_title("Canonical regime threshold"); ax.grid(alpha=.22); ax.legend(frameon=False)
    alphas = [0.3, 0.5, 0.7, 0.9]; gammas = [0.25, 0.5, 0.75, 1.0]
    for col, network in enumerate(("braess", "grid")):
        axh = fig.add_subplot(gs[0, 2] if col == 0 else gs[1, 2]); mat = np.zeros((4, 4))
        for row in cross["states"]:
            if row["network"] == network: mat[alphas.index(row["alpha"]), gammas.index(row["gamma"])] = row["best_K"]
        axh.imshow(mat, vmin=1, vmax=4, cmap=ListedColormap(COLORS), aspect="auto", origin="lower")
        axh.set_xticks(range(4), gammas); axh.set_yticks(range(4), alphas); axh.set_xlabel(r"$\gamma$"); axh.set_ylabel(r"$\alpha$")
        axh.set_title(f"Best K: {network.title()}")
        for i in range(4):
            for j in range(4): axh.text(j, i, str(int(mat[i,j])), ha="center", va="center", color="white", weight="bold")
    axr = fig.add_subplot(gs[1, :2]); labels = ["1-step", "2-step", "Always grand", "Always singleton"]
    keys = ["greedy_merge", "two_step_merge", "always_grand", "always_singleton"]
    means = [design["policy_regret_summary"][k]["mean_recovery_regret_pp"] for k in keys]
    maxima = [design["policy_regret_summary"][k]["max_recovery_regret_pp"] for k in keys]
    x = np.arange(4); w = .36
    axr.bar(x-w/2, means, w, color=COLORS[1], label="Mean regret"); axr.bar(x+w/2, maxima, w, color=COLORS[2], label="Max regret")
    axr.set_xticks(x, labels); axr.set_ylabel("Recovery regret (pp)"); axr.set_title("Sioux Falls policy regret")
    axr.grid(axis="y", alpha=.22); axr.legend(frameon=False)
    fig.savefig(OUT / "fig_regime_validation.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_regime_validation.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "fig_regime_validation.svg", bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__": main()
