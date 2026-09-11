"""Generate the certified Sioux Falls coalition-composition figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "sf_coalition_composition_results.json"
SUMMARY = ROOT / "sf_coalition_composition_summary.json"
OUT = Path(__file__).resolve().parent

COLORS = {
    "singletons": "#7A7A7A",
    "assortative": "#0072B2",
    "mixed": "#D55E00",
    "grand": "#009E73",
}
MARKERS = {
    "singletons": "o",
    "assortative": "s",
    "mixed": "D",
    "grand": "^",
}
LABELS = {
    "singletons": "four singletons",
    "assortative": "assortative pairs",
    "mixed": "mixed-type pairs",
    "grand": "grand coalition",
}


def setup_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.2,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.0,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "grid.linestyle": "-",
        "lines.linewidth": 1.8,
        "lines.markersize": 5.2,
        "savefig.dpi": 450,
        "savefig.bbox": "tight",
    })


def main():
    setup_style()
    raw = json.loads(RAW.read_text())
    summary = json.loads(SUMMARY.read_text())
    rows = raw["profiles"]
    lookup = {
        (row["alpha"], row["gamma"], row["partition"]): row
        for row in rows
    }
    gammas = np.array(raw["protocol"]["gammas"], dtype=float)

    fig = plt.figure(figsize=(7.0, 5.25))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.78], hspace=0.48, wspace=0.28)
    top_axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    bottom = fig.add_subplot(grid[1, :])

    for ax, alpha in zip(top_axes, (0.5, 0.9)):
        for name in ("singletons", "assortative", "mixed", "grand"):
            values = [100.0 * lookup[(alpha, float(g), name)]["recovery"] for g in gammas]
            linestyle = "--" if name == "singletons" else "-"
            ax.plot(
                gammas, values, label=LABELS[name], color=COLORS[name],
                marker=MARKERS[name], linestyle=linestyle,
            )
        ax.set_title(rf"AV penetration $\alpha={alpha:.1f}$")
        ax.set_xlabel(r"cross-member coordination $\gamma$")
        ax.set_ylabel("recoverable congestion gap closed (%)")
        ax.set_xticks(gammas)
        ax.margins(x=0.03)

    handles, labels = top_axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.01),
        handlelength=2.3, columnspacing=1.2,
    )

    contrast = {
        (row["alpha"], row["gamma"]): row
        for row in summary["assortative_minus_mixed"]
    }
    for alpha, color, marker in ((0.5, "#56B4E9", "o"), (0.9, "#CC79A7", "s")):
        values = [contrast[(alpha, float(g))]["delta_recovery_percentage_points"]
                  for g in gammas]
        bottom.plot(
            gammas, values, color=color, marker=marker,
            label=rf"$\alpha={alpha:.1f}$",
        )
    bottom.axhline(0.0, color="#333333", linewidth=0.9)
    bottom.set_xlabel(r"cross-member coordination $\gamma$")
    bottom.set_ylabel("assortative minus mixed recovery (pp)")
    bottom.set_xticks(gammas)
    bottom.legend(loc="lower left", ncol=2)
    fig.subplots_adjust(top=0.90)

    fig.savefig(OUT / "fig_sf_coalition_composition.pdf")
    fig.savefig(OUT / "fig_sf_coalition_composition.png", dpi=450)
    fig.savefig(OUT / "fig_sf_coalition_composition.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
