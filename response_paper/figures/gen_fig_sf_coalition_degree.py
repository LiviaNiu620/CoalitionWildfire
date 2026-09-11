"""Plot coalition degree results from the certified company-count scan."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent


def main():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.2, "axes.titlesize": 10.5, "axes.titleweight": "bold",
        "axes.labelsize": 9.4, "legend.fontsize": 8.0, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.18, "savefig.dpi": 450,
        "savefig.bbox": "tight",
    })
    raw = json.loads((ROOT / "sf_coalition_degree_results.json").read_text())
    rows = raw["profiles"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
    colors = {2: "#0072B2", 4: "#D55E00", 8: "#009E73"}
    labels = {2: "2 companies", 4: "4 companies", 8: "8 companies"}
    for ax, alpha in zip(axes, (0.5, 0.9)):
        for n in (2, 4, 8):
            vals = [
                (r["n_coalitions"], 100.0 * r["recovery"])
                for r in rows
                if r["alpha"] == alpha and r["n_companies"] == n and r["gamma"] == 1.0
            ]
            vals.sort()
            ax.plot([k for k, _ in vals], [v for _, v in vals], marker="o",
                    color=colors[n], label=labels[n])
        ax.set_title(rf"AV penetration $\alpha={alpha:.1f}$")
        ax.set_xlabel("number of active coalitions $K$")
        ax.set_xlim(0.85, 8.4)
        ax.set_xscale("log", base=2)
        ax.set_xticks([1, 2, 4, 8], labels=["1", "2", "4", "8"])
    axes[0].set_ylabel("recoverable congestion gap closed (%)")
    axes[0].legend(loc="lower left")
    fig.suptitle(r"Coalition degree at full cross-member coordination ($\gamma=1$)", y=1.04, fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "fig_sf_coalition_degree.pdf")
    fig.savefig(OUT / "fig_sf_coalition_degree.png", dpi=450)
    fig.savefig(OUT / "fig_sf_coalition_degree.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
