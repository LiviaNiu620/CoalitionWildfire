#!/usr/bin/env python3
"""Plot rank association and HHI-order reversals across fleet penetration."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "sf_alpha_screening_summary.json"
PDF = Path(__file__).resolve().parent / "fig_sf_alpha_screening.pdf"
PNG = Path(__file__).resolve().parent / "fig_sf_alpha_screening.png"


def main():
    summary = json.loads(SUMMARY.read_text())
    alphas = np.array(sorted(float(value) for value in summary["by_alpha"]))
    groups = [summary["by_alpha"][f"{alpha:g}"] for alpha in alphas]
    if any(group["certified_n"] != 100 for group in groups):
        raise RuntimeError("alpha scan is incomplete or contains uncertified profiles")
    spearman = np.array([group["spearman"] for group in groups])
    inversion = 100 * np.array([group["pairwise"][0]["inversion_rate"]
                                for group in groups])
    wide_inversion = 100 * np.array([group["pairwise"][2]["inversion_rate"]
                                     for group in groups])
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8.5,
        "axes.titlesize": 9.4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.16,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.45), constrained_layout=True)
    axes[0].plot(alphas, spearman, color="#0072B2", marker="o", markersize=4)
    axes[0].set_ylabel("Spearman correlation")
    axes[0].set_title("HHI–recovery rank association")
    axes[1].plot(alphas, inversion, color="#D55E00", marker="o", markersize=4,
                 label="all HHI gaps")
    axes[1].plot(alphas, wide_inversion, color="#009E73", marker="s", markersize=4,
                 linestyle="--", label=r"$|\Delta\mathrm{HHI}|\geq0.05$")
    axes[1].set_ylabel("pairwise reversals (%)")
    axes[1].set_title("Ownership-order reversals")
    axes[1].legend(fontsize=7.2)
    for ax in axes:
        ax.set_xlabel(r"fleet-controlled demand share $\alpha$")
        ax.set_xticks(alphas)
    fig.savefig(PDF)
    fig.savefig(PNG)
    plt.close(fig)
    print(f"wrote {PDF}")
    print(f"wrote {PNG}")


if __name__ == "__main__":
    main()
