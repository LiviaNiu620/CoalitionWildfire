#!/usr/bin/env python3
"""Plot the three-regime Sioux Falls response decomposition from real JSON."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "sf_response_decomposition_results.json"
SUMMARY = ROOT / "sf_response_decomposition_summary.json"
PDF = Path(__file__).resolve().parent / "fig_sf_response_decomposition.pdf"
PNG = Path(__file__).resolve().parent / "fig_sf_response_decomposition.png"

REGIMES = (
    ("mass_linked_prior", "mass-linked + prior", "#D55E00"),
    ("full_information", r"common $\beta=1$", "#0072B2"),
    ("common_half", r"common $\beta=0.5$", "#009E73"),
)


def decile_medians(x, y):
    edges = np.quantile(x, np.linspace(0.0, 1.0, 11))
    bins = np.clip(np.digitize(x, edges[1:-1]), 0, 9)
    points = []
    for index in range(10):
        mask = bins == index
        if np.any(mask):
            points.append((np.median(x[mask]), np.median(y[mask])))
    return np.asarray(points)


def main():
    raw = json.loads(RAW.read_text())
    summary = json.loads(SUMMARY.read_text())
    if len(raw["profiles"]) != 600:
        raise RuntimeError("expected 600 certified regime-profile rows")
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8.5,
        "axes.titlesize": 9.4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.14,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    all_recovery = np.array([row["recovery"] for row in raw["profiles"]]) * 100
    padding = max(0.8, 0.04 * np.ptp(all_recovery))
    limits = (all_recovery.min() - padding, all_recovery.max() + padding)
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.45), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (key, title, color) in zip(axes, REGIMES):
        rows = [row for row in raw["profiles"] if row["regime"] == key]
        if len(rows) != 200 or not all(row["certified"] for row in rows):
            raise RuntimeError(f"incomplete or uncertified regime {key}")
        hhi = np.array([row["HHI"] for row in rows])
        recovery = 100 * np.array([row["recovery"] for row in rows])
        ax.scatter(hhi, recovery, s=10, alpha=0.48, color=color,
                   edgecolors="none", rasterized=True)
        medians = decile_medians(hhi, recovery)
        ax.plot(medians[:, 0], medians[:, 1], color="#222222", linewidth=1.15,
                marker="s", markersize=2.6, label="HHI-decile median")
        stats = summary["by_regime"][key]
        inversion = 100 * stats["pairwise"][0]["inversion_rate"]
        ax.text(0.97, 0.05,
                f"Spearman = {stats['spearman']:.3f}\ninversions = {inversion:.1f}%",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.4,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8,
                      "pad": 1.8})
        ax.set_title(title)
        ax.set_xlabel("ownership concentration (HHI)")
        ax.set_ylim(*limits)
    axes[0].set_ylabel("recoverable congestion gap closed (%)")
    axes[0].legend(loc="upper left", fontsize=7.2)
    fig.savefig(PDF)
    fig.savefig(PNG)
    plt.close(fig)
    print(f"wrote {PDF}")
    print(f"wrote {PNG}")


if __name__ == "__main__":
    main()
