"""Generate the Sioux Falls ownership-simplex scan figure from committed data."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sf_simplex_scan_results.json"
SUMMARY = ROOT / "sf_simplex_scan_summary.json"
OUTPUT = Path(__file__).with_name("fig_sf_simplex_scan")
COLORS = {
    0.2: "#CC3311",
    0.5: "#EE7733",
    1.0: "#009988",
    2.0: "#0077BB",
    5.0: "#33BBEE",
}


def largest_inversion(hhi, recovery, travel_time, threshold=0.05):
    i, j = np.triu_indices(len(hhi), 1)
    eligible = (np.abs(hhi[i] - hhi[j]) >= threshold) & (
        (hhi[i] - hhi[j]) * (recovery[i] - recovery[j]) < 0
    )
    candidates = np.flatnonzero(eligible)
    winner = candidates[np.argmax(np.abs(travel_time[i[candidates]] - travel_time[j[candidates]]))]
    return i[winner], j[winner]


def main():
    raw = json.loads(DATA.read_text())
    summary = json.loads(SUMMARY.read_text())
    rows = raw["profiles"]
    hhi = np.array([row["HHI"] for row in rows])
    recovery = 100 * np.array([row["recovery"] for row in rows])
    travel_time = np.array([row["J"] for row in rows])
    concentration = np.array([row["dirichlet_a"] for row in rows])

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "savefig.dpi": 450,
        "svg.fonttype": "none",
    })

    fig, ax = plt.subplots(figsize=(7.0, 4.35), constrained_layout=True)
    for a in COLORS:
        mask = concentration == a
        ax.scatter(
            hhi[mask], recovery[mask], s=25, alpha=0.76,
            color=COLORS[a], edgecolor="white", linewidth=0.35,
            label=rf"Dirichlet $a={a:g}$",
        )

    # Decile medians describe the strong average association without imposing a fit.
    edges = np.quantile(hhi, np.linspace(0, 1, 11))
    bin_id = np.clip(np.digitize(hhi, edges[1:-1]), 0, 9)
    median_hhi = np.array([np.median(hhi[bin_id == b]) for b in range(10)])
    median_recovery = np.array([np.median(recovery[bin_id == b]) for b in range(10)])
    ax.plot(
        median_hhi, median_recovery, color="#333333", linewidth=1.2,
        marker="s", markersize=3.2, label="HHI-decile median",
    )

    first, second = largest_inversion(hhi, recovery, travel_time)
    ax.plot(
        hhi[[first, second]], recovery[[first, second]], color="#222222",
        linestyle="--", linewidth=1.0, zorder=4,
    )
    ax.scatter(
        hhi[[first, second]], recovery[[first, second]], s=55,
        facecolor="none", edgecolor="#222222", linewidth=1.2, zorder=5,
    )
    midpoint = (
        0.5 * (hhi[first] + hhi[second]),
        0.5 * (recovery[first] + recovery[second]),
    )
    ax.annotate(
        "largest inversion with HHI gap >= 0.05",
        xy=midpoint, xytext=(midpoint[0] + 0.08, midpoint[1] - 3.0),
        arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.8},
        fontsize=8, color="#333333",
    )

    all_pairs = summary["pairwise"][0]
    wide_pairs = summary["pairwise"][2]
    stats_text = (
        f"Spearman rho = {summary['spearman']:.3f}\n"
        f"all-pair inversions: {100*all_pairs['inversion_rate']:.1f}%\n"
        "inversions for HHI gap >= 0.05: "
        f"{100*wide_pairs['inversion_rate']:.1f}%"
    )
    ax.text(
        0.985, 0.025, stats_text, transform=ax.transAxes,
        ha="right", va="bottom", fontsize=8,
        bbox={"boxstyle": "square,pad=0.35", "facecolor": "white",
              "edgecolor": "#888888", "linewidth": 0.6, "alpha": 0.94},
    )

    ax.set_xlabel("ownership concentration (HHI)")
    ax.set_ylabel("recovery of the UE-SO gap (%)")
    ax.set_xlim(0.235, 1.015)
    ax.set_ylim(68.0, 92.0)
    ax.legend(loc="upper left", ncol=2, columnspacing=0.9, handletextpad=0.4)

    fig.savefig(OUTPUT.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
