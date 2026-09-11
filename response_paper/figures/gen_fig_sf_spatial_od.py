#!/usr/bin/env python3
"""Plot fixed-HHI Sioux Falls outcomes against spatial OD specialization."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "sf_spatial_od_results.json"
SUMMARY = ROOT / "sf_spatial_od_summary.json"
PDF = Path(__file__).resolve().parent / "fig_sf_spatial_od.pdf"
PNG = Path(__file__).resolve().parent / "fig_sf_spatial_od.png"

FAMILIES = (
    ("origin", "origin exposure", "#0072B2"),
    ("destination", "destination exposure", "#E69F00"),
    ("joint_od", "joint OD exposure", "#009E73"),
    ("corridor", "corridor incidence", "#CC79A7"),
)
MARKERS = {0.35: "o", 0.75: "s"}


def main():
    raw = json.loads(RAW.read_text())
    summary = json.loads(SUMMARY.read_text())
    if len(raw["profiles"]) != 82:
        raise RuntimeError("expected 82 fixed-HHI profile rows")
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
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 2.65), sharex=True,
                             constrained_layout=True)
    for ax, alpha in zip(axes, (0.5, 0.9)):
        rows = [row for row in raw["profiles"] if row["alpha"] == alpha]
        if len(rows) != 41 or not all(row["certified"] for row in rows):
            raise RuntimeError(f"incomplete or uncertified alpha={alpha}")
        baseline = next(row for row in rows if row["family"] == "proportional")
        for family, label, color in FAMILIES:
            family_rows = [row for row in rows if row["family"] == family]
            for temperature, marker in MARKERS.items():
                selected = [row for row in family_rows
                            if np.isclose(row["temperature"], temperature)]
                ax.scatter(
                    [1 - row["weighted_OD_entropy"] for row in selected],
                    [100 * row["recovery"] for row in selected],
                    color=color, marker=marker, s=27, alpha=0.8,
                    edgecolor="white", linewidth=0.35,
                    label=label if np.isclose(temperature, 0.35) else None,
                )
        ax.axhline(100 * baseline["recovery"], color="#333333", linestyle="--",
                   linewidth=1.0, label="proportional baseline")
        stats = summary["by_alpha"][f"{alpha:g}"]
        spread = stats["spatial_recovery"]["max"] - stats["spatial_recovery"]["min"]
        ax.text(0.97, 0.05,
                f"range = {100 * spread:.2f} pp\nmax |ΔJ| = "
                f"{stats['max_abs_delta_J_pct_UE']:.3f}% of UE",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.3,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                      "pad": 1.8})
        ax.set_title(rf"fleet penetration $\alpha={alpha:.1f}$; HHI $=0.25$")
        ax.set_xlabel("spatial OD specialization (1 − normalized entropy)")
    axes[0].set_ylabel("recoverable congestion gap closed (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels, loc="best", fontsize=6.9)
    fig.savefig(PDF)
    fig.savefig(PNG)
    plt.close(fig)
    print(f"wrote {PDF}")
    print(f"wrote {PNG}")


if __name__ == "__main__":
    main()
