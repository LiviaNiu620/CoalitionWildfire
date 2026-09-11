"""Generate compact figures for the coalition sensitivity and partition scans."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})

def sensitivity():
    data = json.loads((ROOT / "sf_coalition_sensitivity_summary.json").read_text())
    rows = data["assortative_minus_mixed"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    colors = {"narrow": "#4472C4", "baseline": "#ED7D31", "wide": "#70AD47"}
    for alpha, ax in zip((0.5, 0.9), axes):
        for design in ("narrow", "baseline", "wide"):
            for scale in (0.5, 1.0, 2.0):
                subset = [r for r in rows if r["alpha"] == alpha and r["type_design"] == design and r["slope_scale"] == scale]
                subset.sort(key=lambda r: r["gamma"])
                ax.plot([r["gamma"] for r in subset], [r["delta_recovery_percentage_points"] for r in subset],
                        marker="o", ms=2.8, lw=1.0, alpha=0.42 + 0.18 * scale,
                        color=colors[design], label=f"{design}, scale={scale:g}")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_title(rf"$\alpha={alpha}$")
        ax.set_xlabel(r"coordination $\gamma$")
        ax.set_ylabel("assortative minus mixed recovery (pp)")
        ax.grid(alpha=0.2)
    axes[1].legend(fontsize=6, ncol=2, frameon=False, loc="lower left")
    fig.savefig(OUT / "fig_coalition_sensitivity.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_coalition_sensitivity.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

def partitions():
    data = json.loads((ROOT / "sf_coalition_partition_scan_summary.json").read_text())
    groups = data["by_group"]
    labels = [f"a={r['alpha']}, g={r['gamma']}" for r in groups]
    rho = [r["spearman_coalition_HHI_recovery"] for r in groups]
    inv = [100 * r["pairwise_inversion_rate"] for r in groups]
    fig, ax = plt.subplots(figsize=(5.8, 3.0), constrained_layout=True)
    x = np.arange(len(labels)); w = 0.36
    ax.bar(x-w/2, rho, w, label="Spearman", color="#4472C4")
    ax.bar(x+w/2, inv, w, label="inversion rate (%)", color="#C0504D")
    ax.set_xticks(x, labels); ax.set_ylabel("value"); ax.set_ylim(0, max(max(rho), max(inv))*1.2)
    ax.set_title("All 15 n=4 coalition partitions")
    ax.grid(axis="y", alpha=0.2); ax.legend(frameon=False, fontsize=8)
    fig.savefig(OUT / "fig_coalition_partition_scan.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_coalition_partition_scan.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    sensitivity(); partitions()
