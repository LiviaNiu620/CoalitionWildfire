#!/usr/bin/env python3
"""Generate the Paper 1 ownership--response mechanism diagram.

The diagram is deliberately produced as vector graphics with Matplotlib rather
than by a generative image model so every label and relation is reproducible.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUT = Path(__file__).resolve().parent


def box(ax, xy, width, height, text, face, edge="#4B5563", lw=0.9,
        fontsize=9, linestyle="-"):
    patch = FancyBboxPatch(
        xy, width, height,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        facecolor=face, edgecolor=edge, linewidth=lw, linestyle=linestyle,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text,
            ha="center", va="center", fontsize=fontsize, linespacing=1.25)
    return patch


def arrow(ax, start, end, color="#6B7280", lw=1.15, linestyle="-",
          mutation=11):
    a = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=mutation,
                        linewidth=lw, linestyle=linestyle, color=color,
                        shrinkA=2, shrinkB=2, connectionstyle="arc3,rad=0")
    ax.add_patch(a)
    return a


def main():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })

    fig, ax = plt.subplots(figsize=(7.25, 3.35))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    blue = "#DCEAF7"
    amber = "#F7E7C2"
    green = "#DCEFE5"
    violet = "#E8E2F3"
    gray = "#F3F4F6"
    dark_blue = "#0072B2"
    vermillion = "#D55E00"

    box(ax, (0.025, 0.64), 0.17, 0.20,
        "Ownership / control\npartition\n$(m_k,d_w^k)$", blue, edge=dark_blue)
    box(ax, (0.025, 0.25), 0.17, 0.20,
        "Information\narchitecture\n$\\widehat\\theta_e^k$", amber, edge=vermillion)

    box(ax, (0.24, 0.46), 0.19, 0.25,
        "Fleet-specific constrained\nresponse\nactive supports + OD projection",
        green, edge="#009E73")
    box(ax, (0.48, 0.46), 0.16, 0.25,
        "Aggregate fleet\nresponse signature\n$(\\mathcal{R}_F,d_F)$",
        violet, edge="#6B5FA7")
    box(ax, (0.69, 0.62), 0.14, 0.19,
        "HDV Wardrop\nblocks\n$(\\Lambda_H,\\Gamma_H,d^H)$",
        blue, edge=dark_blue)
    box(ax, (0.69, 0.27), 0.14, 0.19,
        "Activation boundary\nHDV backfilling\nvs. corner support",
        amber, edge=vermillion)
    box(ax, (0.87, 0.43), 0.105, 0.24,
        "Joint network\nequilibrium\nand congestion",
        green, edge="#009E73", fontsize=8.5)

    arrow(ax, (0.195, 0.74), (0.24, 0.61), color=dark_blue)
    arrow(ax, (0.195, 0.35), (0.24, 0.55), color=vermillion)
    arrow(ax, (0.43, 0.585), (0.48, 0.585))
    arrow(ax, (0.64, 0.585), (0.69, 0.70))
    arrow(ax, (0.64, 0.54), (0.69, 0.37))
    arrow(ax, (0.83, 0.70), (0.87, 0.58), color=dark_blue)
    arrow(ax, (0.83, 0.37), (0.87, 0.50), color=vermillion)

    box(ax, (0.045, 0.025), 0.13, 0.105,
        "HHI\nshare-only view", gray, edge="#6B7280", fontsize=8,
        linestyle="--")
    hhi_arrow = FancyArrowPatch(
        (0.025, 0.72), (0.08, 0.13), arrowstyle="-|>",
        mutation_scale=9, linewidth=0.9, linestyle="--", color="#6B7280",
        connectionstyle="arc3,rad=0.34", shrinkA=2, shrinkB=2,
    )
    ax.add_patch(hhi_arrow)
    ax.text(0.205, 0.105,
            "HHI observes the ownership branch but omits information,\n"
            "support-dependent substitution, and endogenous HDV adjustment.",
            ha="left", va="center", fontsize=8.1, color="#4B5563")

    ax.text(0.5, 0.94,
            "Ownership and interoperability act on different parts of the routing mechanism",
            ha="center", va="center", fontsize=11, fontweight="bold")

    fig.savefig(OUT / "fig_response_mechanism.pdf")
    fig.savefig(OUT / "fig_response_mechanism.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
