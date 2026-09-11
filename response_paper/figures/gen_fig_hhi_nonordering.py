#!/usr/bin/env python3
"""Generate the Paper 1 two-fleet HHI non-ordering figure.

The plot uses the verified canonical bottleneck closed form.  It is kept
separate from the legacy ``fig_hhi_main.pdf`` because that filename has been
produced by several older scripts with incompatible data and titles.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


OUTPUT_DIR = Path(__file__).resolve().parent
LAMBDA = 2.0
X_UE = 1.0


def beta(mass: float) -> float:
    return 1.0 - np.exp(-LAMBDA * mass)


def inverse_response(mass: float) -> float:
    return 1.0 / beta(mass)


def equilibrium(masses: list[float]) -> tuple[float, int]:
    """Return recovery and the number of interior fleets."""
    masses = [mass for mass in masses if mass > 1e-9]
    interior = list(range(len(masses)))
    for _ in range(80):
        sigma = sum(inverse_response(masses[k]) for k in interior)
        corner_flow = sum(
            masses[k] for k in range(len(masses)) if k not in interior
        )
        x_1 = (corner_flow + sigma * X_UE) / (1.0 + sigma)
        gap = X_UE - x_1
        updated = [
            k
            for k in range(len(masses))
            if beta(masses[k]) * masses[k] > gap + 1e-12
        ]
        if set(updated) == set(interior):
            break
        if updated:
            interior = updated

    sigma = sum(inverse_response(masses[k]) for k in interior)
    corner_flow = sum(
        masses[k] for k in range(len(masses)) if k not in interior
    )
    x_1 = (corner_flow + sigma * X_UE) / (1.0 + sigma)
    recovery = 4.0 * x_1 * (1.0 - x_1)
    return recovery, len(interior)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "legend.frameon": False,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.18,
            "lines.linewidth": 2.0,
        }
    )

    deltas = np.linspace(0.0, 0.49, 10_000)
    hhi = (0.5 + deltas) ** 2 + (0.5 - deltas) ** 2
    recovery = np.empty_like(deltas)
    n_interior = np.empty_like(deltas, dtype=int)
    for index, delta in enumerate(deltas):
        recovery[index], n_interior[index] = equilibrium(
            [0.5 + delta, 0.5 - delta]
        )

    switch = next(
        i for i in range(1, len(n_interior))
        if n_interior[i] < n_interior[i - 1]
    )
    trough = int(np.argmin(recovery[: switch + 3]))

    vermillion = "#D55E00"
    blue = "#0072B2"
    gray = "#6B7280"

    fig, ax = plt.subplots(figsize=(6.4, 3.7))
    ax.plot(
        hhi[: switch + 1],
        recovery[: switch + 1],
        color=vermillion,
        label="fixed interior support",
    )
    ax.plot(
        hhi[switch:],
        recovery[switch:],
        color=blue,
        label="after support change",
    )
    ax.axvline(hhi[switch], color=gray, linestyle="--", linewidth=1.1)
    ax.scatter(
        [hhi[0], hhi[trough], hhi[switch], hhi[-1]],
        [recovery[0], recovery[trough], recovery[switch], recovery[-1]],
        color=[vermillion, vermillion, blue, blue],
        edgecolor="white",
        linewidth=0.6,
        s=[34, 34, 34, 48],
        zorder=4,
    )
    ax.annotate(
        "smaller fleet reaches\na route corner",
        xy=(hhi[switch], recovery[switch]),
        xytext=(0.22, 0.13),
        textcoords="axes fraction",
        fontsize=8.5,
        ha="left",
        arrowprops={"arrowstyle": "->", "color": gray, "lw": 1.0},
    )
    ax.set_xlabel("ownership concentration (HHI)")
    ax.set_ylabel(r"recovery ratio $\rho$")
    ax.set_title("HHI does not order recovery along a two-fleet path")
    ax.set_ylim(0.70, 1.01)
    ax.legend(loc="lower right")

    inset = inset_axes(ax, width="38%", height="35%", loc="upper left", borderpad=1.3)
    inset.plot(
        hhi[: switch + 1], recovery[: switch + 1], color=vermillion, linewidth=1.7
    )
    inset.scatter(
        hhi[trough], recovery[trough], color=vermillion, s=22, zorder=3
    )
    inset.set_title("fixed-support dip", fontsize=8)
    inset.tick_params(labelsize=7)
    inset.grid(alpha=0.18)
    mark_inset(ax, inset, loc1=3, loc2=4, fc="none", ec="#9CA3AF", lw=0.7)

    pdf_path = OUTPUT_DIR / "fig_hhi_nonordering.pdf"
    png_path = OUTPUT_DIR / "fig_hhi_nonordering.png"
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)

    print(f"saved {pdf_path}")
    print(f"saved {png_path}")
    print(
        "symmetric={:.6f}, trough={:.6f}, support_change={:.6f}, near_monopoly={:.6f}".format(
            recovery[0], recovery[trough], recovery[switch], recovery[-1]
        )
    )
    print(
        "hhi_symmetric={:.6f}, hhi_trough={:.6f}, hhi_support_change={:.6f}, hhi_near_monopoly={:.6f}".format(
            hhi[0], hhi[trough], hhi[switch], hhi[-1]
        )
    )


if __name__ == "__main__":
    main()
