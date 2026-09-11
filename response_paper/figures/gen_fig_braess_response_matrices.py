#!/usr/bin/env python3
"""Export the Braess fleet-response and HDV-screened response matrices.

The diagnostic uses the affine Braess network and the three routes reported in
the manuscript.  Two fleets use all three routes with constant internalization
coefficients beta=(0.8, 0.5).  The raw matrix has units of inverse link slope;
the plotted matrices divide by their common scale so that the substitution
pattern, rather than a unit convention, is visible.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


OUTPUT_DIR = Path(__file__).resolve().parent
DATA_PATH = OUTPUT_DIR / "fig_braess_response_matrices_data.json"
PDF_PATH = OUTPUT_DIR / "fig_braess_response_matrices.pdf"
PNG_PATH = OUTPUT_DIR / "fig_braess_response_matrices.png"

EDGE_LABELS = [r"$sv$", r"$wt$", r"$sw$", r"$vt$", r"$vw$"]
ROUTES = ((0, 3), (2, 1), (0, 4, 1))
BETAS = np.array([0.8, 0.5])
SLOPES = np.array([0.01, 0.01, 0.0, 0.0, 0.0])


def incidence() -> np.ndarray:
    matrix = np.zeros((len(EDGE_LABELS), len(ROUTES)))
    for route_index, edges in enumerate(ROUTES):
        matrix[list(edges), route_index] = 1.0
    return matrix


def constrained_block(lam: np.ndarray, beta: float) -> np.ndarray:
    """Return the upper-left inverse block of the one-OD bordered KKT system."""
    hessian = lam.T @ np.diag(beta * SLOPES) @ lam
    gamma = np.ones((1, lam.shape[1]))
    bordered = np.block([
        [hessian, -gamma.T],
        [gamma, np.zeros((1, 1))],
    ])
    response = np.linalg.inv(bordered)[: lam.shape[1], : lam.shape[1]]
    np.testing.assert_allclose(gamma @ response, 0.0, atol=1e-11)
    np.testing.assert_allclose(response @ gamma.T, 0.0, atol=1e-11)
    return response


def screening_operator(lam_h: np.ndarray) -> np.ndarray:
    """Construct A_H from the retained tangent basis of an HDV route support."""
    route_count = lam_h.shape[1]
    if route_count == 3:
        tangent_basis = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]])
    elif route_count == 2:
        tangent_basis = np.array([[1.0], [-1.0]])
    else:
        raise ValueError("this diagnostic expects a two- or three-route HDV support")
    edge_tangent = lam_h @ tangent_basis
    b_matrix = np.diag(SLOPES)
    gram = edge_tangent.T @ b_matrix @ edge_tangent
    return (
        np.eye(len(EDGE_LABELS))
        - edge_tangent @ np.linalg.inv(gram) @ edge_tangent.T @ b_matrix
    )


def matrices():
    lam = incidence()
    route_response = sum(constrained_block(lam, beta) for beta in BETAS)
    fleet_response = lam @ route_response @ lam.T

    screen_three = screening_operator(lam)
    screen_outer = screening_operator(lam[:, :2])
    activated_three = screen_three @ fleet_response
    activated_outer = screen_outer @ fleet_response

    scale = float(np.max(np.abs(fleet_response)))
    normalized = {
        "fleet": fleet_response / scale,
        "three": activated_three / scale,
        "outer": activated_outer / scale,
    }
    np.testing.assert_allclose(activated_three, 0.0, atol=1e-10)
    if np.linalg.matrix_rank(activated_outer, tol=1e-10) != 1:
        raise RuntimeError("two-route activated response should have rank one")
    return {
        "lambda": lam,
        "route_response": route_response,
        "fleet_response": fleet_response,
        "screen_three": screen_three,
        "screen_outer": screen_outer,
        "activated_three": activated_three,
        "activated_outer": activated_outer,
        "scale": scale,
        "normalized": normalized,
    }


def export_data(values) -> None:
    payload = {
        "network": "affine Braess",
        "edge_order": ["sv", "wt", "sw", "vt", "vw"],
        "route_edge_indices": [list(route) for route in ROUTES],
        "fleet_betas": BETAS.tolist(),
        "edge_slopes": SLOPES.tolist(),
        "normalization_scale": values["scale"],
        "R_F": values["fleet_response"].tolist(),
        "A_H_three_routes": values["screen_three"].tolist(),
        "A_H_R_F_three_routes": values["activated_three"].tolist(),
        "A_H_outer_routes": values["screen_outer"].tolist(),
        "A_H_R_F_outer_routes": values["activated_outer"].tolist(),
        "rank_R_F": int(np.linalg.matrix_rank(values["fleet_response"], tol=1e-10)),
        "rank_A_H_R_F_three_routes": int(
            np.linalg.matrix_rank(values["activated_three"], tol=1e-10)
        ),
        "rank_A_H_R_F_outer_routes": int(
            np.linalg.matrix_rank(values["activated_outer"], tol=1e-10)
        ),
    }
    DATA_PATH.write_text(json.dumps(payload, indent=2) + "\n")


def annotate_heatmap(ax, matrix: np.ndarray) -> None:
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            color = "white" if abs(value) >= 0.62 else "#222222"
            label = "0" if abs(value) < 5e-12 else f"{value:.2g}"
            ax.text(column, row, label, ha="center", va="center", fontsize=7.3,
                    color=color)


def plot(values) -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8.5,
        "axes.titlesize": 9.2,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    matrices_to_plot = [
        values["normalized"]["fleet"],
        values["normalized"]["three"],
        values["normalized"]["outer"],
    ]
    titles = [
        r"Fleet response $\mathcal{R}_F$" + "\n(rank 2)",
        r"$\mathcal{A}_H\mathcal{R}_F$: three HDV routes" + "\n(rank 0)",
        r"$\mathcal{A}_H\mathcal{R}_F$: outer HDV routes" + "\n(rank 1)",
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.55), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    image = None
    for ax, matrix, title in zip(axes, matrices_to_plot, titles):
        image = ax.imshow(matrix, cmap="RdBu_r", norm=norm, interpolation="nearest")
        annotate_heatmap(ax, matrix)
        ax.set_title(title, pad=5)
        ax.set_xticks(range(len(EDGE_LABELS)), EDGE_LABELS)
        ax.set_yticks(range(len(EDGE_LABELS)), EDGE_LABELS)
        ax.set_xlabel("cost perturbation edge")
        ax.tick_params(length=0)
    axes[0].set_ylabel("fleet-flow response edge")
    for ax in axes[1:]:
        ax.set_yticklabels([])
    colorbar = fig.colorbar(image, ax=axes, shrink=0.82, pad=0.012, aspect=22)
    colorbar.set_label(r"entry divided by $\max |\mathcal{R}_F|$")
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH)
    plt.close(fig)


def main() -> None:
    values = matrices()
    export_data(values)
    plot(values)
    print(f"wrote {DATA_PATH}")
    print(f"wrote {PDF_PATH}")
    print(f"wrote {PNG_PATH}")
    print(f"normalization scale = {values['scale']:.6g}")


if __name__ == "__main__":
    main()
