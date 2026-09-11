"""Finite-difference certificate for the fixed-reference coalition potential."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sf_equilibrium import SFGame, RouteSet


OUT = Path("potential_structure_check.json")


def run_check() -> dict:
    sf = {
        "edges": 2,
        "t0": np.array([1.0, 1.4]),
        "cap": np.array([100.0, 100.0]),
        "od_list": [{"demand": 10.0, "paths": [[0], [1]]}],
    }
    game = SFGame(
        sf,
        alpha=0.6,
        shares=[0.5, 0.5],
        lam_beta=0.1,
        routes=RouteSet(sf["od_list"]),
        mode="eq",
        coalitions=((0, 1),),
        objective_types=[0.5, 1.5],
        coordination_gamma=0.7,
        management_slope=[0.8, 1.1],
        public_edge_penalty=[0.2, 0.1],
    )
    f_h, f_f = game.initial(seed=0)
    f_h[:] = [2.0, 2.0]
    f_f[0, :] = [1.5, 1.5]
    f_f[1, :] = [1.5, 1.5]

    def pack(h, fleets):
        return np.concatenate([h, fleets.ravel()])

    def unpack(z):
        return z[: game.P].copy(), z[game.P :].reshape(game.K, game.P).copy()

    def operator(z):
        h, fleets = unpack(z)
        g_h, g_f, *_ = game.operator(h, fleets)
        return pack(g_h, g_f)

    def potential(z):
        h, fleets = unpack(z)
        x = game.edge_flow(h, fleets)
        physical = np.sum(
            game.t0 * x + game.t0 * 0.15 * x**5 / (5.0 * game.cap**4)
            + game.public_edge_penalty * x
        )
        own = fleets @ game.A.T
        management = 0.0
        for block in game.coalitions:
            idx = np.asarray(block)
            governance = game.coalition_governance_matrix(block)
            management += 0.5 * np.sum(
                game.management_slope
                * np.einsum("ie,ij,je->e", own[idx], governance, own[idx])
            )
        return float(physical + management)

    z = pack(f_h, f_f)
    eps = 1e-6
    n = z.size
    jacobian = np.empty((n, n))
    potential_gradient = np.empty(n)
    for j in range(n):
        direction = np.zeros(n)
        direction[j] = eps
        jacobian[:, j] = (operator(z + direction) - operator(z - direction)) / (2 * eps)
        potential_gradient[j] = (potential(z + direction) - potential(z - direction)) / (2 * eps)

    symmetry_error = float(np.max(np.abs(jacobian - jacobian.T)))
    gradient_error = float(np.max(np.abs(potential_gradient - operator(z))))
    min_symmetric_eigenvalue = float(np.linalg.eigvalsh((jacobian + jacobian.T) / 2.0).min())
    result = {
        "model": "fixed-reference symmetric PSD coalition governance",
        "instance": "two parallel edges, one OD, two companies, one coalition",
        "jacobian_shape": list(jacobian.shape),
        "max_jacobian_symmetry_error": symmetry_error,
        "max_potential_gradient_error": gradient_error,
        "min_symmetric_part_eigenvalue": min_symmetric_eigenvalue,
        "passed": bool(symmetry_error <= 1e-7 and gradient_error <= 1e-7),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    print(json.dumps(run_check(), indent=2))
