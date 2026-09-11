"""Company-level accounting costs for coalition-formation experiments.

The equilibrium operator contains coalition cross terms, but those terms do not
have a unique member-level allocation.  Formation payoffs therefore use the
declared accounting cost

    l_i = sum_e f_i,e c_e(x) + 0.5 sum_e bbar_e tau_i f_i,e^2,

and never split a coalition cross term among members.
"""
from __future__ import annotations

import numpy as np

import sf_equilibrium as SE


def company_accounting_costs(sol, game):
    """Return auditable company travel, management, and total costs."""
    if game.coalitions is None:
        raise ValueError("company accounting costs require coalition governance")
    own_edge = np.asarray(sol["fF"], dtype=float) @ game.A.T
    x = np.asarray(sol["x"], dtype=float)
    physical_cost = SE.bpr_cost(x, game.t0, game.cap)
    travel = own_edge @ physical_cost
    slope = game.management_slope
    if game.management_slope_mode == "state":
        slope = SE.bpr_deriv(x, game.t0, game.cap)
    management = 0.5 * np.sum(
        slope[None, :]
        * game.objective_types[:, None]
        * own_edge ** 2,
        axis=1,
    )
    total = travel + management
    if not (
        np.all(np.isfinite(travel))
        and np.all(np.isfinite(management))
        and np.all(np.isfinite(total))
        and np.all(total >= 0.0)
    ):
        raise FloatingPointError("invalid company accounting cost")
    return {
        "company_travel_cost": travel.tolist(),
        "company_diagonal_management_cost": management.tolist(),
        "company_accounting_cost": total.tolist(),
        "total_av_travel_cost": float(np.sum(travel)),
        "total_company_accounting_cost": float(np.sum(total)),
    }


__all__ = ["company_accounting_costs"]
