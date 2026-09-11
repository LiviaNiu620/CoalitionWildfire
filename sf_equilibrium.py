"""
Sioux Falls equilibrium solver with convergence-certified numerics.

Replaces the fixed-iteration projected-gradient loop in ``sf_report.py``.
Three deficiencies of the old code are addressed (reviewer Major-9):

  9.1  The route set was frozen at the three shortest free-flow paths per OD.
       An omitted route can become cheaper once flows build up, so the old
       solution was an equilibrium of the restricted game only.  We now run
       column generation with an exact reduced-cost certificate.

  9.2  The iteration was a fixed-step simultaneous projected-gradient scheme
       run for a fixed 6000 iterations with no stopping test, described in the
       paper as "projected best response ... stop when the update falls below
       1e-6".  Neither the method nor the stopping rule matched.  We now use
       Korpelevich extragradient with an adaptive step and a gap-based test.

  9.3  The reported "VI gap" only checked the HDV population; the fleets'
       perceived costs were never tested.  We now report the gap over every
       player, normalised, plus the omitted-route certificate.

--------------------------------------------------------------------------
Reduced-cost certificate (why one plain shortest path suffices)
--------------------------------------------------------------------------
Operator k's route gradient is

    PC_r^k = sum_{e in r} [ c_e(x_e) + beta_k * f_e^k * c'_e(x_e) ].

Every term beta_k * f_e^k * c'_e is nonnegative, hence for ANY route r,

    PC_r^k >= sum_{e in r} c_e(x_e) >= SP_w := shortest path under c_e.

So if SP_w >= lambda_{k,w} (operator k's minimum perceived cost over the
routes it currently holds), no omitted route can be profitable for k.  The
same bound with beta = 0 is exactly the Wardrop test for HDVs.  One
shortest-path computation per OD therefore certifies every player at once,
and when it fails it hands us the column to add.
"""
from __future__ import annotations

import numpy as np
import networkx as nx

from sioux_falls_loader import build_sioux_falls, parse_net, parse_trips, DATA_DIR
import os

BPR_B, BPR_ETA = 0.15, 4.0


# ----------------------------------------------------------------------
# cost functions
# ----------------------------------------------------------------------
def bpr_cost(x, t0, cap):
    return t0 * (1.0 + BPR_B * (x / cap) ** BPR_ETA)


def bpr_deriv(x, t0, cap):
    return t0 * BPR_B * BPR_ETA * (x ** (BPR_ETA - 1.0)) / (cap ** BPR_ETA)


def bpr_second_deriv(x, t0, cap):
    return (
        t0 * BPR_B * BPR_ETA * (BPR_ETA - 1.0)
        * (x ** (BPR_ETA - 2.0)) / (cap ** BPR_ETA)
    )


# ----------------------------------------------------------------------
# route-set container
# ----------------------------------------------------------------------
class RouteSet:
    """Per-OD list of edge-index paths, growable by column generation."""

    def __init__(self, od_list):
        self.od_list = od_list
        self.paths = [list(od["paths"]) for od in od_list]

    @property
    def n_od(self):
        return len(self.od_list)

    def total(self):
        return sum(len(p) for p in self.paths)

    def add(self, oi, edge_path):
        key = tuple(edge_path)
        if any(tuple(p) == key for p in self.paths[oi]):
            return False
        self.paths[oi].append(list(edge_path))
        return True

    def incidence(self, n_edges):
        """Return (A, path_od) with A[e, p] = 1 if path p uses edge e."""
        flat, path_od = [], []
        for oi, plist in enumerate(self.paths):
            for pe in plist:
                flat.append(pe)
                path_od.append(oi)
        P = len(flat)
        A = np.zeros((n_edges, P))
        for p, pe in enumerate(flat):
            for e in pe:
                A[e, p] = 1.0
        return A, np.array(path_od)

    def keys(self):
        """Flat list of (od_index, path) identifiers in incidence order."""
        return [(oi, tuple(pe)) for oi, plist in enumerate(self.paths)
                for pe in plist]


# ----------------------------------------------------------------------
# simplex projection
# ----------------------------------------------------------------------
def project_simplex(v, total):
    """Euclidean projection of v onto {f >= 0, sum f = total}."""
    if total <= 0 or v.size == 0:
        return np.zeros_like(v)
    if v.size == 1:
        return np.array([total])
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - total
    idx = np.arange(1, u.size + 1)
    cond = u - css / idx > 0
    pos = np.nonzero(cond)[0]
    r = pos[-1] if pos.size else 0
    theta = css[r] / (r + 1.0)
    return np.maximum(v - theta, 0.0)


def project_simplex_batch(V, totals):
    """
    Vectorised Euclidean projection of each row of V onto its own scaled
    simplex {f >= 0, sum f = totals[i]}.

    V      : (B, n)
    totals : (B,)

    The scalar version above is called once per (player, OD) block, which on
    Sioux Falls is 528 * (K+1) Python calls per projection and two projections
    per iteration -- roughly ten thousand interpreter round-trips per
    iteration, which dominated the runtime.  This routine does the whole batch
    with array operations.  The sorted condition u_j - (cumsum_j - total)/j > 0
    holds on a prefix, so the breakpoint index is just its count minus one.
    """
    B, n = V.shape
    if n == 1:
        return totals.reshape(B, 1).astype(float)
    U = -np.sort(-V, axis=1)
    css = np.cumsum(U, axis=1) - totals[:, None]
    idx = np.arange(1, n + 1)
    cond = U - css / idx > 0
    r = cond.sum(axis=1) - 1
    r = np.maximum(r, 0)
    theta = css[np.arange(B), r] / (r + 1.0)
    out = np.maximum(V - theta[:, None], 0.0)
    out[totals <= 0] = 0.0
    return out


# ----------------------------------------------------------------------
# the game
# ----------------------------------------------------------------------
class SFGame:
    """
    Players: K fleets (atomic, perceived cost c + beta_k f^k c') and one
    non-atomic HDV population (perceived cost c).

    ``mode`` selects the benchmark being computed:
      'eq' : the information-limited equilibrium of Def. 1
      'ue' : user equilibrium  (every player uses c_e)
      'so' : system optimum    (every player uses c_e + x_e c'_e)
    """

    def __init__(self, sf, alpha, shares, lam_beta, routes, mode="eq",
                 kappa=None, theta_bar=None, beta_override=None,
                 dem_firm_override=None, coalitions=None,
                 objective_types=None, coordination_gamma=0.0,
                 management_slope=None, management_slope_mode="reference",
                 public_edge_penalty=None, perceived_edge_penalty=None,
                 operational_pooling=False):
        """
        kappa      : prior slope multiplier.  The prior mean slope on edge e is
                     theta_bar_e = kappa * c'_e(x^UE_e), so beta0_e =
                     theta_bar_e/theta_e equals kappa when flows are at the user
                     equilibrium.  kappa < 1 is regime (R-), kappa > 1 is (R+).
                     kappa=None reproduces the old theta_bar=0 behaviour, which
                     the model does not actually admit (theta_bar>0 is required)
                     and which is retained only for regression tests.
        theta_bar  : the prior slope vector itself, if supplied directly.
        beta_override : optional scalar or length-K vector replacing the
                     mass-linked coverage coefficients.  This is used to
                     separate atomic control from the mass--response mapping.
        dem_firm_override : optional (n_od, K) OD-demand matrix.  Its row sums
                     must equal alpha times OD demand and its column sums must
                     match the fleet masses implied by ``shares``.
        coalitions : optional partition of company indices.  When supplied,
                     the companies in each block are one joint VI player but
                     retain their member-specific OD constraints.  Coalition
                     coupling enters through a positive-semidefinite governance
                     matrix; projection remains separable by member and OD.
        objective_types : positive length-K vector of company management-
                     objective curvature multipliers.
        coordination_gamma : cross-member coordination in [0, 1].  Within a
                     coalition C the member governance matrix is
                     (1-gamma) diag(type_C) +
                     gamma sqrt(type_C) sqrt(type_C)^T.
        management_slope : nonnegative length-E reference edge curvature for
                     the coalition management objective.  It is fixed during
                     the equilibrium solve, so the added objective is convex
                     quadratic even though physical BPR costs remain nonlinear.
        management_slope_mode : ``reference`` uses management_slope as a fixed
                     curvature (the paper's baseline); ``state`` uses the
                     current BPR derivative c'(x) in the coalition gradient.
                     The latter is an optional nonlinear atomic-internalization
                     diagnostic and is not a potential-game guarantee.
        public_edge_penalty : optional nonnegative public route-cost addition
                     applied to every population. It represents a declared
                     common risk/information penalty and preserves the common
                     potential when fixed during a solve.
        perceived_edge_penalty : explicit alias for public_edge_penalty used
                     when separating a stale perceived state from a true
                     physical evaluation state. At most one alias may be set.
        """
        self.sf, self.mode = sf, mode
        self.E = sf["edges"]
        self.t0, self.cap = sf["t0"], sf["cap"]
        if public_edge_penalty is not None and perceived_edge_penalty is not None:
            raise ValueError("set only one of public_edge_penalty or perceived_edge_penalty")
        selected_penalty = public_edge_penalty if public_edge_penalty is not None else perceived_edge_penalty
        if selected_penalty is None:
            self.public_edge_penalty = np.zeros(self.E)
        else:
            penalty = np.asarray(selected_penalty, dtype=float)
            if penalty.shape != (self.E,) or np.any(penalty < 0):
                raise ValueError("public_edge_penalty must be a nonnegative length-E vector")
            self.public_edge_penalty = penalty
        self.routes = routes
        self.od_list = sf["od_list"]
        self.n_od = len(self.od_list)

        total_dem = sum(od["demand"] for od in self.od_list)
        self.total_dem = total_dem
        shares = np.asarray(shares, dtype=float)
        assert abs(shares.sum() - 1.0) < 1e-9, "shares must sum to 1"
        self.masses = shares * alpha * total_dem
        self.K = len(shares)
        self.alpha = alpha

        if beta_override is None:
            self.betas = 1.0 - np.exp(-lam_beta * self.masses)
        else:
            beta = np.asarray(beta_override, dtype=float)
            if beta.ndim == 0:
                beta = np.full(self.K, float(beta))
            if beta.shape != (self.K,):
                raise ValueError("beta_override must be scalar or length K")
            if np.any(beta <= 0):
                raise ValueError("beta_override must be strictly positive")
            self.betas = beta

        # A6: operator k serves m_k * d_w / N on OD w; HDVs serve (1-alpha) d_w.
        # Demand is conserved: sum_k m_k/N + (1-alpha) = alpha + (1-alpha) = 1.
        proportional_demands = np.array(
            [[self.masses[k] * od["demand"] / total_dem for k in range(self.K)]
             for od in self.od_list]
        )
        if dem_firm_override is None:
            self.dem_firm = proportional_demands
        else:
            demand_matrix = np.asarray(dem_firm_override, dtype=float)
            if demand_matrix.shape != (self.n_od, self.K):
                raise ValueError("dem_firm_override must have shape (n_od, K)")
            if np.any(demand_matrix < -1e-10):
                raise ValueError("dem_firm_override contains negative demand")
            target_rows = np.array([alpha * od["demand"] for od in self.od_list])
            if not np.allclose(demand_matrix.sum(axis=1), target_rows,
                               rtol=0, atol=1e-7):
                raise ValueError("dem_firm_override does not conserve OD demand")
            if not np.allclose(demand_matrix.sum(axis=0), self.masses,
                               rtol=0, atol=1e-6):
                raise ValueError("dem_firm_override does not match fleet masses")
            self.dem_firm = np.maximum(demand_matrix, 0.0)
        self.dem_hdv = np.array([(1.0 - alpha) * od["demand"] for od in self.od_list])
        self.kappa = kappa
        self.theta_bar = theta_bar        # (E,) or None

        self.coalitions = None
        self.objective_types = None
        self.coordination_gamma = float(coordination_gamma)
        self.operational_pooling = bool(operational_pooling)
        if management_slope_mode not in {"reference", "state"}:
            raise ValueError("management_slope_mode must be 'reference' or 'state'")
        self.management_slope_mode = management_slope_mode
        self.management_slope = None
        if coalitions is not None:
            if mode != "eq":
                raise ValueError("coalition objectives are available only in eq mode")
            blocks = [tuple(int(i) for i in block) for block in coalitions]
            flat = [i for block in blocks for i in block]
            if sorted(flat) != list(range(self.K)) or len(flat) != len(set(flat)):
                raise ValueError("coalitions must partition company indices 0..K-1")
            if not 0.0 <= self.coordination_gamma <= 1.0:
                raise ValueError("coordination_gamma must lie in [0, 1]")
            types = np.asarray(objective_types, dtype=float)
            if types.shape != (self.K,) or np.any(types <= 0):
                raise ValueError("objective_types must be a positive length-K vector")
            if management_slope_mode == "reference":
                if management_slope is None:
                    raise ValueError("management_slope is required in reference mode")
                slope = np.asarray(management_slope, dtype=float)
                if slope.shape != (self.E,) or np.any(slope < 0):
                    raise ValueError("management_slope must be a nonnegative length-E vector")
            else:
                slope = None if management_slope is None else np.asarray(management_slope, dtype=float)
                if slope is not None and (slope.shape != (self.E,) or np.any(slope < 0)):
                    raise ValueError("management_slope must be a nonnegative length-E vector")
            self.coalitions = tuple(blocks)
            self.objective_types = types
            self.management_slope = slope
            if self.operational_pooling:
                self.pool_demands = np.array([
                    [float(np.sum(self.dem_firm[oi, list(block)]))
                     for block in self.coalitions]
                    for oi in range(self.n_od)
                ])
        elif self.operational_pooling:
            raise ValueError("operational_pooling requires coalitions")

        self.rebuild()

    def coalition_management_edge_gradient(self, own_edge_flows, cp=None):
        """Return the convex quadratic management gradient for every company.

        Companies remain separate flow coordinates and preserve their own OD
        obligations.  A coalition is represented by symmetric PSD cross-member
        curvature.  At gamma=0 the result is independent of the partition;
        positive gamma activates only within-coalition cross terms.
        """
        if self.coalitions is None:
            raise RuntimeError("coalition governance is not configured")
        if self.management_slope_mode == "state":
            if cp is None:
                raise ValueError("state slope mode requires current BPR derivatives")
            slope = np.asarray(cp, dtype=float)
        else:
            slope = self.management_slope
        out = np.zeros_like(own_edge_flows)
        for block in self.coalitions:
            idx = np.asarray(block, dtype=int)
            governance = self.coalition_governance_matrix(block)
            out[idx] = (governance @ own_edge_flows[idx]) * slope
        return out

    def coalition_governance_matrix(self, block):
        """Dimensionless PSD member-coupling matrix for one coalition."""
        idx = np.asarray(block, dtype=int)
        types = self.objective_types[idx]
        root = np.sqrt(types)
        gamma = self.coordination_gamma
        return ((1.0 - gamma) * np.diag(types)
                + gamma * np.outer(root, root))

    def perceived_slope(self, cp, xk):
        """
        theta_hat^k_e = theta_bar_e + beta(m_k) * (c'_e - theta_bar_e).

        With theta_bar = 0 this collapses to beta(m_k)*c'_e, the specification
        used in the first version of the numerical section.  That corresponds to
        beta0 = 0, i.e. a prior mean slope of zero, which contradicts the
        maintained assumption theta_bar_e > 0; it is kept only so that the old
        numbers can be reproduced.
        """
        if self.theta_bar is None:
            return self.betas[:, None] * cp[None, :]
        tb = self.theta_bar[None, :]
        return tb + self.betas[:, None] * (cp[None, :] - tb)

    def rebuild(self):
        self.A, self.path_od = self.routes.incidence(self.E)
        self.P = self.A.shape[1]
        self.masks = [self.path_od == oi for oi in range(self.n_od)]
        self._build_projection_groups()

    def _build_projection_groups(self):
        """
        Group ODs by their number of routes so the projection can be done with
        one batched call per group instead of one call per (player, OD).
        Each group stores the column indices of its ODs as an (n_od_in_group,
        n_routes) index array.
        """
        sizes = {}
        for oi in range(self.n_od):
            cols = np.nonzero(self.masks[oi])[0]
            sizes.setdefault(len(cols), []).append((oi, cols))
        self._groups = []
        for n, items in sizes.items():
            ods = np.array([oi for oi, _ in items])
            cols = np.array([c for _, c in items])          # (G, n)
            self._groups.append((ods, cols))

    # -- flows ---------------------------------------------------------
    def edge_flow(self, fH, fF):
        return self.A @ (fH + fF.sum(axis=0))

    def initial(self, seed=0):
        rng = np.random.default_rng(seed)
        fH = np.zeros(self.P)
        fF = np.zeros((self.K, self.P))
        for oi in range(self.n_od):
            m = self.masks[oi]
            n = int(m.sum())
            if seed == 0:
                w = np.full(n, 1.0 / n)
            else:
                w = rng.random(n)
                w /= w.sum()
            fH[m] = self.dem_hdv[oi] * w
            for k in range(self.K):
                wk = w if seed == 0 else rng.dirichlet(np.ones(n))
                fF[k, m] = self.dem_firm[oi, k] * wk
        return fH, fF

    # -- the VI operator F --------------------------------------------
    def operator(self, fH, fF):
        """Return (grad_H, grad_F) : perceived marginal route costs."""
        x = self.edge_flow(fH, fF)
        ce = bpr_cost(x, self.t0, self.cap) + self.public_edge_penalty
        cp = bpr_deriv(x, self.t0, self.cap)
        if self.mode == "so":
            eff = ce + x * cp
            gH = self.A.T @ eff
            gF = np.tile(gH, (self.K, 1))
            return gH, gF, x, ce, cp
        pc = self.A.T @ ce
        gH = pc.copy()
        gF = np.empty((self.K, self.P))
        if self.mode == "ue":
            gF[:] = pc
        else:
            xk = fF @ self.A.T                      # (K, E) own edge flows
            if self.coalitions is None:
                th = self.perceived_slope(cp, xk)   # (K, E)
                marginal_management = th * xk
            else:
                marginal_management = self.coalition_management_edge_gradient(xk, cp)
            gF = pc[None, :] + (marginal_management @ self.A)
        return gH, gF, x, ce, cp

    def edge_weights(self, fH, fF):
        """
        Per-player edge weights whose route sums are exactly that player's route
        gradient.  These are the weights the reduced-cost oracle must use.

          HDV      : w_e = c_e
          fleet k  : w_e^k = c_e + beta_k * f_e^k * c'_e     (f_e^k known)
          mode 'so': w_e = c_e + x_e c'_e for everyone

        The naive bound w_e = c_e is NOT usable for fleets: their own-flow term
        is strictly positive on the routes they hold, so lambda_k always exceeds
        the plain shortest path and the test fires spuriously.
        """
        x = self.edge_flow(fH, fF)
        ce = bpr_cost(x, self.t0, self.cap) + self.public_edge_penalty
        cp = bpr_deriv(x, self.t0, self.cap)
        if self.mode == "so":
            w = ce + x * cp
            return w, np.tile(w, (self.K, 1))
        if self.mode == "ue":
            return ce, np.tile(ce, (self.K, 1))
        xk = fF @ self.A.T                          # (K, E)
        if self.coalitions is None:
            marginal_management = self.perceived_slope(cp, xk) * xk
        else:
            marginal_management = self.coalition_management_edge_gradient(xk, cp)
        wF = ce[None, :] + marginal_management
        return ce, wF

    def project(self, fH, fF):
        """Batched projection onto the product of demand-scaled simplices."""
        for ods, cols in self._groups:
            fH[cols] = project_simplex_batch(fH[cols], self.dem_hdv[ods])
            if not self.operational_pooling:
                for k in range(self.K):
                    fF[k][cols] = project_simplex_batch(
                        fF[k][cols], self.dem_firm[ods, k])
            else:
                for block_index, block in enumerate(self.coalitions):
                    idx = np.asarray(block, dtype=int)
                    # The pooled coalition simplex spans every member-route
                    # coordinate for this OD.  Assign each projected slice
                    # back explicitly because advanced indexing returns a copy.
                    values = np.stack([fF[k, cols] for k in idx], axis=1)
                    shape = values.shape
                    projected = project_simplex_batch(
                        values.reshape(shape[0], -1),
                        np.full(shape[0], self.pool_demands[ods, block_index]))
                    projected = projected.reshape(shape)
                    for member_pos, k in enumerate(idx):
                        fF[k, cols] = projected[:, member_pos, :]
        return fH, fF

    # -- normalised VI gap over ALL players ----------------------------
    def vi_gap(self, fH, fF):
        """
        gap = max over players p and OD w of
              [ max_{r used by p} PC_r^p  -  min_{r in set} PC_r^p ] / cost_scale
        This is the complementarity residual: zero iff every player puts flow
        only on its cheapest available routes.
        """
        gH, gF, x, ce, cp = self.operator(fH, fF)
        scale = max(float(np.max(self.A.T @ ce)), 1e-12)
        worst = 0.0
        for oi in range(self.n_od):
            m = self.masks[oi]
            if self.dem_hdv[oi] > 0:
                used = fH[m] > 1e-9 * self.dem_hdv[oi]
                if used.any():
                    g = gH[m]
                    worst = max(worst, float(g[used].max() - g.min()))
            if not self.operational_pooling:
                for k in range(self.K):
                    dk = self.dem_firm[oi, k]
                    if dk <= 0:
                        continue
                    used = fF[k, m] > 1e-9 * dk
                    if used.any():
                        g = gF[k][m]
                        worst = max(worst, float(g[used].max() - g.min()))
            else:
                for block_index, block in enumerate(self.coalitions):
                    idx = np.asarray(block, dtype=int)
                    demand = self.pool_demands[oi, block_index]
                    if demand <= 0:
                        continue
                    flat_flow = np.concatenate([fF[k, m] for k in idx])
                    flat_grad = np.concatenate([gF[k][m] for k in idx])
                    used = flat_flow > 1e-9 * demand
                    if used.any():
                        worst = max(worst, float(flat_grad[used].max()
                                                 - flat_grad.min()))
        return worst / scale

    # -- extragradient --------------------------------------------------
    def potential(self, fH, fF):
        """Convex potential value for the fixed-reference governance family."""
        if self.coalitions is None or self.management_slope_mode != "reference":
            raise RuntimeError("potential is defined here only for reference-mode coalitions")
        x = self.edge_flow(fH, fF)
        physical = np.sum(
            self.t0 * x + self.t0 * BPR_B * x**(BPR_ETA + 1.0)
            / ((BPR_ETA + 1.0) * self.cap**BPR_ETA)
            + self.public_edge_penalty * x
        )
        own = fF @ self.A.T
        management = 0.0
        for block in self.coalitions:
            idx = np.asarray(block, dtype=int)
            governance = self.coalition_governance_matrix(block)
            management += 0.5 * np.sum(
                self.management_slope
                * np.einsum("ie,ij,je->e", own[idx], governance, own[idx])
            )
        return float(physical + management)

    def solve_potential(self, seed=0, tol=1e-8, max_iter=20000,
                        verbose=False, warm=None, report_every=2000):
        """Projected descent for the convex fixed-reference potential."""
        fH, fF = (warm[0].copy(), warm[1].copy()) if warm is not None else self.initial(seed)
        fH, fF = self.project(fH, fF)
        gH, gF, *_ = self.operator(fH, fF)
        step = 0.1 * self.total_dem / max(np.abs(gH).max(), np.abs(gF).max(), 1e-12)
        current = self.potential(fH, fF)
        best = np.inf
        for it in range(max_iter):
            accepted = False
            for _ in range(50):
                yH = fH - step * gH
                yF = fF - step * gF
                yH, yF = self.project(yH, yF)
                trial = self.potential(yH, yF)
                if trial <= current + 1e-12 * max(1.0, abs(current)):
                    accepted = True
                    break
                step *= 0.5
            if not accepted:
                break
            fH, fF = yH, yF
            current = trial
            gH, gF, *_ = self.operator(fH, fF)
            if it % 25 == 0:
                gap = self.vi_gap(fH, fF)
                best = min(best, gap)
                if verbose and it % report_every == 0:
                    print(f"      it={it:6d} gap={gap:.3e} step={step:.2e}", flush=True)
                if gap < tol:
                    break
            step = min(step * 1.05, 0.1 * self.total_dem / max(np.abs(gH).max(), np.abs(gF).max(), 1e-12))
        gap = self.vi_gap(fH, fF)
        x = self.edge_flow(fH, fF)
        J = float(np.dot(x, bpr_cost(x, self.t0, self.cap)))
        return dict(fH=fH, fF=fF, x=x, J=J, gap=gap, iters=it + 1,
                    gamma=step, best_gap=best)

    def solve(self, seed=0, tol=1e-8, max_iter=20000, gamma0=None,
              verbose=False, warm=None, report_every=2000, method="extragradient"):
        """
        Korpelevich extragradient with an adaptive step.  For a monotone
        Lipschitz operator this converges to a solution of the VI; for the
        BPR case monotonicity is not certified, so we report the achieved gap
        and let the caller judge.

        Step control: the trial step must satisfy
            gamma * ||F(y) - F(z)||  <=  delta * ||y - z||        (delta < 1)
        which is the standard extragradient condition.  We shrink by 0.5 on
        failure and grow by 1.3 on a first-try success, so the step recovers
        quickly after a transient; a slow growth factor was the reason an
        earlier version stalled on the eight-operator case.
        """
        if method == "potential":
            if self.coalitions is None or self.management_slope_mode != "reference":
                raise ValueError("potential method requires reference-mode coalition governance")
            return self.solve_potential(seed=seed, tol=tol, max_iter=max_iter,
                                        verbose=verbose, warm=warm,
                                        report_every=report_every)
        if warm is not None:
            fH, fF = warm[0].copy(), warm[1].copy()
        else:
            fH, fF = self.initial(seed)
        fH, fF = self.project(fH, fF)
        gH, gF, *_ = self.operator(fH, fF)
        nrm = max(np.abs(gH).max(), np.abs(gF).max(), 1e-12)
        gamma = gamma0 if gamma0 else 0.1 * self.total_dem / nrm

        best = np.inf
        for it in range(max_iter):
            tries = 0
            while True:
                yH = fH - gamma * gH
                yF = fF - gamma * gF
                yH, yF = self.project(yH.copy(), yF.copy())
                gyH, gyF, *_ = self.operator(yH, yF)
                dz = np.sqrt(np.sum((yH - fH) ** 2) + np.sum((yF - fF) ** 2))
                dg = np.sqrt(np.sum((gyH - gH) ** 2) + np.sum((gyF - gF) ** 2))
                if dz < 1e-15 or gamma * dg <= 0.9 * dz or tries >= 40:
                    break
                gamma *= 0.5
                tries += 1
            if tries == 0:
                gamma *= 1.3                      # recover the step quickly
            zH = fH - gamma * gyH
            zF = fF - gamma * gyF
            zH, zF = self.project(zH, zF)
            fH, fF = zH, zF
            gH, gF, *_ = self.operator(fH, fF)

            if it % 25 == 0:
                g = self.vi_gap(fH, fF)
                best = min(best, g)
                if verbose and it % report_every == 0:
                    print(f"      it={it:6d} gap={g:.3e} gamma={gamma:.2e}",
                          flush=True)
                if g < tol:
                    break
        gap = self.vi_gap(fH, fF)
        x = self.edge_flow(fH, fF)
        J = float(np.dot(x, bpr_cost(x, self.t0, self.cap)))
        return dict(fH=fH, fF=fF, x=x, J=J, gap=gap, iters=it + 1,
                    gamma=gamma, best_gap=best)


# ----------------------------------------------------------------------
# column generation with reduced-cost certificate
# ----------------------------------------------------------------------
def build_graph(sf_links):
    G = nx.DiGraph()
    edge_idx = {}
    for i, (u, v, cap, t0) in enumerate(sf_links):
        edge_idx[(u, v)] = i
        G.add_edge(u, v, eidx=i)
    return G, edge_idx


def certificate(game, sol, G, edge_idx, tol_abs):
    """
    Exact reduced-cost test for every player.

    Player p's route gradient is the route sum of a *player-specific* edge
    weight w_e^p (see ``SFGame.edge_weights``).  A route omitted from the set is
    profitable for p exactly when

        min_r sum_{e in r} w_e^p   <   lambda_{p,w} - tol,

    and the left-hand side is delivered by one Dijkstra run per (player,
    origin) under w^p.  This is exact, not a bound: the gradient is evaluated at
    the current profile, which is what the KKT test requires.

    Return (violations, worst_slack) with violations = [(oi, edge_path, slack)].
    """
    fH, fF = sol["fH"], sol["fF"]
    gH, gF, x, ce, cp = game.operator(fH, fF)
    wH, wF = game.edge_weights(fH, fF)

    # group ODs by origin so each player needs one Dijkstra per origin
    by_origin = {}
    for oi, od in enumerate(game.od_list):
        by_origin.setdefault(od["o"], []).append(oi)

    # Players: -1 denotes HDVs.  Under pooling, a coalition is one joint
    # player even though its member route flows remain separate coordinates.
    players = ([-1] if game.dem_hdv.max() > 0 else [])
    if game.operational_pooling:
        players += [tuple(block) for block in game.coalitions]
    else:
        players += list(range(game.K))
    coalition_index = {
        tuple(block): index for index, block in enumerate(game.coalitions or ())
    }

    violations = {}
    worst = -np.inf
    for pid in players:
        members = [-1] if pid == -1 else (
            list(pid) if game.operational_pooling else [pid])
        # Store a shortest-path oracle for each member because pooled members
        # may have different marginal governance costs on the same edge.
        distances, paths_by_member = {}, {}
        for member in members:
            w = wH if member == -1 else wF[member]
            for u, v in list(G.edges()):
                G[u][v]["w"] = float(w[edge_idx[(u, v)]])
            for origin in by_origin:
                dist, paths = nx.single_source_dijkstra(G, origin, weight="w")
                distances[(member, origin)] = dist
                paths_by_member[(member, origin)] = paths
        for origin, ods in by_origin.items():
            for oi in ods:
                od = game.od_list[oi]
                if od["d"] not in distances[(members[0], origin)]:
                    continue
                m = game.masks[oi]
                if pid == -1:
                    dem = game.dem_hdv[oi]
                    lam = float(gH[m].min())
                elif game.operational_pooling:
                    block_index = coalition_index[tuple(pid)]
                    dem = float(game.pool_demands[oi, block_index])
                    lam = float(np.min(np.concatenate([gF[k][m] for k in members])))
                else:
                    dem = game.dem_firm[oi, pid]
                    lam = float(gF[pid][m].min())
                if dem <= 0:
                    continue
                for member in members:
                    dist = distances[(member, origin)]
                    slack = lam - float(dist[od["d"]])
                    worst = max(worst, slack)
                    if slack > tol_abs:
                        nodes = paths_by_member[(member, origin)][od["d"]]
                        epath = tuple(edge_idx[(nodes[i], nodes[i + 1])]
                                      for i in range(len(nodes) - 1))
                        prev = violations.get((oi, epath), -np.inf)
                        violations[(oi, epath)] = max(prev, slack)
    viol = [(oi, list(ep), s) for (oi, ep), s in violations.items()]
    return viol, worst


def solve_with_column_generation(sf, alpha, shares, lam_beta, mode="eq",
                                 seed=0, tol=1e-8, cg_rounds=25,
                                 cg_tol_rel=1e-6, verbose=True):
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    G, edge_idx = build_graph(links)
    routes = RouteSet(sf["od_list"])
    game = SFGame(sf, alpha, shares, lam_beta, routes, mode=mode)
    warm = None

    for rnd in range(cg_rounds):
        sol = game.solve(seed=seed, tol=tol, warm=warm, verbose=verbose)
        wH, wF = game.edge_weights(sol["fH"], sol["fF"])
        scale = float(np.max(game.A.T @ wH))
        viol, worst = certificate(game, sol, G, edge_idx, cg_tol_rel * scale)
        if verbose:
            print(f"  [CG {rnd}] paths={routes.total()} gap={sol['gap']:.2e} "
                  f"worst_reduced_cost={worst:+.3e} (scale={scale:.1f}) "
                  f"violations={len(viol)}", flush=True)
        if not viol:
            sol["cg_rounds"] = rnd + 1
            sol["n_paths"] = routes.total()
            sol["max_reduced_cost"] = worst
            sol["certified"] = True
            return sol, game
        old_keys = routes.keys()
        added = sum(routes.add(oi, ep) for oi, ep, _ in viol)
        if added == 0:
            sol["cg_rounds"] = rnd + 1
            sol["n_paths"] = routes.total()
            sol["max_reduced_cost"] = worst
            sol["certified"] = False
            return sol, game
        # warm start: carry old route flows into the enlarged index set
        pos = {k: i for i, k in enumerate(old_keys)}
        game.rebuild()
        new_keys = routes.keys()
        fH = np.zeros(game.P)
        fF = np.zeros((game.K, game.P))
        for j, k in enumerate(new_keys):
            i = pos.get(k)
            if i is not None:
                fH[j] = sol["fH"][i]
                fF[:, j] = sol["fF"][:, i]
        warm = (fH, fF)

    sol["certified"] = False
    sol["cg_rounds"] = cg_rounds
    sol["n_paths"] = routes.total()
    sol["max_reduced_cost"] = worst
    return sol, game


if __name__ == "__main__":
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    print(f"Sioux Falls: {sf['edges']} links, {len(sf['od_list'])} OD pairs, "
          f"total demand {sum(o['demand'] for o in sf['od_list']):.0f}")
    print("\nUser equilibrium (column generation):")
    sol_ue, _ = solve_with_column_generation(sf, 0.9, [1.0], 3e-4, mode="ue")
    print(f"  J_UE = {sol_ue['J']:.1f}  certified={sol_ue['certified']}  "
          f"paths={sol_ue['n_paths']}")
