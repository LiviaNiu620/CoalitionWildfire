"""Deterministic wildfire-state snapshots for the Sioux Falls adapter.

The scenarios are transport-state counterfactuals: they do not model fire
spread or household evacuation. Closed edges are removed from the route set
and shortest-path oracle; open edges may receive capacity/free-flow changes.
"""
from __future__ import annotations

import copy
from pathlib import Path

import networkx as nx
import numpy as np

from sioux_falls_loader import build_sioux_falls, parse_net


SCENARIO_SCHEMA = "wildfire_snapshot_v1"


def _edge_frequency(sf, limit=20):
    freq = np.zeros(sf["edges"], dtype=int)
    for od in sf["od_list"][:limit]:
        for path in od["paths"]:
            for edge in path:
                freq[edge] += 1
    return freq


def scenario_catalog(sf, od_limit=20):
    """Return deterministic smoke/pilot scenarios based on route incidence."""
    freq = _edge_frequency(sf, limit=od_limit)
    ranked = np.argsort(-freq)
    focus = [int(e) for e in ranked if freq[e] > 0]
    if not focus:
        raise ValueError("no route-incidence edges available for scenarios")
    # Prefer an edge with an alternate route for every selected OD.
    closure = next(
        edge for edge in focus
        if all(any(edge not in path for path in od["paths"]) for od in sf["od_list"])
    )
    degraded = focus[: min(5, len(focus))]
    return [
        {
            "schema": SCENARIO_SCHEMA,
            "scenario_id": "no_hazard",
            "closed_edges": [],
            "capacity_multiplier": {},
            "free_flow_multiplier": {},
            "risk_score": {},
            "critical_edges": degraded[: min(3, len(degraded))],
            "public_information_delay": 0,
        },
        {
            "schema": SCENARIO_SCHEMA,
            "scenario_id": "mild_capacity_degradation",
            "closed_edges": [],
            "capacity_multiplier": {str(e): 0.75 for e in degraded},
            "free_flow_multiplier": {},
            "risk_score": {str(e): 0.25 for e in degraded},
            "critical_edges": degraded,
            "public_information_delay": 0,
        },
        {
            "schema": SCENARIO_SCHEMA,
            "scenario_id": "critical_edge_closure",
            "closed_edges": [closure],
            "capacity_multiplier": {},
            "free_flow_multiplier": {},
            "risk_score": {str(e): 1.0 for e in degraded},
            "critical_edges": degraded,
            "public_information_delay": 0,
        },
        {
            "schema": SCENARIO_SCHEMA,
            "scenario_id": "multiple_corridor_degradation",
            "closed_edges": [],
            "capacity_multiplier": {str(e): 0.55 for e in degraded},
            "free_flow_multiplier": {str(e): 1.15 for e in degraded},
            "risk_score": {str(e): 0.75 for e in degraded},
            "critical_edges": degraded,
            "public_information_delay": 1,
        },
    ]


def apply_snapshot(sf, scenario):
    """Apply one scenario and remove paths containing closed edges."""
    if scenario.get("schema") != SCENARIO_SCHEMA:
        raise ValueError("unsupported wildfire scenario schema")
    closed = {int(e) for e in scenario.get("closed_edges", [])}
    out = copy.deepcopy(sf)
    out["t0"] = np.asarray(sf["t0"], dtype=float).copy()
    out["cap"] = np.asarray(sf["cap"], dtype=float).copy()
    for edge, multiplier in scenario.get("capacity_multiplier", {}).items():
        out["cap"][int(edge)] *= float(multiplier)
    for edge, multiplier in scenario.get("free_flow_multiplier", {}).items():
        out["t0"][int(edge)] *= float(multiplier)
    for od in out["od_list"]:
        od["paths"] = [list(path) for path in od["paths"] if not (set(path) & closed)]
        if not od["paths"]:
            raise ValueError(
                f"scenario {scenario['scenario_id']} disconnects OD {od['o']}->{od['d']}"
            )
    out["hazard"] = {
        "schema": SCENARIO_SCHEMA,
        "scenario_id": scenario["scenario_id"],
        "closed_edges": sorted(closed),
        "risk_score": {str(k): float(v) for k, v in scenario.get("risk_score", {}).items()},
        "critical_edges": [int(e) for e in scenario.get("critical_edges", [])],
        "public_information_delay": int(scenario.get("public_information_delay", 0)),
    }
    return out


def hazard_oracle(scenario, net_path=None):
    """Build a shortest-path oracle with closed roads removed."""
    links = parse_net(net_path) if net_path else parse_net(
        str(Path(__file__).with_name("sioux_falls_data") / "SiouxFalls_net.tntp")
    )
    closed = {int(e) for e in scenario.get("closed_edges", [])}
    graph = nx.DiGraph()
    edge_idx = {}
    for index, (u, v, cap, t0) in enumerate(links):
        if index in closed:
            continue
        edge_idx[(u, v)] = index
        graph.add_edge(u, v, eidx=index)
    return graph, edge_idx


def build_snapshot(od_limit=20, demand_scale=0.1):
    """Return base Sioux Falls data and its deterministic scenario catalog."""
    base = build_sioux_falls(k_paths=3, top_od=od_limit, demand_scale=demand_scale)
    return base, scenario_catalog(base, od_limit=od_limit)


__all__ = ["SCENARIO_SCHEMA", "scenario_catalog", "apply_snapshot", "hazard_oracle", "build_snapshot"]
