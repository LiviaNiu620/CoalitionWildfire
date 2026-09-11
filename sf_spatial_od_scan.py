"""Fixed-HHI Sioux Falls experiment with heterogeneous spatial OD portfolios."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import networkx as nx
import numpy as np

from sf_revision_common import (
    K,
    LAM_TIMES_M,
    SE,
    build_sioux_falls,
    network_oracle,
    profile_record,
    solve_benchmarks,
    solve_profile,
    write_json,
)
from sioux_falls_loader import DATA_DIR, parse_net


ALPHAS = (0.5, 0.9)
FAMILIES = ("origin", "destination", "joint_od", "corridor")
TEMPERATURES = (0.35, 0.75)
SEEDS = (20260801, 20260802, 20260803, 20260804, 20260805)
RAW_PATH = Path("sf_spatial_od_results.json")
SUMMARY_PATH = Path("sf_spatial_od_summary.json")


def chunk_path(index):
    return Path(f"sf_spatial_od_chunk_{index:02d}.json")


def standardize(features):
    features = np.asarray(features, dtype=float)
    scale = features.std(axis=0)
    keep = scale > 1e-12
    if not np.any(keep):
        raise RuntimeError("spatial feature family has no variation")
    centered = features[:, keep] - features[:, keep].mean(axis=0)
    return centered / scale[keep]


def spatial_features(sf):
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    graph = nx.DiGraph()
    for u, v, _, t0 in links:
        graph.add_edge(u, v, weight=t0)
    nodes = sorted(graph.nodes())
    node_pos = {node: i for i, node in enumerate(nodes)}
    forward = dict(nx.all_pairs_dijkstra_path_length(graph, weight="weight"))
    reverse = dict(nx.all_pairs_dijkstra_path_length(graph.reverse(), weight="weight"))

    origin_rows = []
    destination_rows = []
    corridor_rows = []
    for od in sf["od_list"]:
        origin_rows.append([forward[od["o"]][node] for node in nodes])
        destination_rows.append([reverse[od["d"]][node] for node in nodes])
        incidence = np.zeros(sf["edges"])
        incidence[od["paths"][0]] = 1.0
        corridor_rows.append(incidence)
    origin = standardize(origin_rows)
    destination = standardize(destination_rows)
    return {
        "origin": origin,
        "destination": destination,
        "joint_od": standardize(np.hstack([origin, destination])),
        "corridor": standardize(corridor_rows),
    }


def farthest_anchors(features, seed):
    rng = np.random.default_rng(seed)
    selected = [int(rng.integers(0, len(features)))]
    while len(selected) < K:
        delta = features[:, None, :] - features[np.array(selected)][None, :, :]
        min_distance = np.min(np.sum(delta * delta, axis=2), axis=1)
        min_distance[np.array(selected)] = -1.0
        candidates = np.flatnonzero(np.isclose(min_distance, min_distance.max()))
        selected.append(int(rng.choice(candidates)))
    return selected


def affinity_matrix(features, temperature, seed):
    anchors = farthest_anchors(features, seed)
    delta = features[:, None, :] - features[np.array(anchors)][None, :, :]
    squared = np.sum(delta * delta, axis=2)
    positive = np.sqrt(squared[squared > 1e-12])
    scale = float(np.median(positive))
    bandwidth = max(temperature * scale, 1e-8)
    affinity = np.exp(-squared / (2.0 * bandwidth * bandwidth))
    return affinity + 1e-8, anchors


def balance_portfolio(affinity, row_targets):
    """IPF with exact OD and equal-fleet margins to numerical tolerance."""
    row_targets = np.asarray(row_targets, dtype=float)
    col_targets = np.full(K, row_targets.sum() / K)
    matrix = np.asarray(affinity, dtype=float).copy()
    matrix *= row_targets[:, None] / matrix.sum(axis=1)[:, None]
    for _ in range(20_000):
        matrix *= col_targets[None, :] / matrix.sum(axis=0)[None, :]
        matrix *= row_targets[:, None] / matrix.sum(axis=1)[:, None]
        row_error = np.max(np.abs(matrix.sum(axis=1) - row_targets))
        col_error = np.max(np.abs(matrix.sum(axis=0) - col_targets))
        if max(row_error, col_error) < 1e-8:
            break
    else:
        raise RuntimeError("IPF did not converge")
    if not np.allclose(matrix.sum(axis=0), col_targets, rtol=0, atol=1e-6):
        raise RuntimeError("IPF fleet margins failed")
    return matrix


def portfolio_design(sf, alpha):
    demand = np.array([od["demand"] for od in sf["od_list"]])
    row_targets = alpha * demand
    proportional = np.tile(row_targets[:, None] / K, (1, K))
    designs = [{
        "id": "proportional",
        "family": "proportional",
        "temperature": None,
        "seed": None,
        "anchors": None,
        "demands": proportional,
    }]
    features = spatial_features(sf)
    for family in FAMILIES:
        for temperature in TEMPERATURES:
            for seed in SEEDS:
                affinity, anchors = affinity_matrix(features[family], temperature, seed)
                matrix = balance_portfolio(affinity, row_targets)
                designs.append({
                    "id": f"{family}-t{temperature:.2f}-s{seed}",
                    "family": family,
                    "temperature": temperature,
                    "seed": seed,
                    "anchors": anchors,
                    "demands": matrix,
                })
    if len(designs) != 41:
        raise RuntimeError(f"expected 41 portfolios, found {len(designs)}")
    return designs


def portfolio_metadata(design, alpha, total_demand):
    matrix = design["demands"]
    row_totals = matrix.sum(axis=1)
    shares = matrix / row_totals[:, None]
    entropy = -np.sum(shares * np.log(np.maximum(shares, 1e-300)), axis=1) / np.log(K)
    weighted_entropy = float(np.dot(row_totals, entropy) / row_totals.sum())
    proportional = row_totals[:, None] / K
    reassigned_fraction = float(np.sum(np.abs(matrix - proportional)) / (2 * alpha * total_demand))
    masses = matrix.sum(axis=0)
    fleet_shares = masses / masses.sum()
    return {
        "id": design["id"],
        "family": design["family"],
        "temperature": design["temperature"],
        "seed": design["seed"],
        "anchor_od_indices": design["anchors"],
        "alpha": alpha,
        "shares": fleet_shares.tolist(),
        "HHI": float(np.dot(fleet_shares, fleet_shares)),
        "weighted_OD_entropy": weighted_entropy,
        "reassigned_demand_fraction": reassigned_fraction,
        "max_fleet_mass_error": float(np.max(np.abs(masses - masses.mean()))),
    }


def distribution_summary(values):
    values = np.asarray(values, dtype=float)
    return {
        "min": float(values.min()),
        "q25": float(np.quantile(values, 0.25)),
        "median": float(np.median(values)),
        "q75": float(np.quantile(values, 0.75)),
        "max": float(values.max()),
        "mean": float(values.mean()),
        "sd": float(np.std(values, ddof=1)),
    }


def summarize(raw):
    by_alpha = {}
    for alpha in ALPHAS:
        rows = [row for row in raw["profiles"] if row["alpha"] == alpha]
        certified = [row for row in rows if row["certified"]]
        baseline = next(row for row in certified if row["family"] == "proportional")
        spatial = [row for row in certified if row["family"] != "proportional"]
        by_family = {}
        for family in FAMILIES:
            family_rows = [row for row in spatial if row["family"] == family]
            by_family[family] = {
                "n": len(family_rows),
                "J": distribution_summary([row["J"] for row in family_rows]),
                "recovery": distribution_summary([row["recovery"] for row in family_rows]),
            }
        by_alpha[f"{alpha:g}"] = {
            "n": len(rows),
            "certified_n": len(certified),
            "HHI_range": [
                float(min(row["HHI"] for row in certified)),
                float(max(row["HHI"] for row in certified)),
            ],
            "proportional_baseline": {"J": baseline["J"], "recovery": baseline["recovery"]},
            "spatial_J": distribution_summary([row["J"] for row in spatial]),
            "spatial_recovery": distribution_summary([row["recovery"] for row in spatial]),
            "max_abs_delta_J_from_proportional": float(max(
                abs(row["J"] - baseline["J"]) for row in spatial
            )),
            "max_abs_delta_J_pct_UE": float(100 * max(
                abs(row["J"] - baseline["J"]) for row in spatial
            ) / raw["J_UE"]),
            "profile_psd_pass_n": int(sum(row["profile_psd_certified"] for row in certified)),
            "edge_psd_pass_rate": float(
                sum(row["positive_slope_edges"] * row["edge_psd_pass_fraction"]
                    for row in certified)
                / sum(row["positive_slope_edges"] for row in certified)
            ),
            "max_exact_violation_ratio": float(max(
                row["max_exact_violation_ratio"] for row in certified
            )),
            "profile_coarse_pass_n": int(sum(
                row["profile_coarse_certified"] for row in certified
            )),
            "edge_coarse_pass_rate": float(
                sum(row["positive_slope_edges"] * row["coarse_bound_pass_fraction"]
                    for row in certified)
                / sum(row["positive_slope_edges"] for row in certified)
            ),
            "max_coarse_violation_ratio": float(max(
                row["max_coarse_violation_ratio"] for row in certified
            )),
            "by_family": by_family,
        }
    return {"protocol": raw["protocol"], "rows": len(raw["profiles"]), "by_alpha": by_alpha}


def merge_chunks(count):
    chunks = [json.loads(chunk_path(i).read_text()) for i in range(count)]
    rows = [row for chunk in chunks for row in chunk["profiles"]]
    keys = [(row["id"], row["alpha"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate spatial profile/alpha keys")
    expected = 41 * len(ALPHAS)
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} rows, found {len(rows)}")
    first = chunks[0]
    for chunk in chunks[1:]:
        if not np.isclose(chunk["J_UE"], first["J_UE"], rtol=0, atol=1e-6):
            raise RuntimeError("inconsistent J_UE across spatial chunks")
        if not np.isclose(chunk["J_SO"], first["J_SO"], rtol=0, atol=1e-6):
            raise RuntimeError("inconsistent J_SO across spatial chunks")
    if not all(np.isclose(row["HHI"], 0.25, rtol=0, atol=1e-10) for row in rows):
        raise RuntimeError("fixed-HHI design drifted from 0.25")
    if max(row["max_fleet_mass_error"] for row in rows) > 1e-6:
        raise RuntimeError("fixed-HHI fleet-mass margins failed")
    protocol = dict(first["protocol"])
    protocol.pop("chunk_index", None)
    protocol.pop("chunks", None)
    protocol["parallel_chunks"] = count
    raw = {
        "protocol": protocol,
        "J_UE": first["J_UE"],
        "J_SO": first["J_SO"],
        "profiles": sorted(rows, key=lambda row: (row["alpha"], row["id"])),
    }
    write_json(RAW_PATH, raw)
    summary = summarize(raw)
    write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2), flush=True)
    if any(value["certified_n"] != value["n"] for value in summary["by_alpha"].values()):
        raise SystemExit("one or more spatial profiles failed certification")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunks", type=int, default=1)
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--merge-chunks", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    if args.merge_chunks:
        merge_chunks(args.chunks)
        return
    if not 0 <= args.chunk_index < args.chunks:
        parser.error("invalid chunk index")

    output = RAW_PATH if args.chunks == 1 else chunk_path(args.chunk_index)
    sf = build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total = sum(item["demand"] for item in sf["od_list"])
    ue, _, so = solve_benchmarks(sf)
    oracle = network_oracle()
    routes = SE.RouteSet(sf["od_list"])
    designs = {alpha: portfolio_design(sf, alpha) for alpha in ALPHAS}
    indexed = [
        (alpha, design) for alpha in ALPHAS for design in designs[alpha]
    ]
    selected = [item for i, item in enumerate(indexed) if i % args.chunks == args.chunk_index]

    prior = {}
    if args.resume and output.exists():
        old = json.loads(output.read_text())
        prior = {(row["id"], row["alpha"]): row for row in old.get("profiles", [])}

    raw = {
        "protocol": {
            "alphas": list(ALPHAS),
            "K": K,
            "HHI": 0.25,
            "beta_override": 1.0,
            "families": list(FAMILIES),
            "temperatures": list(TEMPERATURES),
            "seeds": list(SEEDS),
            "portfolios_per_alpha": 41,
            "margin_method": "strictly positive affinity plus iterative proportional fitting",
            "VI_tolerance": 1e-6,
            "omitted_route_relative_tolerance": 1e-6,
            "chunks": args.chunks,
            "chunk_index": args.chunk_index,
        },
        "J_UE": float(ue["J"]),
        "J_SO": float(so["J"]),
        "profiles": [],
    }

    print("Building spatial-portfolio route basis", flush=True)
    anchor_ids = {"proportional", "origin-t0.35-s20260801", "destination-t0.35-s20260802",
                  "joint_od-t0.35-s20260803", "corridor-t0.35-s20260804"}
    for alpha in ALPHAS:
        lam_beta = LAM_TIMES_M / (alpha * total)
        for design in designs[alpha]:
            if design["id"] not in anchor_ids:
                continue
            sol, _ = solve_profile(
                sf, routes, oracle, alpha, [0.25] * K, lam_beta,
                beta_override=1.0, dem_firm_override=design["demands"],
                verbose=args.verbose,
            )
            if not sol["certified"]:
                raise RuntimeError(
                    "spatial route-basis anchor failed: "
                    f"alpha={alpha:g}, id={design['id']}, "
                    f"gap={sol['gap']:.3e}, "
                    f"relative_omitted={sol['max_reduced_cost'] / sol['route_cost_scale']:.3e}, "
                    f"cg_rounds={sol['cg_rounds']}"
                )
    print(f"Shared basis contains {routes.total()} paths", flush=True)

    warm = {alpha: None for alpha in ALPHAS}
    for index, (alpha, design) in enumerate(selected, 1):
        key = (design["id"], alpha)
        if key in prior:
            raw["profiles"].append(prior[key])
            continue
        lam_beta = LAM_TIMES_M / (alpha * total)
        sol, game = solve_profile(
            sf, routes, oracle, alpha, [0.25] * K, lam_beta,
            beta_override=1.0, dem_firm_override=design["demands"],
            warm=warm[alpha], verbose=args.verbose,
        )
        metadata = portfolio_metadata(design, alpha, total)
        record = profile_record(metadata, sol, game, ue["J"], so["J"], True)
        raw["profiles"].append(record)
        warm[alpha] = (sol["fH"], sol["fF"])
        write_json(output, raw)
        print(
            f"[{index:03d}/{len(selected)}] alpha={alpha:.1f} {design['id']:<34} "
            f"rho={record['recovery']:.4f} entropy={record['weighted_OD_entropy']:.3f} "
            f"cert={record['certified']}",
            flush=True,
        )

    raw["profiles"].sort(key=lambda row: (row["alpha"], row["id"]))
    write_json(output, raw)
    if args.chunks == 1:
        summary = summarize(raw)
        write_json(SUMMARY_PATH, summary)
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
