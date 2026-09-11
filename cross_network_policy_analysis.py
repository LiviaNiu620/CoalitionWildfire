"""Evaluate greedy merge and fixed policies on certified cross-network states."""
from __future__ import annotations
import itertools, json
from pathlib import Path

RAW = Path("cross_network_regime_results.json")
SUMMARY = Path("cross_network_regime_summary.json")

def canon(p): return tuple(sorted(tuple(sorted(b)) for b in p))
def coarsenings(p):
    blocks = [tuple(b) for b in p]
    for i, j in itertools.combinations(range(len(blocks)), 2):
        yield canon([blocks[k] for k in range(len(blocks)) if k not in (i,j)] + [tuple(sorted(blocks[i]+blocks[j]))])

def main():
    raw = json.loads(RAW.read_text()); states = []; policy = {k: [] for k in ("greedy_merge", "two_step_merge", "always_grand", "always_singleton")}
    for network, bench in raw["benchmarks"].items():
        for alpha in raw["protocol"]["alphas"]:
            for gamma in raw["protocol"]["gammas"]:
                rows = [r for r in raw["profiles"] if r["network"] == network and r["alpha"] == alpha and r["gamma"] == gamma]
                by = {canon(r["partition"]): r for r in rows}; best = min(rows, key=lambda r:r["J"]); eps = 1e-8 * bench["J_UE"]
                cur = canon([(i,) for i in range(4)]); path = [by[cur]["partition_id"]]
                while len(cur) > 1:
                    choice = min((by[p] for p in coarsenings(cur)), key=lambda r:r["J"])
                    if choice["J"] >= by[cur]["J"] - eps: break
                    cur = canon(choice["partition"]); path.append(choice["partition_id"])
                greedy = by[cur]; grand = next(r for r in rows if r["n_coalitions"]==1); singleton = next(r for r in rows if r["n_coalitions"]==4)
                cur2 = canon([(i,) for i in range(4)]); path2 = [by[cur2]["partition_id"]]
                while len(cur2) > 1:
                    candidates=[]
                    for first in coarsenings(cur2):
                        candidates.append((by[first]["J"],first))
                        for second in coarsenings(first): candidates.append((by[second]["J"],first))
                    future_j,nxt=min(candidates,key=lambda item:item[0])
                    if future_j >= by[cur2]["J"]-eps: break
                    cur2=nxt;path2.append(by[cur2]["partition_id"])
                lookahead=by[cur2]
                def reg(r): return {"recovery_regret_pp":100*(best["recovery"]-r["recovery"]), "delta_J_pct_UE":100*(r["J"]-best["J"])/bench["J_UE"]}
                for key,row in (("greedy_merge",greedy),("two_step_merge",lookahead),("always_grand",grand),("always_singleton",singleton)): policy[key].append(reg(row))
                states.append({"network":network,"alpha":alpha,"gamma":gamma,"best_partition_id":best["partition_id"],"best_K":best["n_coalitions"],"best_recovery":best["recovery"],"gain_vs_singleton_pp":100*(best["recovery"]-singleton["recovery"]),"grand_regret_pp":100*(best["recovery"]-grand["recovery"]),"greedy_path":path,"greedy_partition_id":greedy["partition_id"],"greedy_regret_pp":reg(greedy)["recovery_regret_pp"],"two_step_path":path2,"two_step_partition_id":lookahead["partition_id"],"two_step_regret_pp":reg(lookahead)["recovery_regret_pp"]})
    stats={}
    for key, vals in policy.items():
        stats[key]={"mean_recovery_regret_pp":sum(v["recovery_regret_pp"] for v in vals)/len(vals),"max_recovery_regret_pp":max(v["recovery_regret_pp"] for v in vals),"exact_state_count":sum(abs(v["recovery_regret_pp"])<=1e-8 for v in vals),"states":len(vals)}
    output={"states":states,"policy_regret_summary":stats,"certified_rows":sum(r["certified"] for r in raw["profiles"]),"rows":len(raw["profiles"]),"max_VI_gap":max(r["VI_gap"] for r in raw["profiles"]),"max_relative_omitted_route_slack":max(r["omitted_route_slack"]/r["route_cost_scale"] for r in raw["profiles"])}
    SUMMARY.write_text(json.dumps(output,indent=2)+"\n"); print(json.dumps(output,indent=2))

if __name__ == "__main__": main()
