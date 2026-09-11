"""
Sioux Falls 真实数据加载 + networkx 生成 k-最短路径集。
数据: sioux_falls_data/SiouxFalls_net.tntp (76 links), SiouxFalls_trips.tntp (OD矩阵)。
输出: 可直接喂给 BPRNetwork 的 (edges, t0, cap, ods) + 各OD的k条候选路径。
"""
import numpy as np
import networkx as nx
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "sioux_falls_data")


def parse_net(path):
    """解析 SiouxFalls_net.tntp → 边列表。返回 links=[(u,v,cap,t0),...] (节点1-indexed)"""
    links = []
    with open(path) as f:
        started = False
        for line in f:
            line = line.strip()
            if line.startswith("<END OF METADATA>"):
                started = True; continue
            if not started or not line or line.startswith("~"):
                continue
            parts = line.replace(";", "").split()
            if len(parts) < 6:
                continue
            u, v = int(parts[0]), int(parts[1])
            cap, t0 = float(parts[2]), float(parts[4])  # capacity, free_flow_time
            links.append((u, v, cap, t0))
    return links


def parse_trips(path):
    """解析 SiouxFalls_trips.tntp → OD需求字典 {(o,d): demand}"""
    od = {}
    with open(path) as f:
        origin = None
        for line in f:
            line = line.strip()
            if line.startswith("Origin"):
                origin = int(line.split()[1]); continue
            if origin is not None and ";" in line:
                for tok in line.split(";"):
                    tok = tok.strip()
                    if ":" in tok:
                        d, dem = tok.split(":")
                        d = int(d.strip()); dem = float(dem.strip())
                        if dem > 0:
                            od[(origin, d)] = dem
    return od


def build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0):
    """
    构建 Sioux Falls: 边索引化 + 各OD的k最短路径(按边索引)。
    k_paths: 每OD候选路径数; top_od: 只取需求最大的前N个OD对(None=全部,76*可能太多);
    返回: dict(edges,t0,cap, od_list) 其中 od_list=[{'o','d','demand','paths':[[e,..],..]},..]
    """
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    od = parse_trips(os.path.join(DATA_DIR, "SiouxFalls_trips.tntp"))

    # 边索引化
    edge_idx = {}
    t0_list = []; cap_list = []
    G = nx.DiGraph()
    for i, (u, v, cap, t0) in enumerate(links):
        edge_idx[(u, v)] = i
        t0_list.append(t0); cap_list.append(cap)
        G.add_edge(u, v, weight=t0, eidx=i)
    E = len(links)

    # 选OD: 按需求排序取前 top_od (减小规模)
    od_sorted = sorted(od.items(), key=lambda kv: -kv[1])
    if top_od is not None:
        od_sorted = od_sorted[:top_od]

    od_list = []
    for (o, d), dem in od_sorted:
        if o == d or not nx.has_path(G, o, d):
            continue
        # k最短简单路径(按free-flow time)
        paths_nodes = []
        try:
            gen = nx.shortest_simple_paths(G, o, d, weight="weight")
            for j, p in enumerate(gen):
                paths_nodes.append(p)
                if j + 1 >= k_paths:
                    break
        except nx.NetworkXNoPath:
            continue
        # 节点路径 → 边索引路径
        paths_edges = []
        for pn in paths_nodes:
            pe = [edge_idx[(pn[i], pn[i+1])] for i in range(len(pn)-1)]
            paths_edges.append(pe)
        od_list.append({'o': o, 'd': d, 'demand': dem*demand_scale, 'paths': paths_edges})

    return dict(edges=E, t0=np.array(t0_list), cap=np.array(cap_list), od_list=od_list)


if __name__ == "__main__":
    print("解析 Sioux Falls...")
    links = parse_net(os.path.join(DATA_DIR, "SiouxFalls_net.tntp"))
    od = parse_trips(os.path.join(DATA_DIR, "SiouxFalls_trips.tntp"))
    print(f"  边数={len(links)}, OD对数={len(od)}, 总需求={sum(od.values()):.0f}")
    print(f"  自由流时间范围=[{min(l[3] for l in links):.1f},{max(l[3] for l in links):.1f}]")
    print(f"  容量范围=[{min(l[2] for l in links):.0f},{max(l[2] for l in links):.0f}]")
    print()
    print("构建(k=3路径, 取需求最大前20个OD):")
    sf = build_sioux_falls(k_paths=3, top_od=20)
    print(f"  边数={sf['edges']}, OD数={len(sf['od_list'])}")
    print("  示例OD:")
    for od_item in sf['od_list'][:3]:
        print(f"    ({od_item['o']}→{od_item['d']}) demand={od_item['demand']:.0f}, "
              f"{len(od_item['paths'])}条路径, 路径长度={[len(p) for p in od_item['paths']]}")
