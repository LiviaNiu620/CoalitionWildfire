"""
通用网络混合原子-非原子均衡求解器。
支持: 任意网络(边成本函数) + 多 OD + 任意数量公司(各自 θ) + 非原子 HDV。

主体:
  - HDV: 非原子自私 (Wardrop), 每个 OD 一份需求
  - 公司 k: 原子, 质量 m_k, θ_k. 感知边际 = 路径成本 + 内部化项 - θ_k*对手项
    内部化项(路径 r) = Σ_{e∈r} f^k_e * c'_e(x_e)
    对手项(路径 r)   = Σ_{k'≠k} (1/m_{k'}) Σ_{e∈r} f^{k'}_e * c'_e(x_e)

网络用: edges 成本 c_e(x)=free_e + slope_e * x  (线性, 便于 c'_e=slope_e)
        paths: 每个 OD 的路径列表, 每条路径是边索引列表
"""
import numpy as np

EPS = 1e-12


class Network:
    def __init__(self, n_edges, free, slope, ods):
        """
        n_edges: 边数
        free[e], slope[e]: 边成本 c_e(x)=free[e]+slope[e]*x
        ods: list of dict, 每个 OD = {'demand_hdv':.., 'paths':[[e,..],..]}
             公司需求单独给
        """
        self.n_edges = n_edges
        self.free = np.array(free, float)
        self.slope = np.array(slope, float)
        self.ods = ods

    def edge_cost(self, x):
        return self.free + self.slope * x

    def edge_cprime(self, x):
        return self.slope  # 线性


def solve(net, company_demands, thetas, iters=15000, lam=0.01, seed=0, verbose=False):
    """
    company_demands: (n_od, K) 每个 OD 每家公司的需求量
    thetas: (K,) 各公司 θ
    返回 均衡流量与系统成本 J
    """
    n_od = len(net.ods)
    K = len(thetas)
    rng = np.random.default_rng(seed)

    # 决策: 每个 OD, HDV 路径流 + 每公司路径流
    # 初始化: 均分
    fH = []   # fH[o] = array over paths
    F = []    # F[o] = (K, n_paths_o)
    for o, od in enumerate(net.ods):
        npath = len(od['paths'])
        fH.append(np.full(npath, od['demand_hdv'] / npath))
        Fo = np.zeros((K, npath))
        for k in range(K):
            Fo[k] = company_demands[o][k] / npath
        F.append(Fo)

    def proj_simplex(v, total):
        if total <= EPS:
            return np.zeros_like(v)
        if len(v) == 1:
            return np.array([total])
        u = np.sort(v)[::-1]
        css = np.cumsum(u) - total
        idx = np.arange(1, len(u) + 1)
        cond = u - css / idx > 0
        r = np.nonzero(cond)[0][-1]
        theta_ = css[r] / (r + 1.0)
        return np.maximum(v - theta_, 0.0)

    def compute_edge_flows():
        x = np.zeros(net.n_edges)          # 总流量
        xk = np.zeros((K, net.n_edges))    # 各公司流量
        for o, od in enumerate(net.ods):
            for p, path in enumerate(od['paths']):
                for e in path:
                    x[e] += fH[o][p] + F[o][:, p].sum()
                    xk[:, e] += F[o][:, p]
        return x, xk

    for it in range(iters):
        x, xk = compute_edge_flows()
        ce = net.edge_cost(x)
        cp = net.edge_cprime(x)

        for o, od in enumerate(net.ods):
            npath = len(od['paths'])
            # 路径成本
            pc = np.array([sum(ce[e] for e in od['paths'][p]) for p in range(npath)])
            # HDV 自私
            demH = od['demand_hdv']
            fH[o] = proj_simplex(fH[o] - lam * pc * demH, demH)
            # 各公司
            for k in range(K):
                mk = company_demands[o][k]
                if mk < EPS:
                    continue
                # 内部化项: Σ_{e∈r} f^k_e * c'_e
                intern = np.array([sum(xk[k, e] * cp[e] for e in od['paths'][p]) for p in range(npath)])
                # 对手项
                opp = np.zeros(npath)
                for kp in range(K):
                    if kp == k:
                        continue
                    mkp = company_demands[o][kp]
                    if mkp < EPS:
                        continue
                    opp += np.array([sum(xk[kp, e] * cp[e] for e in od['paths'][p]) for p in range(npath)]) / mkp
                grad = pc + intern - thetas[k] * opp
                F[o][k] = proj_simplex(F[o][k] - lam * grad * mk, mk)

    # 系统成本
    x, xk = compute_edge_flows()
    ce = net.edge_cost(x)
    J = float(np.dot(x, ce))
    return dict(x=x, J=J, fH=fH, F=F)


# ============ 网络定义 ============
def make_braess_extended():
    """扩展 Braess: s->t, 多条平行 + 捷径。
    边: 0=sv(拥堵),1=wt(拥堵),2=sw(常数1),3=vt(常数1),4=vw(捷径0)
    路径: P1=[0,3], P2=[2,1], P3=[0,4,1]"""
    free  = [0, 0, 1, 1, 0]
    slope = [1, 1, 0, 0, 0]
    ods = [{'demand_hdv': None, 'paths': [[0,3],[2,1],[0,4,1]]}]
    return Network(5, free, slope, ods)


def make_two_route(a=0.6):
    """两路: 拥堵路 c=x, 备选 c=a. 边0=拥堵(slope1),边1=备选(free a)"""
    free=[0, a]; slope=[1, 0]
    ods=[{'demand_hdv':None,'paths':[[0],[1]]}]
    return Network(2, free, slope, ods)


def make_multi_od_grid():
    """简单 2-OD 网络共享瓶颈。
    OD1: s1->t 经 共享拥堵边 e0 或 绕行 e1(常数)
    OD2: s2->t 经 共享拥堵边 e0 或 绕行 e2(常数)
    e0 被两个 OD 共享(耦合)。"""
    free=[0, 1.0, 1.0]; slope=[1,0,0]
    ods=[
        {'demand_hdv':None,'paths':[[0],[1]]},  # OD1: 拥堵 or 绕行1
        {'demand_hdv':None,'paths':[[0],[2]]},  # OD2: 拥堵 or 绕行2
    ]
    return Network(3, free, slope, ods)


if __name__ == "__main__":
    print("通用求解器自检: 扩展Braess, α=0.5, 两公司对称θ")
    net = make_braess_extended()
    alpha=0.5
    net.ods[0]['demand_hdv'] = 1-alpha
    for th in [-1,0,1]:
        cd = [[alpha/2, alpha/2]]
        eq = solve(net, cd, np.array([th,th]))
        print(f"  θ={th:+d}: J={eq['J']:.4f}")
