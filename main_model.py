"""
完备版主模型（§9-10）验证：统一 θ 轴 + AV 容量效应 + 三类主体混合均衡。

主体：
  - HDV: 非原子自私，Wardrop（感知成本 = 实际路径成本）
  - 公司 k（原子，质量 m_k）: 目标 = 自己车队成本 - θ_k·(对手车队平均成本)
      θ_k < 0 利他 / θ_k = 0 规模内部化 / θ_k > 0 竞争抢占
    感知边际(路径 r) = C_r + Σ_{e∈r} f^k_e·(∂c_e/∂x^AV) - θ_k·Σ_{k'≠k}(1/m_{k'})Σ_{e∈r} f^{k'}_e·(∂c_e/∂x^AV)

链路成本（AV 容量效应, §10）:
    c_e = c0_e * ( (x^HDV_e + κ(ρ_e)·x^AV_e) / Cap_e ) + free_e   [线性 latency 便于解析]
    ρ_e = x^AV_e / x_e ,  κ(ρ) = 1 - b·ρ  (b∈[0,1), AV 越密越省空间)
    b=0 时退化为无容量效应（AV 与 HDV 同等占空间）。

网络：Braess（含捷径），单 OD，需求=1。
路径: P1=[sv,vt], P2=[sw,wt], P3(捷径)=[sv,vw,wt]
拥堵边 sv,wt: c0=1, Cap=1, free=0 ; 常数边 sw,vt: 恒=1 ; 捷径 vw: 恒=0
"""

import numpy as np

DEMAND = 1.0
EPS = 1e-12

# ---- 边定义：只有 sv, wt 是拥堵边（有容量效应），sw,vt 常数, vw=0 ----
# 路径-边关联
# 边索引: 0=sv, 1=wt, 2=sw, 3=vt, 4=vw
# P1 = sv + vt ; P2 = sw + wt ; P3 = sv + vw + wt
PATH_EDGES = {
    0: [0, 3],      # P1: sv, vt
    1: [2, 1],      # P2: sw, wt
    2: [0, 4, 1],   # P3: sv, vw, wt
}
CONGESTED = {0, 1}  # sv, wt 有容量效应
CONST_EDGE = {2: 1.0, 3: 1.0, 4: 0.0}  # sw=1, vt=1, vw=0


def edge_flows(fH, F_av):
    """给定 HDV 路径流 fH(3,) 和 AV 各公司路径流 F_av(K,3)，返回每条边的 HDV 流、AV 流。
       返回 xH[5], xA[5]."""
    xH = np.zeros(5)
    xA = np.zeros(5)
    for p in range(3):
        for e in PATH_EDGES[p]:
            xH[e] += fH[p]
            xA[e] += F_av[:, p].sum()
    return xH, xA


def edge_cost_and_grad(xH_e, xA_e, b):
    """拥堵边成本 c_e 及 ∂c_e/∂x^AV。线性: c = (xH + κ(ρ)·xA)/Cap, Cap=1, c0=1.
       κ(ρ)=1-b·ρ, ρ=xA/(xH+xA)."""
    x = xH_e + xA_e
    if x < EPS:
        return 0.0, 1.0  # 空边，边际=κ(0)=1
    rho = xA_e / x
    kappa = 1.0 - b * rho
    c = xH_e + kappa * xA_e            # = effective load (Cap=1,c0=1)
    # ∂c/∂xA: 对 xA 求导，注意 κ 依赖 ρ 依赖 xA
    # c = xH + (1 - b*xA/x)*xA = xH + xA - b*xA^2/x
    # dc/dxA = 1 - b*(2*xA*x - xA^2)/x^2 = 1 - b*(2ρ - ρ^2) = 1 - b*ρ*(2-ρ)
    dc_dxA = 1.0 - b * rho * (2.0 - rho)
    return c, dc_dxA


def path_actual_costs(xH, xA, b):
    """每条路径的实际成本（用于 HDV Wardrop、系统成本）。"""
    ce = np.zeros(5)
    for e in range(5):
        if e in CONGESTED:
            ce[e], _ = edge_cost_and_grad(xH[e], xA[e], b)
        else:
            ce[e] = CONST_EDGE[e]
    pc = np.array([sum(ce[e] for e in PATH_EDGES[p]) for p in range(3)])
    return pc, ce


def path_marginal_AV(xH, xA, b):
    """每条路径上 Σ_{e∈r,congested} ∂c_e/∂x^AV （用于原子公司的内部化项）。
       常数边 ∂/∂xA=0。返回 g[3]."""
    dce = np.zeros(5)
    for e in CONGESTED:
        _, dce[e] = edge_cost_and_grad(xH[e], xA[e], b)
    g = np.array([sum(dce[e] for e in PATH_EDGES[p]) for p in range(3)])
    return g


def solve_equilibrium(m, theta, b, iters=40000, lam=0.005, seed=0, f_init=None):
    """
    m: (K,) 各公司质量; theta: (K,) 各公司 θ; b: 容量斜率.
    返回均衡与系统成本。投影梯度（到单纯形）迭代。
    """
    K = len(m)
    alpha = m.sum()
    hdv = DEMAND - alpha
    rng = np.random.default_rng(seed)

    # 初始化
    if f_init is None:
        fH = np.full(3, hdv / 3)
        F = np.array([np.full(3, mk / 3) for mk in m])  # (K,3)
    else:
        fH, F = f_init

    def proj_simplex(v, total):
        if total <= EPS:
            return np.zeros_like(v)
        u = np.sort(v)[::-1]
        css = np.cumsum(u) - total
        idx = np.arange(1, len(u) + 1)
        cond = u - css / idx > 0
        rho_k = np.nonzero(cond)[0][-1]
        theta_ = css[rho_k] / (rho_k + 1.0)
        return np.maximum(v - theta_, 0.0)

    for it in range(iters):
        xH, xA = edge_flows(fH, F)
        pc, _ = path_actual_costs(xH, xA, b)          # 路径实际成本
        gAV = path_marginal_AV(xH, xA, b)             # 路径 AV 边际拥堵

        # 每个公司在各边上自己的 AV 流量（用于内部化项 Σ f^k_e · ∂c/∂xA）
        # 近似：用路径级 gAV × 公司自己路径流的边贡献。精确做边级：
        # 公司 k 路径 r 的内部化 = Σ_{e∈r} f^k_e · dc_e/dxA
        # 先算各公司边级 AV 流
        # 边级：fk_e = Σ_{r∋e} F[k,r]
        dce = np.zeros(5)
        for e in CONGESTED:
            dce[e] = edge_cost_and_grad(xH[e], xA[e], b)[1]

        # HDV: 自私 Wardrop 梯度 = pc
        fH = proj_simplex(fH - lam * pc * hdv * 3, hdv)  # 缩放步长

        # 各公司
        Fnew = F.copy()
        for k in range(K):
            # 公司 k 边级自身 AV 流
            fk_e = np.zeros(5)
            for p in range(3):
                for e in PATH_EDGES[p]:
                    fk_e[e] += F[k, p]
            # 内部化项 per path
            intern = np.array([sum(fk_e[e] * dce[e] for e in PATH_EDGES[p] if e in CONGESTED) for p in range(3)])
            # 对手项 per path: θ_k · Σ_{k'≠k}(1/m_{k'}) Σ_{e∈r} f^{k'}_e·dc_e
            opp = np.zeros(3)
            for kp in range(K):
                if kp == k or m[kp] < EPS:
                    continue
                fkp_e = np.zeros(5)
                for p in range(3):
                    for e in PATH_EDGES[p]:
                        fkp_e[e] += F[kp, p]
                opp_kp = np.array([sum(fkp_e[e] * dce[e] for e in PATH_EDGES[p] if e in CONGESTED) for p in range(3)])
                opp += opp_kp / m[kp]
            grad_k = pc + intern - theta[k] * opp
            Fnew[k] = proj_simplex(F[k] - lam * grad_k * m[k] * 3, m[k])
        F = Fnew

    xH, xA = edge_flows(fH, F)
    pc, ce = path_actual_costs(xH, xA, b)
    # 系统总成本 = Σ_edge x_e·c_e  (实际行驶成本)
    xtot = xH + xA
    J = 0.0
    for e in range(5):
        if e in CONGESTED:
            J += xtot[e] * ce[e]
        else:
            J += xtot[e] * CONST_EDGE[e]
    return dict(fH=fH, F=F, xH=xH, xA=xA, pc=pc, J=J)


def system_optimum(alpha, b):
    """网格搜索社会最优 J（所有流量由规划者分配，含容量效应）。"""
    best = None
    grid = np.linspace(0, 1, 121)
    # 简化：假设 AV 全部集中最优分配，HDV 也被规划。用总流量在 P1/P2/P3 分配 + AV占比。
    # 这里只做粗略下界参考：全部车按最优路径分配（避免捷径），AV 占比 alpha。
    for f3 in np.linspace(0, 1, 61):
        rem = 1 - f3
        for f1 in np.linspace(0, rem, 61):
            f2 = rem - f1
            fH = np.array([f1, f2, f3]) * (1 - alpha)
            # AV 也按同比例（规划者对称分配）
            FA = np.array([[f1, f2, f3]]) * alpha
            xH, xA = edge_flows(fH, FA)
            _, ce = path_actual_costs(xH, xA, b)
            xtot = xH + xA
            J = sum(xtot[e] * (ce[e] if e in CONGESTED else CONST_EDGE[e]) for e in range(5))
            if best is None or J < best:
                best = J
    return best


if __name__ == "__main__":
    alpha = 0.5
    b = 0.0  # 先无容量效应，隔离 θ 轴效果
    print(f"=== Braess, alpha={alpha}, b(容量斜率)={b} ===")
    print("Wardrop(全自私) 参考: 全体 θ=0 且 m→0")
    eqW = solve_equilibrium(np.array([1e-6]*1), np.array([0.0]), b)  # 近似
    print(f"  (退化) J≈{eqW['J']:.4f}")
    print()
    # 两家公司，各 alpha/2，扫描 θ 从 全规模型(0,0) 到 全竞争型
    print("两家公司各 m=alpha/2, 扫描竞争强度 θ (对称)")
    print("  θ      xA_shortcut   J")
    for th in [-1.0, -0.5, 0.0, 0.5, 1.0, 2.0]:
        m = np.array([alpha/2, alpha/2])
        theta = np.array([th, th])
        eq = solve_equilibrium(m, theta, b)
        av_short = eq['F'][:, 2].sum()
        print(f"  {th:+.2f}   {av_short:.4f}       {eq['J']:.5f}")
