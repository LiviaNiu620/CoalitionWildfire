"""
两子群体模型 (正确微观基础) 通用求解器 —— Wardrop 一致版本.

单瓶颈 canonical: 拥堵链路 c1(x1)=t0+t1*x1, 备选链路 c2=a (constant, ample capacity).
  - 运营商 k, 质量 m_k, Σ m_k = α N. 协调子队容量 y_k = β(m_k) m_k, 内部化"自身协调流量".
  - 未协调子队 (1-β_k)m_k 与 HDV (1-α)N 合并为自私体 W, 走 Wardrop.

均衡两个 regime (第四轮审稿修正: 角点方向)
-------------------------------------------------
自私体是非原子的, 面对 c1(x1) vs a:
  (S) 自私体内点(两条链路都用) <=> c1(x1)=a <=> x1=x^UE. 此时 μ=0, 全部 g_k=0.
      自洽条件: 自私体质量足够填满瓶颈到 UE 水平, 即 W >= x^UE. 结构完全无效, ρ=0.
  (C) 自私体角点. 关键: 角点方向是"自私体全部在**瓶颈**上"(便宜的那条),
      不是全部在备选上. 若 x1 < x^UE 则 c1(x1) < a, 便宜链路上放零流量不可能是
      Wardrop 均衡. 自洽条件 W < x^UE. 此时
          x1 = W + Σ_k min(μ, y_k),   μ = x^UE - x1.

      旧版 (technote v1, Lemma 1) 写的是 x1 = Σ_k g_k, 相容条件 "iff c1(x1) <= a",
      方向反了: 那个 regime 除退化点外是空集 (要求同时 x1>=x^UE 与 x1<=x^UE, 而
      x1=x^UE 处 μ=0 => g_k=0 => x1=0). 因此自私体质量 W 必须进入闭式.

闭式 (regime C)
-------------------------------------------------
    x1 = (H + n_int · x^UE) / (1 + n_int),
    H  = (1-α)N + Σ_{k∉D} m_k + Σ_{k∈D} (1-β(m_k)) m_k
       = N - Σ_{k∈D} β(m_k) m_k          (等价化简, 由 Σ m_k = αN)
其中 D = {k : y_k > μ} 为内点集, n_int = |D|.
注意角点运营商对 x1 的贡献是 (1-β_k)m_k + y_k = m_k, 与 β_k 无关.

分类不动点 (Prop 1 的存在唯一性)
-------------------------------------------------
    G(μ) := μ + W + Σ_k min(μ, y_k) - x^UE
G 连续、分段线性、斜率 >= 1 故严格递增; G(0)=W-x^UE, G(x^UE)=W+Σ min(x^UE,y_k) >= 0.
所以 W < x^UE 时在 (0, x^UE] 上有唯一根, 用二分求解 (非 damped iteration).
W >= x^UE 时根落在 μ<0, 边界解即 regime (S).

适定性前提
-------------------------------------------------
恢复率 ρ = (J^UE - J^eq)/(J^UE - J^SO) 的标准解释要求 **N >= x^UE**, 否则 UE 下备选
链路根本不被使用 (全部需求都在瓶颈上且成本仍低于 a), J^UE / x^UE 失去意义.
本求解器对 N < x^UE 报错而非静默返回.
"""
import numpy as np

__all__ = ["beta_exp", "solve_two_subgroup", "critical_scale", "shutdown_K",
           "alpha_activation", "solve_legacy_note_v1"]


def beta_exp(m, lam):
    """代表性协调度 β(m) = 1 - exp(-λ m): 满足 (B1) β(0)=0, (B2) β'>0, (B3) β''<=0."""
    return 1 - np.exp(-lam * np.asarray(m, dtype=float))


def _cost_J(x1, N, t0, t1, a):
    """总系统成本 J(x1) = x1 c1(x1) + (N - x1) a."""
    return x1 * (t0 + t1 * x1) + (N - x1) * a


def solve_two_subgroup(masses, lam=None, alpha=1.0, N=None, t0=0.0, t1=1.0, a=1.0,
                       beta_func=None, check_demand=True):
    """求解两子群体混合 Cournot-Nash / Wardrop 均衡.

    参数
    ----
    masses      : 各运营商车队质量 m_k. 应满足 Σ m_k = α N (不满足时按 alpha, N 给定值
                  处理 HDV 质量, 并在返回值里给出 mass_residual 供调用方检查).
    lam         : β(m)=1-exp(-λm) 的 λ (当 beta_func 为 None 时使用).
    alpha       : AV 渗透率. AV 质量 = α N, HDV 质量 = (1-α) N.
    N           : 总需求. 必须显式给出, 且需 N >= x^UE (见模块 docstring).
    t0, t1, a   : 瓶颈 c1(x)=t0+t1 x, 备选常数成本 a. x^UE=(a-t0)/t1, x^SO=x^UE/2.
    beta_func   : 自定义协调度函数 (向量化, 接受 ndarray). 给出时忽略 lam.
    check_demand: 是否强制 N >= x^UE.

    返回 dict
    ----
    x1, rho, regime('S'/'C'), mu, n_int, D(内点下标), H, W, B_cap, y, g,
    betas, xUE, xSO, mass_residual
    """
    masses = np.asarray([m for m in np.atleast_1d(masses) if m > 1e-14], dtype=float)
    if beta_func is None:
        if lam is None:
            raise ValueError("需给出 lam 或 beta_func")
        beta_func = lambda m: beta_exp(m, lam)
    if N is None:
        raise ValueError("N 必须显式给出 (旧版默认 N=1.0 会落在 N<x^UE 的失效区)")

    xUE = (a - t0) / t1
    xSO = xUE / 2.0                       # argmin J: t0 + 2 t1 x = a  =>  x = x^UE/2
    if check_demand and N < xUE - 1e-12:
        raise ValueError(
            f"N={N} < x^UE={xUE}: UE 下备选链路未被使用, ρ=(J^UE-J)/(J^UE-J^SO) 无标准解释. "
            "请取 N >= x^UE, 或显式传 check_demand=False 并自行解释结果.")

    betas = beta_func(masses)
    y = betas * masses                                     # 协调子队容量
    W = (1 - alpha) * N + float(np.sum(masses - y))        # 自私体: HDV + 全部未协调子队

    if W >= xUE - 1e-15:
        # regime (S): 自私体内点, Wardrop 钉 c1=a, 结构完全无效
        mu, x1, regime = 0.0, xUE, 'S'
        g = np.zeros_like(y)
        D = np.array([], dtype=int)
    else:
        # regime (C): 自私体全在瓶颈; 对 μ 单调二分 G(μ)=μ+W+Σmin(μ,y_k)-x^UE
        G = lambda m_: m_ + W + float(np.sum(np.minimum(m_, y))) - xUE
        lo, hi = 0.0, xUE
        assert G(lo) < 0.0 <= G(hi) + 1e-15
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if G(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        mu = 0.5 * (lo + hi)
        g = np.minimum(mu, y)
        x1 = W + float(np.sum(g))
        regime = 'C'
        D = np.flatnonzero(y > mu + 1e-12)                  # 内点集

    n_int = int(D.size)
    B_cap = float(np.sum(y[y <= mu + 1e-12])) if y.size else 0.0
    # H 的两种等价写法, 互为交叉校验
    corner = np.setdiff1d(np.arange(masses.size), D)
    H_3term = (1 - alpha) * N + float(np.sum(masses[corner])) \
        + float(np.sum(masses[D] - y[D]))
    H_short = N - float(np.sum(y[D]))
    resid = abs(N - (1 - alpha) * N - float(np.sum(masses)))   # Σm_k 是否 = αN
    if resid < 1e-9:
        assert abs(H_3term - H_short) < 1e-9, (H_3term, H_short)
    H = H_3term

    J_UE, J_SO, J_eq = _cost_J(xUE, N, t0, t1, a), _cost_J(xSO, N, t0, t1, a), _cost_J(x1, N, t0, t1, a)
    rho = (J_UE - J_eq) / (J_UE - J_SO) if J_UE - J_SO > 1e-15 else 0.0

    return dict(x1=x1, rho=rho, regime=regime, mu=mu, n_int=n_int, D=D,
                H=H, W=W, B_cap=B_cap, y=y, g=g, betas=betas,
                xUE=xUE, xSO=xSO, mass_residual=resid)


# --------------------------------------------------------------------------
# 修正模型下的三个阈值 (替代 v1 的 m^\star / Θ(1/K))
# --------------------------------------------------------------------------
def _bisect(f, lo, hi, iters=200):
    assert f(lo) <= 0 <= f(hi), (f(lo), f(hi))
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def critical_scale(lam, N, t0=0.0, t1=1.0, a=1.0, beta_func=None, hi=1e3):
    """单运营商真 kink 的临界规模 m†: β(m†) m† = N - x^UE.

    m < m† 时 W >= x^UE, 结构完全失效, ρ ≡ 0 (不是"打折");
    m > m† 时 ρ 从 0 光滑上升. 与 v1 的 β(m*)m* = x^UE/2 是不同的方程 —— 二者仅在
    N - x^UE = x^UE/2 (即 N = 1.5 x^UE) 时数值巧合相同.
    """
    if beta_func is None:
        beta_func = lambda m: beta_exp(m, lam)
    xUE = (a - t0) / t1
    target = N - xUE
    if target <= 0:
        return 0.0
    return _bisect(lambda m: beta_func(np.array([m]))[0] * m - target, 1e-12, hi)


def shutdown_K(lam, alpha, N, t0=0.0, t1=1.0, a=1.0, beta_func=None, Kmax=100000):
    """碎片化有限关停阈值 K̄: 最小的 K 使等分 K 家时 W(K) >= x^UE, 此后 ρ = 0 精确.

    W(K) = (1-α)N + (1 - β(αN/K)) αN, 随 K 递增并趋于 N. 返回 None 表示 K<=Kmax 内不关停.
    """
    if beta_func is None:
        beta_func = lambda m: beta_exp(m, lam)
    xUE = (a - t0) / t1
    for K in range(1, Kmax + 1):
        m = alpha * N / K
        W = (1 - alpha) * N + (1 - beta_func(np.array([m]))[0]) * alpha * N
        if W >= xUE:
            return K
    return None


def alpha_activation(masses_profile, lam, N, t0=0.0, t1=1.0, a=1.0, beta_func=None):
    """给定车队规模剖面, 结构起作用(H < x^UE)所需的激活渗透率.

    注意: 条件 (1-α)N + Σ(1-β(m_k))m_k < x^UE 中 Σ m_k = αN 本身依赖 α, 所以这是
    **隐式**阈值. 本函数按"以规模剖面为自变量"读法返回
        α_c = 1 - [x^UE - Σ(1-β(m_k))m_k] / N,
    并同时返回 Σ m_k / N 供调用方检查是否自洽 (α_c 应与之相符才是一个真解).
    """
    if beta_func is None:
        beta_func = lambda m: beta_exp(m, lam)
    ms = np.asarray(masses_profile, dtype=float)
    xUE = (a - t0) / t1
    uncoord = float(np.sum((1 - beta_func(ms)) * ms))
    return dict(alpha_c=1.0 - (xUE - uncoord) / N,
                alpha_implied=float(np.sum(ms)) / N,
                uncoordinated_mass=uncoord)


# --------------------------------------------------------------------------
# v1 (technote 第一版) 的错误规则 —— 仅供回归测试记录"改了什么", 不要用于任何新结果
# --------------------------------------------------------------------------
def solve_legacy_note_v1(masses, lam, t0=0.0, t1=1.0, a=1.0, beta_func=None):
    """复现 technote v1 Lemma 1 / Prop 1: x1 = Σ_k min(μ, y_k), 自私体质量 W 不进入.

    保留仅为文档化第四轮审稿修正前后的差异 (例如 A=(0.9,0.4,0.2) 给 ρ=0.918).
    该规则的 regime 是空集, 结论不可用.
    """
    masses = np.asarray([m for m in np.atleast_1d(masses) if m > 1e-14], dtype=float)
    if beta_func is None:
        beta_func = lambda m: beta_exp(m, lam)
    xUE = (a - t0) / t1
    y = beta_func(masses) * masses
    G = lambda m_: m_ + float(np.sum(np.minimum(m_, y))) - xUE
    lo, hi = 0.0, xUE
    if G(lo) >= 0:
        mu = 0.0
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if G(mid) < 0:
                lo = mid
            else:
                hi = mid
        mu = 0.5 * (lo + hi)
    x1 = min(xUE, float(np.sum(np.minimum(mu, y))))
    r = x1 / a
    return dict(x1=x1, rho=4 * r * (1 - r), mu=mu,
                n_int=int(np.sum(y > mu + 1e-12)),
                B_cap=float(np.sum(y[y <= mu + 1e-12])))


if __name__ == "__main__":
    N, alpha, lam = 1.5, 1.0, 2.0
    print(f"canonical: a=t1=1, t0=0 => x^UE=1, x^SO=0.5;  N={N}, alpha={alpha}, lam={lam}")
    print(f"{'config':<20}{'regime':<8}{'W':>9}{'H':>9}{'n_int':>7}{'x1':>10}{'rho':>10}")
    for nm, ms in [("monopoly 1.5", [1.5]), ("2 x 0.75", [0.75] * 2),
                   ("3 x 0.5", [0.5] * 3), ("(0.9,0.4,0.2)", [0.9, 0.4, 0.2]),
                   ("(1.0,0.5)", [1.0, 0.5]), ("8 x 0.1875", [1.5 / 8] * 8)]:
        r = solve_two_subgroup(ms, lam=lam, alpha=alpha, N=N)
        print(f"{nm:<20}{r['regime']:<8}{r['W']:>9.4f}{r['H']:>9.4f}"
              f"{r['n_int']:>7}{r['x1']:>10.6f}{r['rho']:>10.6f}")
    print(f"\nm_dagger (N={N}) = {critical_scale(lam, N):.6f}"
          f"   [v1 的 m* = 0.675 只是 N=1.5 下的数值巧合]")
    print(f"K_bar (alpha={alpha}, N={N}) = {shutdown_K(lam, alpha, N)}  => K >= K_bar 时 rho = 0 精确")
