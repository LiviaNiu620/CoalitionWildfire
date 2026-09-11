"""
向量化 BPR 求解器 (numpy矩阵化, 支持真实Sioux Falls规模)。
用 路径-边关联矩阵 A (边×路径) 批量运算, 避免Python循环。
"""
import numpy as np
from sioux_falls_loader import build_sioux_falls


def solve_bpr_vec(E, t0, cap, od_list, alpha, s, K_I=1, cap_b=0.0,
                  bpr_a=0.15, bpr_b=4.0, iters=4000, lam=1e-6):
    """
    向量化求解混合原子-非原子均衡。
    每个OD: HDV(非原子Wardrop) + K_I家规模型(β=1) + 1家竞争型(β=0)。
    返回 系统总成本 J。
    """
    # 构建全局 路径列表 与 关联矩阵
    # 每条路径归属一个OD和一个"玩家类型群": 我们对每个(OD)分别持有 HDV流 + 各公司流。
    # 为向量化: 把所有OD所有路径拼成一个大路径集, 记录每条路径的边集(关联矩阵A: E×P),
    # 以及每条路径属于哪个OD、哪个玩家。
    paths = []       # 每条路径的边索引列表
    path_od = []     # 路径→OD编号
    for oi, od in enumerate(od_list):
        for pe in od['paths']:
            paths.append(pe); path_od.append(oi)
    P = len(paths)
    path_od = np.array(path_od)
    # 关联矩阵 A: E×P, A[e,p]=1 若边e在路径p
    A = np.zeros((E, P))
    for p, pe in enumerate(paths):
        for e in pe:
            A[e, p] = 1.0

    n_od = len(od_list)
    # 玩家: HDV + K_I规模型 + 1竞争型, 共 K=K_I+1 家公司
    K = K_I + 1
    betas = np.array([1.0]*K_I + [0.0])
    # 各OD各玩家的需求
    dem_H = np.array([(1-alpha)*od['demand'] for od in od_list])           # (n_od,)
    dem_firm = np.zeros((n_od, K))
    for oi, od in enumerate(od_list):
        d = od['demand']
        for k in range(K_I): dem_firm[oi,k] = alpha*s*d/K_I
        dem_firm[oi,K-1] = max(alpha*(1-s)*d, 1e-9)

    # 路径流: HDV路径流 fH(P,), 各公司路径流 fF(K,P)
    # 每OD内均分初始
    fH = np.zeros(P); fF = np.zeros((K,P))
    for oi in range(n_od):
        mask = path_od==oi; npth=mask.sum()
        fH[mask] = dem_H[oi]/npth
        for k in range(K): fF[k,mask] = dem_firm[oi,k]/npth

    # OD→路径 的分组索引
    od_masks = [ (path_od==oi) for oi in range(n_od) ]

    def bpr(x, xav):
        c = cap.copy()
        if cap_b>0:
            rho = xav/np.maximum(x,1e-12)
            c = cap*(1+cap_b*rho)
        return t0*(1+bpr_a*(x/c)**bpr_b)
    def bpr_prime(x, xav):
        c = cap.copy()
        if cap_b>0:
            rho = xav/np.maximum(x,1e-12)
            c = cap*(1+cap_b*rho)
        return t0*bpr_a*bpr_b*(x**(bpr_b-1))/(c**bpr_b)

    def proj_simplex_group(v, tot):
        if tot<=1e-12 or len(v)==0: return np.zeros_like(v)
        if len(v)==1: return np.array([tot])
        u=np.sort(v)[::-1]; css=np.cumsum(u)-tot; idx=np.arange(1,len(u)+1)
        pos=np.nonzero(u-css/idx>0)[0]
        if len(pos)==0: return np.maximum(v,0)
        r=pos[-1]; th=css[r]/(r+1.0)
        return np.maximum(v-th,0.0)

    for it in range(iters):
        # 边流量: x = A @ (fH + Σ fF)
        f_all = fH + fF.sum(axis=0)
        x = A @ f_all
        xav = A @ fF.sum(axis=0)
        ce = bpr(x, xav); cp = bpr_prime(x, xav)
        # 路径成本 pc = A^T @ ce  (P,)
        pc = A.T @ ce
        # 各公司边级流量 xk = A @ fF[k]
        # 内部化项(路径) = A^T @ (cp * (A@fF[k]))
        # HDV 更新
        gradH = pc
        for oi in range(n_od):
            m=od_masks[oi]
            fH[m] = proj_simplex_group(fH[m]-lam*gradH[m]*dem_H[oi], dem_H[oi])
        # 公司更新
        for k in range(K):
            xk = A @ fF[k]                       # (E,)
            intern = A.T @ (cp * xk)             # (P,) 内部化项
            grad = pc + betas[k]*intern
            for oi in range(n_od):
                m=od_masks[oi]
                mk = dem_firm[oi,k]
                if mk<1e-12: continue
                fF[k,m] = proj_simplex_group(fF[k,m]-lam*grad[m]*mk, mk)

    f_all = fH + fF.sum(axis=0)
    x = A @ f_all; xav = A @ fF.sum(axis=0)
    J = float(np.dot(x, bpr(x, xav)))
    return J


if __name__ == "__main__":
    import time
    print("加载 Sioux Falls (k=3, top12 OD, scale=0.05)...")
    sf = build_sioux_falls(k_paths=3, top_od=12, demand_scale=0.05)
    print(f"  边={sf['edges']}, OD={len(sf['od_list'])}")
    t=time.time()
    J = solve_bpr_vec(sf['edges'], sf['t0'], sf['cap'], sf['od_list'],
                      alpha=0.6, s=0.5, iters=3000, lam=1e-6)
    print(f"  单次求解 J={J:.1f}, 用时={time.time()-t:.1f}s")
