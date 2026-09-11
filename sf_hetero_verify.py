"""
一般网络(Sioux Falls)确认: 异质β(m)框架的 HHI-ρ 单调性 + 三regime.
基于bpr_vectorized核心, 改为: 任意规模份额{share_k} → 内生β_k=1-e^{-λ·m_k}.
"""
import numpy as np
from sioux_falls_loader import build_sioux_falls

def solve_hetero_bpr(sf, alpha, shares, lam_beta, bpr_a=0.15, bpr_b=4.0,
                     iters=5000, lam=1e-5):
    """
    shares: 各coalition占AV的份额(和=1). β_k=1-e^{-λ_β·m_k}, m_k=share_k·α·(总需求).
    返回 J.
    """
    E=sf['edges']; t0=sf['t0']; cap=sf['cap']; od_list=sf['od_list']
    paths=[]; path_od=[]
    for oi,od in enumerate(od_list):
        for pe in od['paths']:
            paths.append(pe); path_od.append(oi)
    P=len(paths); path_od=np.array(path_od)
    A=np.zeros((E,P))
    for p,pe in enumerate(paths):
        for e in pe: A[e,p]=1.0
    n_od=len(od_list)
    Kf=len(shares)
    total_dem=sum(od['demand'] for od in od_list)
    masses=np.array(shares)*alpha*total_dem   # 各coalition总质量(跨OD)
    betas=1-np.exp(-lam_beta*masses)
    # 各OD各玩家需求
    dem_H=np.array([(1-alpha)*od['demand'] for od in od_list])
    dem_firm=np.zeros((n_od,Kf))
    for oi,od in enumerate(od_list):
        for k in range(Kf):
            dem_firm[oi,k]=max(alpha*shares[k]*od['demand'],1e-12)
    fH=np.zeros(P); fF=np.zeros((Kf,P))
    for oi in range(n_od):
        m=path_od==oi; npth=m.sum()
        fH[m]=dem_H[oi]/npth
        for k in range(Kf): fF[k,m]=dem_firm[oi,k]/npth
    od_masks=[(path_od==oi) for oi in range(n_od)]
    def bpr(x): return t0*(1+bpr_a*(x/cap)**bpr_b)
    def bpr_p(x): return t0*bpr_a*bpr_b*(x**(bpr_b-1))/(cap**bpr_b)
    def proj(v,tot):
        if tot<=1e-12 or len(v)==0: return np.zeros_like(v)
        if len(v)==1: return np.array([tot])
        u=np.sort(v)[::-1]; css=np.cumsum(u)-tot; idx=np.arange(1,len(u)+1)
        pos=np.nonzero(u-css/idx>0)[0]
        if len(pos)==0: return np.maximum(v,0)
        r=pos[-1]; return np.maximum(v-css[r]/(r+1.0),0.0)
    for it in range(iters):
        f_all=fH+fF.sum(axis=0); x=A@f_all
        ce=bpr(x); cp=bpr_p(x); pc=A.T@ce
        for oi in range(n_od):
            m=od_masks[oi]; fH[m]=proj(fH[m]-lam*pc[m]*dem_H[oi],dem_H[oi])
        for k in range(Kf):
            xk=A@fF[k]; intern=A.T@(cp*xk); grad=pc+betas[k]*intern
            for oi in range(n_od):
                m=od_masks[oi]; mk=dem_firm[oi,k]
                if mk<1e-12: continue
                fF[k,m]=proj(fF[k,m]-lam*grad[m]*mk,mk)
    f_all=fH+fF.sum(axis=0); x=A@f_all
    return float(np.dot(x,bpr(x))), betas

def HHI(shares): s=np.array(shares); return (s**2).sum()
def dist(K,conc):
    if K==1: return [1.0]
    big=1.0/K+conc*(1-1.0/K); rest=(1-big)/(K-1); return [big]+[rest]*(K-1)

if __name__=="__main__":
    import time
    print("加载 Sioux Falls (k=2, top6 OD, scale=1.0 拥堵)...")
    sf=build_sioux_falls(k_paths=2, top_od=6, demand_scale=1.0)
    print(f"  边={sf['edges']}, OD={len(sf['od_list'])}")
    lam_beta=3e-4; alpha=0.9
    t=time.time()
    print(f"\nSioux Falls: HHI-ρ 单调性验证 (α={alpha}, λ_β={lam_beta})")
    print(f"{'配置':>18} {'HHI':>6} {'β_k(首几个)':>20} {'J':>10}")
    cfgs=[('8家均匀',dist(8,0)),('4家均匀',dist(4,0)),('2家均匀',dist(2,0)),
          ('4家:1大3小',dist(4,0.8)),('1家',[1.0])]
    res=[]
    for name,sh in cfgs:
        J,betas=solve_hetero_bpr(sf,alpha,sh,lam_beta,iters=4000)
        res.append((name,HHI(sh),J))
        print(f"{name:>18} {HHI(sh):>6.3f} {str(np.round(betas[:3],2)):>20} {J:>10.0f}")
    print(f"\n用时={time.time()-t:.0f}s")
    # 检查单调
    res_sorted=sorted(res,key=lambda x:x[1])
    Js=[r[2] for r in res_sorted]
    mono=all(Js[i]>=Js[i+1]-1 for i in range(len(Js)-1))
    print(f"按HHI升序J: {[round(j) for j in Js]}")
    print(f"J随HHI单调↓(ρ单调↑)? {'✓是' if mono else '✗否'}")
