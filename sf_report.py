"""
DEPRECATED -- retained for provenance only.  Do not use for paper numbers.

This script produced the Sioux Falls table in earlier drafts.  Three defects
were found and are the reason it was replaced by ``sf_equilibrium.py`` +
``sf_recompute.py``:

1. **Non-convergence.**  The loop runs a fixed 6000 iterations of a fixed-step
   simultaneous projected-gradient scheme with no stopping test.  It does not
   reach an equilibrium: at 6000 iterations with four equal operators it
   returns J = 8,044,810, which is *worse* than the user equilibrium
   (7,480,225) even though the operators internalise.  The J values reported in
   the earlier draft (186,128-190,891) are below the physical lower bound
   3,176,000 = sum_w d_w * (free-flow shortest-path time) and cannot be
   produced by this script at all.

2. **Wrong lambda_beta.**  ``lam_beta = 3e-4`` gives lambda_beta * alpha * N =
   97.4, under which every operator has beta ~ 1 and Sigma degenerates to K.
   The Sigma column reported in the draft (14.0 / 4.9 / 2.1 / 7.1 / 1.0)
   corresponds to lambda_beta * alpha * N = 6.78, i.e. lambda_beta = 2.089e-5.
   The script therefore never reproduced its own table.

3. **Incomplete VI residual.**  The ``gap`` computed below loops only over the
   HDV flows ``fH`` and the plain path costs ``pc``; the fleets' perceived
   costs are never tested.  It is a Wardrop residual for the human population,
   not an equilibrium certificate for the game.

A fourth issue is structural rather than a bug: the route set is frozen at the
three shortest free-flow paths per OD pair, which is not an equilibrium route
set for Sioux Falls (the largest profitable omitted-route improvement is 24.1
against a cost scale of 86.8, on 101 of 528 OD pairs).

Use instead::

    python3 sf_recompute.py      # certified table
    python3 sf_equilibrium.py    # solver and column generation

"""
import numpy as np
from sioux_falls_loader import build_sioux_falls
from sf_hetero_verify import solve_hetero_bpr, dist

def solve_with_edgeflow(sf, alpha, shares, lam_beta, iters=6000, seed=0, lam=1e-5):
    """复刻 solve_hetero_bpr 但返回边流x和VI残差, 支持随机初值."""
    E=sf['edges']; t0=sf['t0']; cap=sf['cap']; od_list=sf['od_list']
    paths=[]; path_od=[]
    for oi,od in enumerate(od_list):
        for pe in od['paths']:
            paths.append(pe); path_od.append(oi)
    P=len(paths); path_od=np.array(path_od)
    A=np.zeros((E,P))
    for p,pe in enumerate(paths):
        for e in pe: A[e,p]=1.0
    n_od=len(od_list); Kf=len(shares)
    total_dem=sum(od['demand'] for od in od_list)
    masses=np.array(shares)*alpha*total_dem
    betas=1-np.exp(-lam_beta*masses)
    dem_H=np.array([(1-alpha)*od['demand'] for od in od_list])
    dem_firm=np.zeros((n_od,Kf))
    for oi,od in enumerate(od_list):
        for k in range(Kf): dem_firm[oi,k]=max(alpha*shares[k]*od['demand'],1e-12)
    rng=np.random.default_rng(seed)
    fH=np.zeros(P); fF=np.zeros((Kf,P))
    for oi in range(n_od):
        m=path_od==oi; npth=int(m.sum())
        # 随机初值(seed>0)或均匀(seed=0)
        if seed==0:
            fH[m]=dem_H[oi]/npth
            for k in range(Kf): fF[k,m]=dem_firm[oi,k]/npth
        else:
            w=rng.random(npth); w/=w.sum(); fH[m]=dem_H[oi]*w
            for k in range(Kf):
                w=rng.random(npth); w/=w.sum(); fF[k,m]=dem_firm[oi,k]*w
    od_masks=[(path_od==oi) for oi in range(n_od)]
    def bpr(x): return t0*(1+0.15*(x/cap)**4)
    def bpr_p(x): return t0*0.15*4*(x**3)/(cap**4)
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
    # VI残差(gap): 每个玩家每OD, 感知路径成本的 max-min 应->0 (只在used paths上)
    ce=bpr(x); cp=bpr_p(x); pc=A.T@ce
    gap=0.0
    for oi in range(n_od):
        m=od_masks[oi]
        used=fH[m]>1e-8*dem_H[oi] if dem_H[oi]>0 else np.zeros(m.sum(),bool)
        if used.any():
            pcm=pc[m][used]; gap=max(gap,pcm.max()-pcm.min())
    J=float(np.dot(x,bpr(x)))
    # For BPR exponent 4, c''/c'=3/x on positive-flow edges.
    tau_max=0.0
    min_sym_eigenvalue=float("inf")
    tau_by_operator=[]
    positive=x>1e-12
    for k in range(Kf):
        xk=A@fF[k]
        tau=np.zeros_like(x)
        tau[positive]=betas[k]*xk[positive]*3.0/x[positive]
        tau_max=max(tau_max,float(tau.max()))
        tau_by_operator.append(tau)
    for e in np.where(positive)[0]:
        one=np.ones(Kf+1)
        sym=np.outer(one,one)
        sym[:Kf,:Kf]+=np.diag(betas)
        for k in range(Kf):
            ek=np.zeros(Kf+1); ek[k]=1.0
            sym+=0.5*tau_by_operator[k][e]*(np.outer(ek,one)+np.outer(one,ek))
        min_sym_eigenvalue=min(min_sym_eigenvalue,float(np.linalg.eigvalsh(sym)[0]))
    return J,x,betas,gap,masses,tau_max,min_sym_eigenvalue

if __name__=="__main__":
    sf=build_sioux_falls(k_paths=3, top_od=None, demand_scale=1.0)
    total_dem=sum(od['demand'] for od in sf['od_list'])
    n_paths=sum(len(od['paths']) for od in sf['od_list'])
    lam_beta=3e-4; alpha=0.9
    print(f"# Sioux Falls 参数报告")
    print(f"links={sf['edges']}, OD_pairs={len(sf['od_list'])}, paths(k=3/OD)={n_paths}")
    print(f"total_demand={total_dem:.0f}, demand_scale=1.0, α={alpha}, αN(AV mass)={alpha*total_dem:.0f}")
    print(f"λ_β={lam_beta}, BPR(a=0.15,b=4)")
    print()
    configs=[("8 equal",[1/8]*8),("4 equal",[1/4]*4),("2 equal",[1/2]*2),
             ("1 large+3 small",[0.7,0.1,0.1,0.1]),("monopoly",[1.0])]
    # 基准
    J_UE,_,_,_,_,_,_=solve_with_edgeflow(sf,alpha,[1e6*[1.0][0]] if False else [1.0],lam_beta,iters=1)  # placeholder
    print(f"{'config':<18}{'Σ=Σ1/β':>10}{'J':>12}{'VI gap':>10}{'β_max':>8}")
    for nm,sh in configs:
        J,x,betas,gap,masses,tau,mineig=solve_with_edgeflow(sf,alpha,sh,lam_beta,iters=6000,seed=0)
        Sig=float(np.sum(1/betas))
        print(f"{nm:<18}{Sig:>10.3f}{J:>12.0f}{gap:>10.2e}{betas.max():>8.3f} tau={tau:.3f} mineig={mineig:.3f}")
    # M4: 多初值边流偏差 (以 4 equal 为例)
    print("\n# M4 多初值边流唯一性验证 (4 equal, 5个随机初值 vs 均匀初值):")
    J0,x0,_,_,_,_,_=solve_with_edgeflow(sf,alpha,[1/4]*4,lam_beta,iters=8000,seed=0)
    maxdev=0
    for s in range(1,6):
        Js,xs,_,_,_,_,_=solve_with_edgeflow(sf,alpha,[1/4]*4,lam_beta,iters=8000,seed=s)
        dev=np.max(np.abs(xs-x0))/max(np.max(x0),1e-9)
        maxdev=max(maxdev,dev)
        print(f"  seed={s}: max|Δx|/max(x)={dev:.2e}, ΔJ/J={abs(Js-J0)/J0:.2e}")
    print(f"  => 最大边流相对偏差={maxdev:.2e} (edge flow 数值唯一)")
