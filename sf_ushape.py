"""
在真实Sioux Falls制造U型: 增强竞争型触发过载R3。
机制: 竞争型β_C < 0 (负内部化=过度抢占共享瓶颈), 或降备选容量放大稀缺。
扫规模型占比s, 看J是否呈U型(两端差,中间好)。
"""
import numpy as np
from sioux_falls_loader import build_sioux_falls

# 复制向量化求解器, 但竞争型β可调(含负值)
def solve(E,t0,cap,od_list,alpha,s,beta_scale=1.0,beta_comp=0.0,
          bpr_a=0.15,bpr_b=4.0,iters=6000,lam=1e-5,cap_scale=1.0):
    cap=cap*cap_scale  # 降容量放大稀缺
    paths=[];path_od=[]
    for oi,od in enumerate(od_list):
        for pe in od['paths']: paths.append(pe);path_od.append(oi)
    P=len(paths);path_od=np.array(path_od)
    A=np.zeros((E,P))
    for p,pe in enumerate(paths):
        for e in pe: A[e,p]=1.0
    n_od=len(od_list)
    # 2类公司: 规模型(β=beta_scale) + 竞争型(β=beta_comp可负)
    betas=np.array([beta_scale, beta_comp])
    dem_H=np.array([(1-alpha)*od['demand'] for od in od_list])
    dem_firm=np.zeros((n_od,2))
    for oi,od in enumerate(od_list):
        d=od['demand']; dem_firm[oi,0]=alpha*s*d; dem_firm[oi,1]=max(alpha*(1-s)*d,1e-9)
    fH=np.zeros(P);fF=np.zeros((2,P))
    for oi in range(n_od):
        m=path_od==oi;npth=m.sum();fH[m]=dem_H[oi]/npth
        for k in range(2):fF[k,m]=dem_firm[oi,k]/npth
    om=[path_od==oi for oi in range(n_od)]
    def bpr(x):return t0*(1+bpr_a*(x/cap)**bpr_b)
    def bprp(x):return t0*bpr_a*bpr_b*(x**(bpr_b-1))/(cap**bpr_b)
    def proj(v,tot):
        if tot<=1e-12 or len(v)==0:return np.zeros_like(v)
        if len(v)==1:return np.array([tot])
        u=np.sort(v)[::-1];css=np.cumsum(u)-tot;idx=np.arange(1,len(u)+1)
        pos=np.nonzero(u-css/idx>0)[0]
        if len(pos)==0:return np.maximum(v,0)
        r=pos[-1];th=css[r]/(r+1.0);return np.maximum(v-th,0.0)
    for it in range(iters):
        f=fH+fF.sum(0);x=A@f;ce=bpr(x);cp=bprp(x);pc=A.T@ce
        for oi in range(n_od):
            m=om[oi];fH[m]=proj(fH[m]-lam*pc[m]*dem_H[oi],dem_H[oi])
        for k in range(2):
            xk=A@fF[k];intern=A.T@(cp*xk);grad=pc+betas[k]*intern
            for oi in range(n_od):
                m=om[oi];mk=dem_firm[oi,k]
                if mk<1e-12:continue
                fF[k,m]=proj(fF[k,m]-lam*grad[m]*mk,mk)
    x=A@(fH+fF.sum(0));return float(np.dot(x,bpr(x)))

if __name__=="__main__":
    sf=build_sioux_falls(k_paths=2,top_od=10,demand_scale=1.0)
    S=[0.0,0.2,0.4,0.6,0.8,1.0]
    print("在真实SF制造U型: 扫规模型占比s, 变竞争型β_C(负=过度抢占)")
    print("="*66)
    for bc in [0.0, -0.5, -1.0, -2.0]:
        Js=[solve(sf['edges'],sf['t0'],sf['cap'],sf['od_list'],0.6,s,
                  beta_comp=bc,iters=5000,lam=1e-5) for s in S]
        Js=np.array(Js);imin=Js.argmin()
        shape="U型!" if 0<imin<len(Js)-1 else ("单调降" if imin==len(Js)-1 else "单调升")
        print(f"  β_C={bc:>5}: minJ@s={S[imin]:.1f} {shape:>7} J={[round(j/1000,1) for j in Js]}")
    print()
    print("再试: 竞争型过度抢占(β_C=-1.5) + 降备选容量(cap_scale=0.5)放大稀缺")
    for cs in [1.0,0.7,0.5]:
        Js=[solve(sf['edges'],sf['t0'],sf['cap'],sf['od_list'],0.6,s,
                  beta_comp=-1.5,cap_scale=cs,iters=5000,lam=1e-5) for s in S]
        Js=np.array(Js);imin=Js.argmin()
        shape="U型!" if 0<imin<len(Js)-1 else ("单调降" if imin==len(Js)-1 else "单调升")
        print(f"  cap_scale={cs}: minJ@s={S[imin]:.1f} {shape:>7} J={[round(j/1000,1) for j in Js]}")
