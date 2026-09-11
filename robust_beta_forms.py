"""
稳健性1: β(m)换形式. 验证核心结论(HHI-ρ单调, 碎片化坍缩)不依赖具体β(m)形式.
四种β(m): 指数饱和/线性截断/幂/logistic. 全部 β:[0,∞)→[0,1] 单调凹增.
Braess载体, 全AV.
"""
import numpy as np

PATH_EDGES={0:[0,3],1:[2,1],2:[0,4,1]}
CONGESTED={0,1}; CONST={2:45.0,3:45.0,4:0.0}; CAP=100.0

def costs(x): return np.array([x[e]/CAP if e in CONGESTED else CONST[e] for e in range(5)])
def pcost(ce): return np.array([sum(ce[e] for e in PATH_EDGES[p]) for p in range(3)])
def sys_cost(x): ce=costs(x); return sum(x[e]*ce[e] for e in range(5))
def proj(v,tot):
    if tot<=1e-12: return np.zeros_like(v)
    u=np.sort(v)[::-1]; css=np.cumsum(u)-tot; idx=np.arange(1,len(u)+1)
    rho=np.nonzero(u-css/idx>0)[0]
    if len(rho)==0: return np.maximum(v,0)
    r=rho[-1]; return np.maximum(v-css[r]/(r+1.0),0.0)

# 四种β(m)形式, 都标定成在 m=500(参考)处β≈0.5
def beta_forms(name, m):
    if name=='指数饱和': return 1-np.exp(-1.4e-3*m)
    if name=='线性截断': return np.minimum(1.0, m/1000.0)
    if name=='幂饱和':   return 1-1.0/(1+m/700.0)   # m/(m+700)
    if name=='logistic': return 1.0/(1+np.exp(-(m-500)/300.0))*2-1 if m>0 else 0
    return 0

def solve(D, shares, bname, iters=40000, step=2e-6):
    AV=D  # 全AV
    masses=np.array(shares)*AV
    K=len(masses)
    betas=np.array([max(0,beta_forms(bname,m)) for m in masses])
    F=[np.full(3,masses[k]/3) for k in range(K)]
    for it in range(iters):
        x=np.zeros(5)
        for p in range(3):
            tot=sum(F[k][p] for k in range(K))
            for e in PATH_EDGES[p]: x[e]+=tot
        ce=costs(x); pc=pcost(ce); cp=np.array([1/CAP if e in CONGESTED else 0 for e in range(5)])
        for k in range(K):
            fke=np.zeros(5)
            for p in range(3):
                for e in PATH_EDGES[p]: fke[e]+=F[k][p]
            intern=np.array([betas[k]*sum(fke[e]/CAP for e in PATH_EDGES[p] if e in CONGESTED) for p in range(3)])
            F[k]=proj(F[k]-step*(pc+intern)*masses[k]*3,masses[k])
    x=np.zeros(5)
    for p in range(3):
        tot=sum(F[k][p] for k in range(K))
        for e in PATH_EDGES[p]: x[e]+=tot
    return sys_cost(x),betas

def wardrop(D,iters=40000,step=2e-6):
    f=np.full(3,D/3)
    for it in range(iters):
        x=np.zeros(5)
        for p in range(3):
            for e in PATH_EDGES[p]: x[e]+=f[p]
        f=proj(f-step*pcost(costs(x))*D*3,D)
    x=np.zeros(5)
    for p in range(3):
        for e in PATH_EDGES[p]: x[e]+=f[p]
    return sys_cost(x)
def social(D,ng=121):
    best=None
    for a in np.linspace(0,D,ng):
        for b in np.linspace(0,D-a,ng):
            c=D-a-b
            if c<-1e-9: continue
            x=np.zeros(5)
            for p,fp in zip(range(3),[a,b,c]):
                for e in PATH_EDGES[p]: x[e]+=fp
            J=sys_cost(x)
            if best is None or J<best: best=J
    return best
def HHI(sh): s=np.array(sh); return (s**2).sum()
def dist(K,c):
    if K==1: return [1.0]
    big=1.0/K+c*(1-1.0/K); return [big]+[(1-big)/(K-1)]*(K-1)

if __name__=="__main__":
    D=4000.0
    Jw=wardrop(D); Jso=social(D)
    def rho(J): return (Jw-J)/(Jw-Jso) if Jw-Jso>1e-6 else 0
    print(f"稳健性1: β(m)换形式 (Braess, D={D}, J_UE={Jw:.0f}, J_SO={Jso:.0f})")
    print("验证HHI-ρ单调是否不依赖β形式:")
    print(f"{'配置':>14} {'HHI':>6} |"+"".join(f'{n:>11}' for n in ['指数饱和','线性截断','幂饱和']))
    print("-"*66)
    cfgs=[('8家均匀',dist(8,0)),('4家均匀',dist(4,0)),('2家均匀',dist(2,0)),('4家:1大3小',dist(4,0.8)),('1家',[1.0])]
    for name,sh in cfgs:
        row=f"{name:>14} {HHI(sh):>6.3f} |"
        for bn in ['指数饱和','线性截断','幂饱和']:
            J,b=solve(D,sh,bn); row+=f"{rho(J):>10.1%} "
        print(row)
    print()
    print("★ 若三种形式下ρ都随HHI单调↑ → 核心结论不依赖β(m)具体形式 ✓")
