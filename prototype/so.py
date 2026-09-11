import numpy as np, json
from bpr import *
from core import solve_qp
def ddc(net,x): return net.t0*1.8*x**2/net.cap**4
def solve_so(net):
    inst=make_inst(net,0.0,[1.0],0.0,np.zeros(len(net.edges)))  # 1 dummy company with zero demand + HDV
    L=inst.Lam_tot; y=np.zeros(inst.N)
    for k,row in enumerate(inst.A):
        idx=np.where(row>0)[0]; y[idx]=inst.d[k]/len(idx)
    for it in range(200):
        x=L@y; g=L.T@(net.c(x)+x*net.dc(x)); Hh=L.T@np.diag(2*net.dc(x)+x*ddc(net,x))@L
        yn,_,_,_=solve_qp(Hh+1e-12*np.eye(inst.N),g-Hh@y,inst.A,inst.d,y)
        f=lambda z: float((L@z)@net.c(L@z)); s=1.0
        while s>1e-8 and f(y+s*(yn-y))>f(y): s*=0.5
        dy=s*(yn-y); y=y+dy
        if np.abs(dy).max()<1e-11*np.abs(y).max(): break
    x=L@y; return float(x@net.c(x))
if __name__=="__main__":
    tau=[0.5,0.5,1.5,1.5]
    rows=json.load(open('exp_bpr.json'))
    for netname,net in [('grid',grid_bpr(1.0)),('braess',braess_bpr())]:
        ue=solve_bpr(net,0.0,tau,0.0,[[0],[1],[2],[3]],'frozen',np.zeros(len(net.edges)))['J']; so=solve_so(net)
        print(netname,'J_UE=%.1f J_SO=%.1f gap=%.2f%%'%(ue,so,100*(ue-so)/ue))
        for r in rows:
            if r['net']==netname and r['cross_regret']>0:
                # cross_regret is relative to J*_self; convert to recovery pp approx using J*_self ~ J
                print('   a=%.1f g=%.1f cross-regret = %.2f recovery pp'%(r['alpha'],r['gamma'], 100*r['cross_regret']*ue/(ue-so)))
