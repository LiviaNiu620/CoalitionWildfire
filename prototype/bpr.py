"""BPR version: frozen reference slope (paper's model) vs self-consistent slope (= true atomic internalization).
Equilibria solved by sequential QP (Newton) on the convex potential; self-consistent mode updates bbar=c'(x) each iterate."""
import numpy as np, itertools
from core import Instance, E_of_partition, solve_qp, set_partitions, canon

class BPRNet:
    def __init__(self, edges, t0, cap, ods, demand):
        self.edges=edges; self.t0=np.asarray(t0,float); self.cap=np.asarray(cap,float); self.ods=ods; self.demand=np.asarray(demand,float)
    def c(self,x):  return self.t0*(1+0.15*(x/self.cap)**4)
    def dc(self,x): return self.t0*0.6*x**3/self.cap**4
    def cint(self,x): return self.t0*(x+0.03*x**5/self.cap**4)

def make_inst(net, alpha, tau, gamma, bbar):
    # affine placeholders a,b are overwritten inside Newton; Instance used for structure only
    return Instance(net.edges, np.zeros(len(net.edges)), np.zeros(len(net.edges)), net.ods, net.demand, alpha, tau, gamma, bbar=bbar)

def potential_parts(inst, net, E_mat):
    n=inst.n; R=inst.R
    Q=inst.Q(E_mat)
    Hman=np.zeros((inst.N,inst.N))
    for i in range(n):
        for j in range(n):
            Hman[i*R:(i+1)*R, j*R:(j+1)*R]=Q[i,j]*(inst.Delta.T@np.diag(inst.bbar)@inst.Delta)
    return Hman

def solve_bpr(net, alpha, tau, gamma, part, mode='frozen', bbar_ref=None, y0=None, tol=1e-11, maxit=400, omega=0.5):
    inst=make_inst(net, alpha, tau, gamma, bbar_ref if bbar_ref is not None else np.zeros(len(net.edges)))
    E_mat=E_of_partition(part, inst.n); L=inst.Lam_tot
    y=y0
    if y is None:
        y=np.zeros(inst.N)
        for k,row in enumerate(inst.A):   # spread each population-OD demand evenly
            idx=np.where(row>0)[0]; y[idx]=inst.d[k]/len(idx)
    for it in range(maxit):
        x=L@y
        if mode=='self': inst.bbar=net.dc(x) if it==0 else (1-omega)*inst.bbar+omega*net.dc(x)
        Hman=potential_parts(inst,net,E_mat)
        grad=L.T@net.c(x)+Hman@y
        Hk=L.T@np.diag(net.dc(x))@L+Hman
        qk=grad-Hk@y
        ynew,lam,S,rc=solve_qp(Hk+1e-12*np.eye(inst.N),qk,inst.A,inst.d,y)
        # damped step on the (frozen-at-this-iterate) potential
        def Pfun(z): return float(net.cint(L@z).sum()+0.5*z@Hman@z)
        step=1.0; P0=Pfun(y)
        while step>1e-6 and Pfun(y+step*(ynew-y))>P0+1e-12*abs(P0): step*=0.5
        dy=step*(ynew-y); y=y+dy
        if np.abs(dy).max()<tol*max(1,np.abs(y).max()) and (mode!='self' or np.abs(inst.bbar-net.dc(L@y)).max()<1e-10*max(1e-12,np.abs(inst.bbar).max())): break
    x=L@y
    if mode=='self': inst.bbar=net.dc(x)
    Hman=potential_parts(inst,net,E_mat)
    F=L.T@net.c(x)+Hman@y
    # VI residual of the TRUE state-dependent model (bbar=c'(x)) -- reported for both modes
    inst_true_bbar=net.dc(x); inst2=make_inst(net,alpha,tau,gamma,inst_true_bbar)
    Ftrue=L.T@net.c(x)+potential_parts(inst2,net,E_mat)@y
    def gap(Fv):   # normalised VI gap: sum_p,w [F.y - d*min_r F]
        g=0.0
        for k,row in enumerate(inst.A):
            idx=np.where(row>0)[0]; g+=Fv[idx]@y[idx]-inst.d[k]*Fv[idx].min()
        return g/(Fv@y)
    return dict(y=y,x=x,J=float(x@net.c(x)),gap_model=gap(F),gap_true=gap(Ftrue),it=it)

def grid_bpr(scale=1.0, seed=3):
    rng=np.random.default_rng(seed)
    edges=[]
    for r in range(3):
        for c in range(3):
            u=3*r+c
            if c<2: edges.append((u,u+1))
            if r<2: edges.append((u,u+3))
    edges += [(4,2),(6,4),(4,1)]
    E=len(edges); t0=rng.uniform(5,15,E); cap=rng.uniform(300,700,E)
    return BPRNet(edges,t0,cap,[(0,8),(0,5),(3,8)],np.array([1200,700,700])*scale)

def braess_bpr(D=4000.0):
    edges=[(0,1),(1,3),(0,2),(2,3),(1,2)]
    return BPRNet(edges,[20,45,45,20,1],[2000,1e5,1e5,2000,1e5],[(0,3)],[D])
