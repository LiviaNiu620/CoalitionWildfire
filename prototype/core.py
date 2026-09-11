"""Prototype: coalition-structured mixed equilibrium (frozen-slope / affine => exact potential QP),
adjoint merge scores, and partition search. Affine costs c_e = a_e + b_e x_e."""
import numpy as np, itertools
from scipy.optimize import minimize

# ---------------- network utilities ----------------
def simple_paths(edges, o, d, maxlen=8):
    adj = {}
    for k,(u,v) in enumerate(edges): adj.setdefault(u, []).append((v,k))
    out=[]
    def dfs(u, seen, path):
        if len(path)>maxlen: return
        if u==d: out.append(list(path)); return
        for v,k in adj.get(u,[]):
            if v not in seen:
                seen.add(v); path.append(k); dfs(v,seen,path); path.pop(); seen.remove(v)
    dfs(o,{o},[])
    return out

class Instance:
    """n companies + 1 HDV population; proportional OD portfolios unless portfolio matrix given."""
    def __init__(self, edges, a, b, ods, demand, alpha, tau, gamma, bbar=None, portfolio=None, q=None):
        self.E=len(edges); self.a=np.asarray(a,float); self.b=np.asarray(b,float)
        self.bbar = self.b.copy() if bbar is None else np.asarray(bbar,float)
        self.tau=np.asarray(tau,float); self.n=len(tau); self.gamma=gamma
        self.routes=[]; self.route_od=[]
        for w,(o,d) in enumerate(ods):
            for p in simple_paths(edges,o,d): self.routes.append(p); self.route_od.append(w)
        self.R=len(self.routes); self.W=len(ods)
        self.Delta=np.zeros((self.E,self.R))
        for r,p in enumerate(self.routes): self.Delta[p,r]=1
        demand=np.asarray(demand,float)
        n=self.n; P=n+1
        if portfolio is None: portfolio=np.full((n,self.W),1.0/n)   # share of AV demand per company per OD
        self.dem=np.zeros((P,self.W))
        self.dem[:n]=alpha*demand[None,:]*portfolio
        self.dem[n]=(1-alpha)*demand
        # equality constraints: for each population p and od w
        self.N=P*self.R
        rows=[]; rhs=[]
        for p in range(P):
            for w in range(self.W):
                row=np.zeros(self.N)
                idx=[p*self.R+r for r in range(self.R) if self.route_od[r]==w]
                row[idx]=1; rows.append(row); rhs.append(self.dem[p,w])
        self.A=np.array(rows); self.d=np.array(rhs)
        self.q_extra = np.zeros(self.N) if q is None else np.asarray(q,float)  # linear route priorities
        self.Lam_tot=np.hstack([self.Delta]*P)               # E x N
        self.DtBbarD=self.Delta.T@np.diag(self.bbar)@self.Delta
        self.H_road=self.Lam_tot.T@np.diag(self.b)@self.Lam_tot
    def Q(self, E_mat):
        s=np.sqrt(self.tau)
        return (np.outer(s,s))*((1-self.gamma)*np.eye(self.n)+self.gamma*E_mat)
    def hessian(self, E_mat):
        n=self.n; R=self.R
        H=self.H_road.copy()
        Q=self.Q(E_mat)
        for i in range(n):
            for j in range(n):
                H[i*R:(i+1)*R, j*R:(j+1)*R]+=Q[i,j]*self.DtBbarD
        return H
    def qvec(self): return self.Lam_tot.T@self.a + self.q_extra
    def x_of(self,y): return self.Lam_tot@y
    def J(self,y):
        x=self.x_of(y); return float(self.a@x + x@(self.b*x))

def E_of_partition(part, n):
    E=np.zeros((n,n))
    for blk in part:
        for i in blk:
            for j in blk: E[i,j]=1
    return E

# ---------------- convex QP: min 1/2 y'Hy + q'y, Ay=d, y>=0 ----------------
def solve_qp(H,q,A,d,y0=None,tol=1e-10,maxit=200):
    N=len(q)
    scale=max(1.0,np.abs(d).max())
    if y0 is not None:
        try:
            return _active_set(H,q,A,d,set(np.where(y0>1e-7*scale)[0]),tol,maxit)
        except RuntimeError:
            pass
    y0=np.linalg.lstsq(A,d,rcond=None)[0]; y0=np.maximum(y0,0)
    # warm phase: SLSQP
    res=minimize(lambda y:0.5*y@H@y+q@y, y0, jac=lambda y:H@y+q, method='SLSQP',
                 bounds=[(0,None)]*N, constraints=[{'type':'eq','fun':lambda y:A@y-d,'jac':lambda y:A}],
                 options={'ftol':1e-15,'maxiter':2000})
    y=np.maximum(res.x,0)
    return _active_set(H,q,A,d,set(np.where(y>1e-7*scale)[0]),tol,maxit)

def _active_set(H,q,A,d,S,tol,maxit):
    N=len(q); scale=max(1.0,np.abs(d).max())
    for it in range(maxit):
        Sl=sorted(S); m=A.shape[0]
        K=np.block([[H[np.ix_(Sl,Sl)], -A[:,Sl].T],[A[:,Sl], np.zeros((m,m))]])
        sol=np.linalg.lstsq(K, np.concatenate([-q[Sl], d]), rcond=None)[0]
        yS=sol[:len(Sl)]; lam=sol[len(Sl):]
        if yS.min()< -1e-9*scale:
            # drop most negative, step along
            k=Sl[int(np.argmin(yS))]; S.discard(k); continue
        y=np.zeros(N); y[Sl]=np.maximum(yS,0)
        rc=H@y+q-A.T@lam
        viol=[k for k in range(N) if k not in S and rc[k]< -tol*max(1,np.abs(rc).max())]
        if not viol:
            feas=np.abs(A@y-d).max()/scale; stat=np.abs(rc[Sl]).max()/max(1,np.abs(H@y+q).max()) if Sl else 0
            if feas>1e-8 or stat>1e-8: raise RuntimeError('KKT residual too large')
            return y, lam, Sl, rc
        S.add(viol[int(np.argmin(rc[viol]))])
    raise RuntimeError('active set did not converge')

def equilibrium(inst, part, y0=None):
    H=inst.hessian(E_of_partition(part,inst.n)); q=inst.qvec()
    y,lam,S,rc=solve_qp(H,q,inst.A,inst.d,y0)
    return dict(y=y,lam=lam,S=S,rc=rc,H=H,J=inst.J(y))

def vi_gap(inst, eq):
    """relative complementarity residual: sum_r y_r*(rc_r) and min rc on unused."""
    rc=eq['rc']; y=eq['y']
    return float(np.abs(y*rc).sum()/ (np.abs(eq['H']@y+inst.qvec())@y+1e-12)), float(rc.min())

# ---------------- sensitivities / adjoint merge scores ----------------
def dHdeta_pair(inst, i, j):
    """d Hessian / d eta when adding cross-coupling between companies i and j (symmetric)."""
    n=inst.n; R=inst.R
    dQ=inst.gamma*np.sqrt(inst.tau[i]*inst.tau[j])
    G=np.zeros((inst.N,inst.N))
    G[i*R:(i+1)*R, j*R:(j+1)*R]+=dQ*inst.DtBbarD
    G[j*R:(j+1)*R, i*R:(i+1)*R]+=dQ*inst.DtBbarD
    return G

def adjoint(inst, eq):
    S=eq['S']; H=eq['H']; A=inst.A; m=A.shape[0]
    K=np.block([[H[np.ix_(S,S)], -A[:,S].T],[A[:,S], np.zeros((m,m))]])
    x=inst.x_of(eq['y']); c=inst.Lam_tot.T@(inst.a+2*inst.b*x)
    wmu=np.linalg.lstsq(K.T, np.concatenate([c[S], np.zeros(m)]), rcond=None)[0]
    w=np.zeros(inst.N); w[S]=wmu[:len(S)]
    return w

def pair_scores(inst, eq, w=None):
    """m_ij = dJ/deta for adding coupling (i,j) at current equilibrium; ONE adjoint solve for all pairs."""
    if w is None: w=adjoint(inst,eq)
    y=eq['y']; n=inst.n; m=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            g=dHdeta_pair(inst,i,j)@y
            m[i,j]=m[j,i]=-(w@g)
    return m

def merge_score(m, A_blk, B_blk):
    return sum(m[i,j] for i in A_blk for j in B_blk)

def J_along_path(inst, part, A_blk, B_blk, eta, y0=None):
    E=E_of_partition(part,inst.n)
    for i in A_blk:
        for j in B_blk: E[i,j]=E[j,i]=eta
    H=inst.hessian(E); y,lam,S,rc=solve_qp(H,inst.qvec(),inst.A,inst.d,y0)
    return inst.J(y), y

# ---------------- partitions ----------------
def set_partitions(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        yield [[first]]+p
        for k in range(len(p)):
            yield p[:k]+[[first]+p[k]]+p[k+1:]

def canon(part): return tuple(sorted(tuple(sorted(b)) for b in part))
def merge(part, a, b):
    new=[blk for k,blk in enumerate(part) if k not in (a,b)]; new.append(sorted(part[a]+part[b])); return new
