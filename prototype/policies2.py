import numpy as np, itertools
from core import *
def screened_two_step(C,inst,n,k=4):
    """two-step lookahead, but each level only evaluates the k candidates with the most negative first-order scores
    (scores from one adjoint solve at the current / intermediate equilibrium)."""
    p=[[i] for i in range(n)]; eq=C.eq(p); J=eq['J']
    while len(p)>1:
        m=pair_scores(inst,eq)
        c1=sorted((merge_score(m,p[a],p[b]),a,b) for a,b in itertools.combinations(range(len(p)),2))[:k]
        best=(J,None)
        for _,a,b in c1:
            q=merge(p,a,b); eq1=C.eq(q); J1=eq1['J']
            if J1<best[0]*(1-1e-10): best=(J1,q)
            if len(q)>1:
                m1=pair_scores(inst,eq1)
                c2=sorted((merge_score(m1,q[a2],q[b2]),a2,b2) for a2,b2 in itertools.combinations(range(len(q)),2))[:k]
                for _,a2,b2 in c2:
                    J2=C.J(merge(q,a2,b2))
                    if J2<best[0]*(1-1e-10): best=(J2,q)
        if best[1] is None: break
        p=best[1]; eq=C.eq(p); J=eq['J']
    return p

def predict_J(inst, eq, newpart):
    """Theorem-1 style prediction: re-solve the KKT system of the NEW partition on the CURRENT active support
    (one linear solve, no inequality handling). Exact if the support does not change (affine costs)."""
    S=eq['S']; H=inst.hessian(E_of_partition(newpart,inst.n)); A=inst.A; m=A.shape[0]; q=inst.qvec()
    K=np.block([[H[np.ix_(S,S)], -A[:,S].T],[A[:,S], np.zeros((m,m))]])
    sol=np.linalg.lstsq(K,np.concatenate([-q[S],inst.d]),rcond=None)[0]
    y=np.zeros(inst.N); y[S]=sol[:len(S)]
    return inst.J(y)

def predicted_two_step(C,inst,n,k=3):
    """two-step lookahead where candidates are ranked by fixed-support predictions (linear solves),
    and only the top-k at each level are verified with certified equilibrium solves."""
    p=[[i] for i in range(n)]; eq=C.eq(p); J=eq['J']
    while len(p)>1:
        c1=sorted((predict_J(inst,eq,merge(p,a,b)),a,b) for a,b in itertools.combinations(range(len(p)),2))[:k]
        best=(J,None)
        for _,a,b in c1:
            q=merge(p,a,b); eq1=C.eq(q); J1=eq1['J']
            if J1<best[0]*(1-1e-10): best=(J1,q)
            if len(q)>1:
                c2=sorted((predict_J(inst,eq1,merge(q,a2,b2)),a2,b2) for a2,b2 in itertools.combinations(range(len(q)),2))[:k]
                for _,a2,b2 in c2:
                    J2=C.J(merge(q,a2,b2))
                    if J2<best[0]*(1-1e-10): best=(J2,q)
        if best[1] is None: break
        p=best[1]; eq=C.eq(p); J=eq['J']
    return p
