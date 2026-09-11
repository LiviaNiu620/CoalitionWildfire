import numpy as np
from core import Instance
def braess(alpha, tau, gamma, D=4000.0, **kw):
    edges=[(0,1),(1,3),(0,2),(2,3),(1,2)]   # sv, vt, sw, wt, vw
    a=[0,45,45,0,0]; b=[1/100,0,0,1/100,0]
    return Instance(edges,a,b,[(0,3)],[D],alpha,tau,gamma,**kw)
def grid(alpha, tau, gamma, scale=1.0, seed=3, **kw):
    rng=np.random.default_rng(seed)
    edges=[]
    for r in range(3):
        for c in range(3):
            u=3*r+c
            if c<2: edges.append((u,u+1))
            if r<2: edges.append((u,u+3))
    edges += [(4,2),(6,4),(4,0+1)]  # a few back/diagonal-ish links for route diversity
    E=len(edges)
    a=rng.uniform(5,15,E); cap=rng.uniform(300,700,E)
    b=a/cap
    ods=[(0,8),(0,5),(3,8)]; dem=np.array([1200,700,700])*scale
    return Instance(edges,a,b,ods,dem,alpha,tau,gamma,**kw)
