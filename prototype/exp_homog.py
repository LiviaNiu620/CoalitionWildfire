"""Same random n=6 instances but homogeneous tau=1: is grand always optimal? (affine => frozen = true atomic)"""
import numpy as np, json, time
from core import *
from nets import grid, braess
n=6; rng=np.random.default_rng(7); out=[]
for trial in range(6):
    tau_h=np.round(rng.uniform(0.3,2.0,n),2); port=rng.dirichlet(np.ones(3)*0.7,n); port=port/port.sum(0,keepdims=True)
    for tau,label in [(np.ones(n),'homog'),(tau_h,'hetero')]:
        for alpha,gamma in [(0.5,1.0),(0.9,1.0),(0.9,0.5)]:
            inst=grid(alpha,tau,gamma,portfolio=port); y0=equilibrium(inst,[[i] for i in range(n)])['y']
            Js={canon(p):equilibrium(inst,p,y0)['J'] for p in set_partitions(list(range(n)))}
            Jstar=min(Js.values()); Jg=Js[canon([list(range(n))])]; Jsg=Js[canon([[i] for i in range(n)])]
            optK=sorted({len(k) for k,v in Js.items() if v-Jstar<=1e-9*Jstar})
            out.append((trial,label,alpha,gamma,optK,(Jg-Jstar)/max(Jsg-Jstar,1e-12)))
            print(trial,label,alpha,gamma,'optK',optK,'grand regret frac %.3f'%out[-1][-1],flush=True)
