import numpy as np
from core import *
from nets import grid
n=6; rng=np.random.default_rng(7)
for trial in range(6):
    tau=np.round(rng.uniform(0.3,2.0,n),2); port=rng.dirichlet(np.ones(3)*0.7,n); port=port/port.sum(0,keepdims=True)
    for alpha,gamma in [(0.5,0.5),(0.5,1.0),(0.9,0.5),(0.9,1.0)]:
        if (trial,alpha,gamma) not in [(0,0.9,1.0),(2,0.5,1.0),(3,0.9,1.0)]: continue
        inst=grid(alpha,tau,gamma,portfolio=port); y0=equilibrium(inst,[[i] for i in range(n)])['y']
        Js={canon(p):equilibrium(inst,p,y0)['J'] for p in set_partitions(list(range(n)))}
        best=sorted(Js.items(),key=lambda kv:kv[1])[:3]
        print('trial',trial,'alpha',alpha,'gamma',gamma,'tau',tau.tolist())
        for k,v in best: print('   ',[[float(tau[i]) for i in blk] for blk in k], '%.4f'%v)
