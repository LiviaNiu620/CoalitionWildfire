import numpy as np, json, time, sys
from core import *
from nets import grid
from exp_n4 import Counter, one_step, two_step
from exp_n8 import surrogate_greedy_cc
from policies2 import screened_two_step
n=8; rng=np.random.default_rng(11); out=[]
for trial in range(3):
    tau=np.round(rng.uniform(0.3,2.0,n),2); port=rng.dirichlet(np.ones(3)*0.7,n); port=port/port.sum(0,keepdims=True)
    for alpha,gamma in [(0.9,1.0),(0.5,1.0)]:
        if alpha==0.9: continue
        inst=grid(alpha,tau,gamma,portfolio=port); cache={}; y0=equilibrium(inst,[[i] for i in range(n)])['y']
        t=time.time()
        for p in set_partitions(list(range(n))): cache[canon(p)]=equilibrium(inst,p,y0)
        Jstar=min(v['J'] for v in cache.values()); Js=cache[canon([[i] for i in range(n)])]['J']
        res={}
        for pol,fn in [('screened2_k3',lambda C:screened_two_step(C,inst,n,3)),('screened2_k5',lambda C:screened_two_step(C,inst,n,5))]:
            C=Counter(inst,cache); p=fn(C); res[pol]=((cache[canon(p)]['J']-Jstar)/(Js-Jstar),len(C.solved))
        print(trial,alpha,gamma,res,'audit %.0fs'%(time.time()-t),flush=True)
