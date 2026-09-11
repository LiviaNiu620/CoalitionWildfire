import numpy as np, json, time, sys
from core import *
from nets import grid
from exp_n4 import Counter, one_step, two_step, surrogate_cc
from exp_n6 import score_vbs, local_refine
import itertools
def surrogate_greedy_cc(C,inst,n):
    """correlation clustering on pairwise first-order scores via greedy agglomeration on the SURROGATE (no solves)."""
    eq=C.eq([[i] for i in range(n)]); m=pair_scores(inst,eq)
    p=[[i] for i in range(n)]
    while True:
        best=(0,None)
        for a,b in itertools.combinations(range(len(p)),2):
            v=merge_score(m,p[a],p[b])
            if v<best[0]: best=(v,(a,b))
        if best[1] is None: break
        p=merge(p,*best[1])
    return local_refine(C,inst,p)
if __name__=="__main__":
    n=8; rows=[]; rng=np.random.default_rng(11)
    T=int(sys.argv[1]) if len(sys.argv)>1 else 3
    for trial in range(T):
        tau=np.round(rng.uniform(0.3,2.0,n),2); port=rng.dirichlet(np.ones(3)*0.7,n); port=port/port.sum(0,keepdims=True)
        for alpha,gamma in [(0.9,1.0),(0.5,1.0)]:
            inst=grid(alpha,tau,gamma,portfolio=port)
            t=time.time(); cache={}; y0=equilibrium(inst,[[i] for i in range(n)])['y']
            for p in set_partitions(list(range(n))): cache[canon(p)]=equilibrium(inst,p,y0)
            ta=time.time()-t; Jstar=min(v['J'] for v in cache.values())
            Jg=cache[canon([list(range(n))])]['J']; Js=cache[canon([[i] for i in range(n)])]['J']
            opt=[k for k,v in cache.items() if v['J']-Jstar<=1e-9*Jstar]
            rec=dict(trial=trial,alpha=alpha,gamma=gamma,optK=sorted({len(k) for k in opt}),grand=(Jg-Jstar)/(Js-Jstar),audit_s=ta,n_part=len(cache))
            for pol,fn in [('one_step',lambda C:one_step(C,n)),('two_step',lambda C:two_step(C,n)),
                           ('surrogate+refine',lambda C:surrogate_greedy_cc(C,inst,n))]:
                C=Counter(inst,cache); tt=time.time(); p=fn(C)
                rec[pol]=dict(regret=(cache[canon(p)]['J']-Jstar)/(Js-Jstar),solves=len(C.solved))
            rows.append(rec)
            print(trial,alpha,gamma,'optK',rec['optK'],'grand %.3f'%rec['grand'],
                  ' '.join(f"{k}:{rec[k]['regret']:.3f}/{rec[k]['solves']}" for k in ['one_step','two_step','surrogate+refine']),
                  'audit %d partitions %.0fs'%(len(cache),ta),flush=True)
    json.dump(rows,open('exp_n8.json','w'),indent=1)
