"""Random heterogeneous instances, n=6 (203 partitions) exhaustive audit vs policies; counts equilibrium solves."""
import numpy as np, itertools, json, time, sys
from core import *
from nets import grid
from exp_n4 import Counter, one_step, two_step, surrogate_cc
def score_vbs(C,inst,n,k=2):
    """score-ranked merges; verify top-k even if first-order score >= 0 (guards screened/support-change cases)."""
    p=[[i] for i in range(n)]; eq=C.eq(p); J=eq['J']
    while len(p)>1:
        m=pair_scores(inst,eq)
        cands=sorted((merge_score(m,p[a],p[b]),a,b) for a,b in itertools.combinations(range(len(p)),2))[:k]
        best=None
        for sc,a,b in cands:
            c=merge(p,a,b); Jc=C.J(c)
            if Jc<J*(1-1e-10) and (best is None or Jc<best[0]): best=(Jc,c)
        if best is None: break
        J,p=best; eq=C.eq(p)
    return p
def local_refine(C,inst,p):
    """move/merge/split local search with true solves (used after surrogate)."""
    J=C.J(p); improved=True
    while improved:
        improved=False; n=inst.n
        neigh=[]
        for a,b in itertools.combinations(range(len(p)),2): neigh.append(merge(p,a,b))
        for i in range(n):
            src=[k for k,blk in enumerate(p) if i in blk][0]
            for k in range(len(p)+1):
                if k==src: continue
                q=[list(b) for b in p]; q[src].remove(i)
                if k==len(p): q.append([i])
                else: q[k].append(i)
                q=[b for b in q if b]; neigh.append(q)
        for q in neigh:
            Jq=C.J(q)
            if Jq<J*(1-1e-10): p,J=q,Jq; improved=True; break
    return p
if __name__=='__main__':
    n=6; rows=[]; rng=np.random.default_rng(7)
    for trial in range(int(sys.argv[1]) if len(sys.argv)>1 else 6):
        tau=np.round(rng.uniform(0.3,2.0,n),2); port=rng.dirichlet(np.ones(3)*0.7,n)
        port=port/port.sum(0,keepdims=True)   # each OD's AV demand fully split across companies
        for alpha in [0.5,0.9]:
            for gamma in [0.5,1.0]:
                inst=grid(alpha,tau,gamma,scale=1.0,portfolio=port)
                t=time.time(); cache={}; y0=equilibrium(inst,[[i] for i in range(n)])['y']
                for p in set_partitions(list(range(n))): cache[canon(p)]=equilibrium(inst,p,y0)
                ta=time.time()-t; Jstar=min(v['J'] for v in cache.values())
                opt=[k for k,v in cache.items() if v['J']-Jstar<=1e-9*Jstar]
                Jg=cache[canon([list(range(n))])]['J']; Js=cache[canon([[i] for i in range(n)])]['J']
                rec=dict(trial=trial,tau=tau.tolist(),alpha=alpha,gamma=gamma,optK=sorted({len(k) for k in opt}),
                         grand_regret_frac=(Jg-Jstar)/max(Js-Jstar,1e-12), audit_s=ta)
                for pol,fn in [('one_step',lambda C:one_step(C,n)),('two_step',lambda C:two_step(C,n)),
                               ('score_vbs2',lambda C:score_vbs(C,inst,n,2)),
                               ('surrogate+refine',lambda C:local_refine(C,inst,surrogate_cc(C,inst,n)))]:
                    C=Counter(inst,cache); p=fn(C)
                    rec[pol]=dict(regret_frac=(cache[canon(p)]['J']-Jstar)/max(Js-Jstar,1e-12),solves=len(C.solved))
                rows.append(rec)
                print(trial,alpha,gamma,'optK',rec['optK'],'grand %.3f'%rec['grand_regret_frac'],
                      ' '.join(f"{k}:{rec[k]['regret_frac']:.3f}/{rec[k]['solves']}" for k in ['one_step','two_step','score_vbs2','surrogate+refine']),
                      'audit %.0fs'%ta,flush=True)
    json.dump(rows,open('exp_n6.json','w'),indent=1)
