"""n=4 full audit + policies: one-step, two-step, score-guided (adjoint), pairwise surrogate."""
import numpy as np, itertools, json, time, sys
from core import *
from nets import braess, grid
from scipy.stats import spearmanr

def audit(inst):
    n=inst.n; parts=list(set_partitions(list(range(n))))
    res={}; y0=None
    sing=equilibrium(inst,[[i] for i in range(n)]); y0=sing['y']
    for p in parts:
        eq=equilibrium(inst,p,y0); res[canon(p)]=eq
    return res

class Counter:
    def __init__(s,inst,cache): s.inst=inst; s.cache=cache; s.solved=set()
    def J(s,p):
        k=canon(p); s.solved.add(k); return s.cache[k]['J']
    def eq(s,p):
        k=canon(p); s.solved.add(k); return s.cache[k]

def one_step(C,n):
    p=[[i] for i in range(n)]; J=C.J(p)
    while len(p)>1:
        cands=[merge(p,a,b) for a,b in itertools.combinations(range(len(p)),2)]
        Js=[C.J(c) for c in cands]; k=int(np.argmin(Js))
        if Js[k]<J*(1-1e-10): p,J=cands[k],Js[k]
        else: break
    return p

def two_step(C,n):
    p=[[i] for i in range(n)]; J=C.J(p)
    while len(p)>1:
        best=(J,None)
        for a,b in itertools.combinations(range(len(p)),2):
            c1=merge(p,a,b); J1=C.J(c1)
            if J1<best[0]*(1-1e-10): best=(J1,c1)
            for a2,b2 in itertools.combinations(range(len(c1)),2):
                c2=merge(c1,a2,b2); J2=C.J(c2)
                if J2<best[0]*(1-1e-10): best=(J2,c1)   # take first merge of best path
        if best[1] is None: break
        p=best[1]; J=C.J(p)
    return p

def score_greedy(C,inst,n,k=1):
    """rank candidate merges by adjoint first-order score (1 linear solve/step); verify only top-k."""
    p=[[i] for i in range(n)]; eq=C.eq(p); J=eq['J']
    while len(p)>1:
        m=pair_scores(inst,eq)
        cands=[(merge_score(m,p[a],p[b]),a,b) for a,b in itertools.combinations(range(len(p)),2)]
        cands=[c for c in sorted(cands) if c[0]<0][:k]
        best=None
        for sc,a,b in cands:
            c=merge(p,a,b); Jc=C.J(c)
            if Jc<J*(1-1e-10) and (best is None or Jc<best[0]): best=(Jc,c)
        if best is None: break
        J,p=best; eq=C.eq(p)
    return p

def surrogate_cc(C,inst,n):
    """pairwise first-order surrogate at singleton equilibrium -> correlation clustering (exact enum for small n),
    then verify by score-greedy from the surrogate solution."""
    p0=[[i] for i in range(n)]; eq=C.eq(p0); m=pair_scores(inst,eq)
    best=None
    for p in set_partitions(list(range(n))):
        val=sum(m[i,j] for blk in p for i,j in itertools.combinations(blk,2))
        if best is None or val<best[0]: best=(val,p)
    p=best[1]; C.J(p)
    return p

if __name__=="__main__":
    rows=[]
    tau=[0.5,0.5,1.5,1.5]
    for name,mk in [('braess',lambda a,g: braess(a,tau,g)),('grid',lambda a,g: grid(a,tau,g))]:
        for alpha in [0.1,0.3,0.5,0.7,0.9]:
            for gamma in [0.25,0.5,0.75,1.0]:
                inst=mk(alpha,gamma); t=time.time(); cache=audit(inst); ta=time.time()-t
                Jstar=min(v['J'] for v in cache.values())
                # UE / SO benchmarks for normalisation: all-HDV UE and SO on same network
                opt=[k for k,v in cache.items() if v['J']-Jstar<=1e-9*Jstar]
                # sign accuracy of first-order merge score for every partition x candidate merge
                sgn=[];pairs=[]
                for key,eq in cache.items():
                    p=[list(b) for b in key]
                    if len(p)==1: continue
                    m=pair_scores(inst,eq)
                    for a,b in itertools.combinations(range(len(p)),2):
                        sc=merge_score(m,p[a],p[b]); dJ=cache[canon(merge(p,a,b))]['J']-eq['J']
                        if abs(dJ)>1e-9*Jstar: sgn.append(np.sign(sc)==np.sign(dJ)); pairs.append((sc,dJ))
                pairs=np.array(pairs)
                rec={'net':name,'alpha':alpha,'gamma':gamma,'Jstar':Jstar,'optK':sorted({len(k) for k in opt}),
                     'sign_acc':float(np.mean(sgn)) if sgn else None,'n_merges':len(sgn),
                     'spearman':float(spearmanr(pairs[:,0],pairs[:,1])[0]) if len(pairs)>2 else None}
                for pol,fn in [('one_step',lambda C:one_step(C,4)),('two_step',lambda C:two_step(C,4)),
                               ('score_k1',lambda C:score_greedy(C,inst,4,1)),('score_k2',lambda C:score_greedy(C,inst,4,2)),
                               ('surrogate',lambda C:surrogate_cc(C,inst,4)),('grand',lambda C:(C.J([[0,1,2,3]]),[[0,1,2,3]])[1])]:
                    C=Counter(inst,cache); p=fn(C)
                    rec[pol]={'regret_rel':(cache[canon(p)]['J']-Jstar)/Jstar,'solves':len(C.solved),'K':len(p)}
                rows.append(rec)
                print(name,alpha,gamma,'optK',rec['optK'],'signacc %.3f'%rec['sign_acc'] if rec['sign_acc'] is not None else '-',
                      ' '.join(f"{k}:{rec[k]['regret_rel']:.1e}/{rec[k]['solves']}" for k in ['one_step','two_step','score_k1','score_k2','surrogate','grand']),
                      'audit %.1fs'%ta, flush=True)
    json.dump(rows,open('exp_n4.json','w'),indent=1,default=str)
