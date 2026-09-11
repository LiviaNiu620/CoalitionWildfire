import numpy as np, json, time, sys, itertools
from core import *
from nets import grid
from exp_n4 import Counter, one_step, two_step
from policies2 import predicted_two_step, predict_J
n=int(sys.argv[1]); seeds={6:7,8:11}; rng=np.random.default_rng(seeds[n]); T=int(sys.argv[2])
states=[(0.5,0.5),(0.5,1.0),(0.9,0.5),(0.9,1.0)] if n==6 else [(0.5,0.5),(0.9,0.5)]
allrows=[]
for trial in range(T):
    tau=np.round(rng.uniform(0.3,2.0,n),2); port=rng.dirichlet(np.ones(3)*0.7,n); port=port/port.sum(0,keepdims=True)
    for alpha,gamma in states:
        inst=grid(alpha,tau,gamma,portfolio=port); cache={}; y0=equilibrium(inst,[[i] for i in range(n)])['y']
        for p in set_partitions(list(range(n))): cache[canon(p)]=equilibrium(inst,p,y0)
        Jstar=min(v['J'] for v in cache.values()); Js=cache[canon([[i] for i in range(n)])]['J']
        # prediction quality over all single merges from all partitions (sampled)
        errs=[]; agree=[]
        keys=list(cache.keys()); rs=np.random.default_rng(0).choice(len(keys),min(60,len(keys)),replace=False)
        for idx in rs:
            p=[list(b) for b in keys[idx]]; eq=cache[keys[idx]]
            for a,b in itertools.combinations(range(len(p)),2):
                q=merge(p,a,b); Jt=cache[canon(q)]['J']; Jp=predict_J(inst,eq,q)
                if abs(Jt-eq['J'])>1e-9*Jstar: agree.append(np.sign(Jp-eq['J'])==np.sign(Jt-eq['J']))
        res={}
        for pol,fn in [('one_step',lambda C:one_step(C,n)),('two_step',lambda C:two_step(C,n)),
                       ('pred2_k2',lambda C:predicted_two_step(C,inst,n,2)),('pred2_k3',lambda C:predicted_two_step(C,inst,n,3))]:
            C=Counter(inst,cache); p=fn(C); res[pol]=(round((cache[canon(p)]['J']-Jstar)/(Js-Jstar),3),len(C.solved))
        print(trial,alpha,gamma,'pred sign acc %.3f'%np.mean(agree),res,flush=True)
        allrows.append(dict(trial=trial,alpha=alpha,gamma=gamma,sign_acc=float(np.mean(agree)),**{k:v for k,v in res.items()}))
json.dump(allrows,open(f'exp_pred_n{n}.json','w'),indent=1)
