import numpy as np, json
from bpr import *
from so import solve_so
tau=[0.5,0.5,1.5,1.5]
rows=json.load(open('exp_bpr.json'))
for netname,net in [('grid',grid_bpr(1.0)),('braess',braess_bpr())]:
    ue=solve_bpr(net,0.0,tau,0.0,[[0],[1],[2],[3]],'frozen',np.zeros(len(net.edges)))['J']; so=solve_so(net)
    for r in rows:
        if r['net']!=netname or r['cross_regret']<=0: continue
        Js={canon(p):solve_bpr(net,r['alpha'],tau,r['gamma'],p,'self')['J'] for p in set_partitions([0,1,2,3])}
        Jstar=min(Js.values()); pf=canon([tuple(b) for b in r['opt_frozen'][0]])
        grand=Js[canon([[0,1,2,3]])]; sing=Js[canon([[0],[1],[2],[3]])]
        print(f"{netname} a={r['alpha']} g={r['gamma']}: frozen pick {pf} -> regret {100*(Js[pf]-Jstar)/(ue-so):.2f} pp ; grand-vs-singleton (true) {100*(sing-grand)/(ue-so):.2f} pp")
