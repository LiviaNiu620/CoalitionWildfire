import numpy as np, json, time
from bpr import *
tau=[0.5,0.5,1.5,1.5]
rows=[]
for netname,net in [('grid',grid_bpr(1.0)),('braess',braess_bpr())]:
    ue=solve_bpr(net,0.0,tau,0.0,[[0],[1],[2],[3]],'frozen',np.zeros(len(net.edges)))
    bbar_ue=net.dc(ue['x'])
    for alpha in [0.5,0.7,0.9]:
        for gamma in [0.5,1.0]:
            out={}
            for mode in ['frozen','self']:
                t=time.time(); Js={}; y0=None; gaps=[]
                for p in set_partitions([0,1,2,3]):
                    r=solve_bpr(net,alpha,tau,gamma,p,mode,bbar_ue,y0); Js[canon(p)]=r['J']; gaps.append(r['gap_model'])
                    if mode=='frozen': y0=None
                Jstar=min(Js.values()); opt=[k for k,v in Js.items() if v-Jstar<=1e-9*Jstar]
                out[mode]=dict(Js=Js,opt=opt,Jstar=Jstar,maxgap=max(gaps),t=time.time()-t)
            # cross-regret: partition chosen under frozen, evaluated under self-consistent (true) behavior
            pf=out['frozen']['opt'][0]; Jtrue=out['self']['Js'][pf]
            reg=(Jtrue-out['self']['Jstar'])/out['self']['Jstar']
            kf=sorted({len(k) for k in out['frozen']['opt']}); ks=sorted({len(k) for k in out['self']['opt']})
            # slope mismatch at frozen grand-coalition equilibrium
            rg=solve_bpr(net,alpha,tau,gamma,[[0,1,2,3]],'frozen',bbar_ue)
            used=rg['x']>1e-6*rg['x'].max(); ratio=net.dc(rg['x'])[used]/np.maximum(bbar_ue[used],1e-12)
            print(f"{netname} a={alpha} g={gamma}: optK frozen={kf} self={ks} | frozen-choice regret under true={reg:.2e} | "
                  f"c'(x*)/bbar on used links: median {np.median(ratio):.2f} min {ratio.min():.2f} max {ratio.max():.2f} | maxgap {out['frozen']['maxgap']:.1e}/{out['self']['maxgap']:.1e}",flush=True)
            rows.append(dict(net=netname,alpha=alpha,gamma=gamma,optK_frozen=kf,optK_self=ks,opt_frozen=[list(map(list,k)) for k in out['frozen']['opt']],
                             opt_self=[list(map(list,k)) for k in out['self']['opt']],cross_regret=reg))
json.dump(rows,open('exp_bpr.json','w'),indent=1)
