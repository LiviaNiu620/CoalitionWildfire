import numpy as np, time
from core import *
from nets import braess, grid
tau=[0.5,0.5,1.5,1.5]
for name,mk in [('braess',lambda a,g: braess(a,tau,g)),('grid',lambda a,g: grid(a,tau,g))]:
    for alpha in [0.5,0.9]:
        inst=mk(alpha,0.75)
        print(f'== {name} alpha={alpha}: routes={inst.R}, vars={inst.N}')
        part=[[0],[1],[2],[3]]
        t=time.time(); eq=equilibrium(inst,part); print(' solve time %.3fs  J=%.4f  compl-res=%.1e minrc=%.1e'%(time.time()-t,eq['J'],*vi_gap(inst,eq)))
        m=pair_scores(inst,eq)
        # finite-difference check for merge {0}+{2} and block {0,1}+{2,3}
        for A,B in [([0],[2]),([0,1],[2,3]),([1],[3])]:
            h=1e-4
            Jp,_=J_along_path(inst,part,A,B,h,eq['y']); Jm=eq['J']
            fd=(Jp-Jm)/h
            print('  merge',A,B,' adjoint=%.6g  FD=%.6g'%(merge_score(m,A,B),fd))
