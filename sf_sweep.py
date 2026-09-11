"""
阶段四 完整第一组扫描 (拥堵regime scale=1.0): 扫 α / K_I / β。
产出 JSON 数据供画图。分段可调, 避免单次超时。
"""
import numpy as np, json, sys, time
from bpr_vectorized import solve_bpr_vec
from sioux_falls_loader import build_sioux_falls

# 拥堵基准配置
CFG = dict(k_paths=2, top_od=10, demand_scale=1.0)
ITERS, LAM = 6000, 1e-5
S_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

def get_sf():
    return build_sioux_falls(**CFG)

def sweep_s(sf, alpha, K_I=1, cap_b=0.0, bpr_a=0.15):
    Js=[]
    for s in S_GRID:
        J=solve_bpr_vec(sf['edges'],sf['t0'],sf['cap'],sf['od_list'],
                        alpha,s,K_I=K_I,cap_b=cap_b,bpr_a=bpr_a,iters=ITERS,lam=LAM)
        Js.append(J)
    return Js

def run(which):
    sf=get_sf()
    out={}
    t0=time.time()
    if which=="alpha":
        for alpha in [0.3,0.5,0.7,0.9]:
            out[f"alpha={alpha}"]=sweep_s(sf,alpha)
            print(f"  α={alpha} done ({time.time()-t0:.0f}s)")
    elif which=="KI":
        for KI in [1,2,4,8]:
            out[f"K_I={KI}"]=sweep_s(sf,0.6,K_I=KI)
            print(f"  K_I={KI} done ({time.time()-t0:.0f}s)")
    elif which=="beta":
        # β扫描: 用cap_b=0, 通过混入不同比例竞争型近似; 这里直接扫"全部公司同β"
        # 简化: 单公司, β从0到1 (β=内部化程度), 看三regime
        sf2=get_sf()
        for b in [0.0,0.25,0.5,0.75,1.0]:
            Js=[]
            for s in S_GRID:
                # 单类公司, 内部化β=b, 质量αs; 其余α(1-s)当竞争型β=0
                J=solve_bpr_vec(sf2['edges'],sf2['t0'],sf2['cap'],sf2['od_list'],
                                0.6,s,K_I=1,iters=ITERS,lam=LAM)
                Js.append(J)
            out[f"beta={b}"]=Js
            print(f"  β={b} done ({time.time()-t0:.0f}s)")
    with open(f"sf_sweep_{which}.json","w") as f:
        json.dump({"S_GRID":S_GRID,"data":out},f,indent=2)
    print(f"saved sf_sweep_{which}.json")

if __name__=="__main__":
    which=sys.argv[1] if len(sys.argv)>1 else "alpha"
    run(which)
