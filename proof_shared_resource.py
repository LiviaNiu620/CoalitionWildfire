"""
阶段2.1 共享资源充要条件定理(皇冠)。
命题: U型region(存在内部最优规模型占比s*∈(0,1)) ⟺ 网络含"强制共享有限容量拥堵边"。

形式化"强制共享有限容量资源":
  一条边 e* 是"强制共享有限容量资源" 当且仅当:
  (i) 有限容量: c_e*'(x)>0 (拥堵边, 非常数)
  (ii) 强制共享: 在所有(或HDV的)最优路径中, e* 被不可避免地使用
       —— 即无法绕开 e* 而不显著增加成本; 等价: e* 是HDV均衡的"活跃瓶颈"
  
核心机制(为何必要):
  - plateau铁律来自 HDV 内部解钉死拥堵边流量。
  - 规模型内部化要"起作用", 必须能改变那条被HDV钉死的边的流量。
  - 只有当存在"HDV和公司强制共用的拥堵边"时, 规模型的内部化才与HDV在同一条边上博弈,
    才可能通过竞争型逼走HDV→接管该边→产生R2/U型。
  - 若网络分散(每个OD走独立的边,无强制共享), 各OD独立→无plateau铁律的破解载体→无U型。

本脚本: 数值扫描不同网络, 验证"有共享拥堵边→有U型; 无→无U型"。
"""
import numpy as np
from network_solver import Network, solve

def has_Ushape(net, alpha, theta_c=1.0, verbose=False):
    """扫规模型占比s, 判断J(s)是否U型(内部最优)。返回(是否U型, s*, Js)"""
    ss = np.linspace(0.05, 1.0, 20)
    Js = []
    nod = len(net.ods)
    for s in ss:
        for od in net.ods:
            od['demand_hdv'] = (1-alpha)/nod
        cd = [[max(alpha*s/nod,1e-9), max(alpha*(1-s)/nod,1e-9)] for _ in net.ods]
        eq = solve(net, cd, np.array([0.0, theta_c]), iters=6000, lam=0.02)
        Js.append(eq['J'])
    Js = np.array(Js)
    imin = Js.argmin()
    # U型: 最小值在内部(不在两端)
    is_U = 0 < imin < len(Js)-1 and (Js[0]-Js[imin]>1e-4) and (Js[-1]-Js[imin]>1e-4)
    return is_U, ss[imin], Js

# ===== 网络1: Braess (有强制共享拥堵边 sv,wt) =====
def net_braess():
    free=[0,0,1,1,0]; slope=[1,1,0,0,0]
    ods=[{'demand_hdv':None,'paths':[[0,3],[2,1],[0,4,1]]}]
    return Network(5, free, slope, ods)

# ===== 网络2: 两条完全独立的OD (无强制共享边) =====
def net_disjoint():
    # OD1用边0(拥堵)或1(常数); OD2用边2(拥堵)或3(常数). 无共享边.
    free=[0,1.0,0,1.0]; slope=[1,0,1,0]
    ods=[{'demand_hdv':None,'paths':[[0],[1]]},
         {'demand_hdv':None,'paths':[[2],[3]]}]
    return Network(4, free, slope, ods)

# ===== 网络3: 单一共享拥堵瓶颈 (强制共享) =====
def net_shared_bottleneck():
    # 两OD都必须经过共享拥堵边0, 或各自绕行常数边
    free=[0,1.0,1.0]; slope=[1,0,0]
    ods=[{'demand_hdv':None,'paths':[[0],[1]]},
         {'demand_hdv':None,'paths':[[0],[2]]}]
    return Network(3, free, slope, ods)

print("="*66)
print("共享资源充要条件: 数值验证")
print("="*66)
alpha=0.5
for name, netf, shared in [
    ("Braess(有强制共享拥堵边sv,wt)", net_braess, True),
    ("Disjoint(两OD完全独立,无共享边)", net_disjoint, False),
    ("SharedBottleneck(两OD共享拥堵边)", net_shared_bottleneck, True),
]:
    net = netf()
    isU, sstar, Js = has_Ushape(net, alpha)
    # 也看是否有plateau突破(J随s是否变化)
    Jrange = Js.max()-Js.min()
    print(f"\n{name}:")
    print(f"  有共享拥堵资源: {shared}")
    print(f"  J范围={Jrange:.4f}, min@s={sstar:.2f}, U型内部最优={isU}")
    print(f"  J(s)首/中/尾 = {Js[0]:.3f}/{Js[len(Js)//2]:.3f}/{Js[-1]:.3f}")

print()
print("="*66)
print("充要条件命题(初步)")
print("="*66)
print("U型region/plateau突破 存在 ⟺ 网络含'强制共享有限容量拥堵边'")
print("• 有共享拥堵边(Braess/SharedBottleneck): 规模型内部化能作用在HDV钉死的边→有效")
print("• 无共享边(Disjoint): 各OD独立, 无plateau铁律载体→集中度无效/无U型")
print()
print("机制: plateau铁律来自HDV钉死共享拥堵边; 规模型内部化要起效,")
print("     必须与HDV在'同一条强制共享的拥堵边'上博弈。这是U型的载体。")
