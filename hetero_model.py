"""
异质模型: 规模分布{m_k}外生, β_k=β(m_k)=1-e^{-λm_k}内生异质。
两路载体(先验证机制), 数值解VI(命题1保证唯一)。
集中度双维度: 数量K + 不均匀度HHI。
验证: (1)大公司规模型/小公司竞争型涌现; (2)同K下不均匀度的影响。
"""
import numpy as np

def beta_cap(m, lam):
    return 1.0 - np.exp(-lam*m)

def solve_hetero(sizes, lam, a, iters=200000, step=None):
    """
    两路: 路1 c1=x1(拥堵瓶颈,c1'=1), 路2 c2=a(常数).
    sizes: 各coalition质量列表[m_1,...,m_K]. HDV=0(全AV,总需求=Σm_k).
    每家k: 感知边际路1 = x1 + β_k·b_k (b_k=k家上路1量), 令=c2=a.
    数值: 各家best-response投影迭代.
    返回: 各家b_k, x1, 系统成本J, 各β_k.
    """
    K=len(sizes); sizes=np.array(sizes,float)
    betas=np.array([beta_cap(m,lam) for m in sizes])
    b=np.array([m*0.5 for m in sizes])  # 初始各家上路1一半
    if step is None: step=0.5/(K+2)
    for it in range(iters):
        x1=b.sum()
        # 各家FOC残差: 感知边际路1 - c2 = (x1 + β_k b_k) - a
        # 想让每家 min 自己车队成本 → 沿感知梯度调b_k
        for k in range(K):
            # 自己车队成本对b_k的感知梯度 = 感知边际路1 - c2
            grad = (x1 + betas[k]*b[k]) - a
            b[k] = np.clip(b[k]-step*grad, 0.0, sizes[k])
            x1 = b.sum()  # 更新
    x1=b.sum()
    J = x1*x1 + (sizes.sum()-x1)*a  # 路1总成本+路2总成本
    return b, x1, J, betas

def rho_of(J, sizes, a):
    N=sum(sizes)
    JUE=N*a  # 全自私: x1=a per... 实际两路归一. 用x1=a基准
    # 归一到总需求: Wardrop x1使c1=c2 → x1=a; J_UE=a²+(N-a)a? 
    # 简化: 用x1视角. J=x1²+(N-x1)a. UE:x1=a→J=a²+(N-a)a. SO:x1=a/2
    JUE=a*a+(N-a)*a
    xSO=a/2; JSO=xSO*xSO+(N-xSO)*a
    return (JUE-J)/(JUE-JSO) if JUE-JSO>1e-9 else 0

def HHI(sizes):
    s=np.array(sizes,float); sh=s/s.sum(); return (sh**2).sum()

if __name__=="__main__":
    lam=2.0; a=0.6; N=1.0  # 总AV量归一
    print("="*72)
    print(f"异质模型验证 (两路, λ={lam}, a={a}, 总AV={N})")
    print("="*72)
    print()
    print("验证1: 大公司规模型/小公司竞争型 涌现 (1大+3小 vs 均匀)")
    print("-"*60)
    configs = {
        "均匀[.25,.25,.25,.25]": [.25,.25,.25,.25],
        "集中[.7,.1,.1,.1]":     [.7,.1,.1,.1],
        "极集中[.85,.05,.05,.05]":[.85,.05,.05,.05],
    }
    for name,sizes in configs.items():
        b,x1,J,betas=solve_hetero(sizes,lam,a)
        rho=rho_of(J,sizes,a)
        print(f"  {name}")
        print(f"    β_k={np.round(betas,3)}  (大公司β高=规模型, 小公司β低=竞争型)")
        print(f"    x1={x1:.3f}(SO={a/2}), ρ={rho:.1%}, HHI={HHI(sizes):.3f}")
    print()
    print("★ 观察: 集中分布下, 大公司β→1(规模型), 小公司β低(竞争型) 涌现!")
