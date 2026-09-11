"""
推导: 公司碎片化 (m_k→0) 时 regime 坍缩。
两路载体: 拥堵路 c1=x1, 备选路 c2=a.
规模型拆成 K_I 家(每家 m=αs/K_I), 竞争型拆成 K_C 家(每家 α(1-s)/K_C).

核心: 原子玩家的内部化项 ∝ 单个公司规模 m_k, 不是总质量.
     一家规模型公司 k 的感知边际(拥堵路) = c1 + m_k * c1'  (只内部化自己这家的流量)
     K_I 家对称 → 每家 m = αs/K_I.
关键区别: 
  - K_I=1: 内部化项 = αs (全部规模型质量)
  - K_I→∞: 每家 m→0, 内部化项→0 → 退化为非原子自私(等同HDV)
"""
import sympy as sp

a, alpha, s, gamma = sp.symbols('a alpha s gamma', positive=True)
KI, KC = sp.symbols('K_I K_C', positive=True, integer=True)

# 单个规模型公司质量
m_I1 = alpha*s/KI     # 每家规模型
m_C1 = alpha*(1-s)/KC # 每家竞争型

print("="*60)
print("规模型每家质量 m_I1 =", m_I1)
print("竞争型每家质量 m_C1 =", m_C1)
print()

# ---- Regime 2 (HDV撤出) 的重推, 含 K_I ----
# 对称: 每家规模型走拥堵路量 uI1, 走备选 vI1; 每家竞争型 uC1.
# 拥堵路总流量 x1 = K_I*uI1 + K_C*uC1  (HDV撤出)
# 一家规模型内部解: 感知拥堵路 = 备选路 a
#   c1 + m_I1 * c1' = a  →  x1 + m_I1 = a   (c1=x1, c1'=1)
#   注意: 内部化项是 m_I1(自己这家), 不是总规模型质量!
# 一家竞争型内部解(spite): x1 + m_C1 - γ*(对手在拥堵路的量) = a
#   对手 = 其他所有公司在拥堵路的量 ≈ x1 - uC1 (自己那份除外), 近似 x1
uI1, uC1 = sp.symbols('u_I1 u_C1', nonnegative=True)
x1 = KI*uI1 + KC*uC1

# 规模型一阶(每家): x1 + m_I1 = a
eqI = sp.Eq(x1 + m_I1, a)
# 竞争型一阶(每家, spite抬对手): x1 + m_C1 - gamma*(x1 - uC1) = a
eqC = sp.Eq(x1 + m_C1 - gamma*(x1 - uC1), a)

sol = sp.solve([eqI, eqC], [uI1, uC1], dict=True)
if sol:
    sB = sol[0]
    uI1s = sp.simplify(sB[uI1]); uC1s = sp.simplify(sB[uC1])
    x1s = sp.simplify((KI*uI1s + KC*uC1s))
    print("Regime2 (含 K_I,K_C):")
    print("  x1 =", x1s)
    J = sp.simplify(x1s**2 + (1-x1s)*a)
    print("  J =", J)
    print()
    # ---- 关键: K_I→∞ 极限 ----
    print("=== K_I → ∞ (规模型碎片化) 极限 ===")
    x1_lim = sp.limit(x1s, KI, sp.oo)
    J_lim = sp.limit(J, KI, sp.oo)
    print("  x1(K_I→∞) =", sp.simplify(x1_lim))
    print("  J(K_I→∞)  =", sp.simplify(J_lim))
    print()
    # 对比 K_I=1
    print("=== K_I = 1 (规模型集中) ===")
    x1_1 = sp.simplify(x1s.subs(KI,1))
    J_1 = sp.simplify(J.subs(KI,1))
    print("  x1(K_I=1) =", x1_1)
    print("  J(K_I=1)  =", J_1)

print()
print("="*60)
print("Plateau 铁律参考: HDV 内部解 x1=a → J=a")
print("  若 J(K_I→∞) → a, 说明规模型碎片化后系统退回 plateau (坍缩)")
print("="*60)

# 数值验证
print()
print("数值(a=0.6,α=0.6,s=0.5,γ=1,K_C=1): K_I 扫描")
subs0 = {a:sp.Rational(6,10), alpha:sp.Rational(6,10), s:sp.Rational(5,10),
         gamma:1, KC:1}
for kival in [1,2,5,20,100]:
    Jv = float(J.subs({**subs0, KI:kival}))
    x1v = float(x1s.subs({**subs0, KI:kival}))
    print(f"  K_I={kival:4d}: x1={x1v:.4f}, J={Jv:.4f}")
print(f"  plateau J=a={float(subs0[a]):.4f}")
