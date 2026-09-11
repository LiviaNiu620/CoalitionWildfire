"""
严格推 Regime Index：单 OD 两平行路载体，三类主体一阶条件联立求闭式。

网络：路 1(拥堵) 成本 c1(x1)=x1 ; 路 2(线性备选) 成本 c2(x2)=a + x2 (a>0常数项).
总需求=1. 三类主体决策"走路1的量"。
  I(规模型 mI, 内部化组内边际): 感知边际(路1) = c1 + mI·c1' = x1 + mI  (c1'=1)
  C(竞争型 mC, spite γ): 感知边际(路1) = c1 + mC·c1' - γ·(对手在路1的边际拥堵)
       = x1 + mC - γ·(∂对手成本/∂自己走路1) ; 对手在路1流量的边际 = c1'·(对手路1量)
  HDV(1-α, 自私): 感知 = c1 = x1
每类在"路1 vs 路2"感知成本相等(内部解)。

变量: I走路1量 uI, C走路1量 uC, HDV走路1量 uH.
x1 = uI+uC+uH ; x2 = 1 - x1.
c1=x1 ; c2 = a + (1-x1).
"""
import sympy as sp

a, gamma, alpha, s = sp.symbols('a gamma alpha s', positive=True)
uI, uC, uH = sp.symbols('u_I u_C u_H', nonnegative=True)
mI, mC = sp.symbols('m_I m_C', positive=True)

x1 = uI + uC + uH
x2 = 1 - x1
c1 = x1
c2 = a + x2
c1p = 1  # c1'

# 各类"走路1"与"走路2"的感知边际相等(内部解)。路2线性 c2'=1,自身边际项对称,简化只在路1加内部化。
# I: 路1感知 = c1 + mI*c1p ; 路2感知 = c2 + mI*1(路2也内部化自身,对称) 
#   内部解: c1 + mI = c2 + mI  → c1=c2? 那内部化抵消了。
#   ⇒ 必须两路拥堵程度不同才有内部化净效应。设路2不拥堵(c2'=0, 常数容量大): c2=a(常数).
# 重设: 路2 = 常数 a (不拥堵备选,类似"逃生路"但这次有限). c2=a, c2'=0.
c2 = a
# I 内部解: c1 + mI*c1p = c2  → x1 + mI = a
# C 内部解(spite): c1 + mC*c1p - γ*(对手路1边际) = c2
#   对手在路1的边际拥堵对 C 走路1的敏感 = c1p*(uI+uH) [I和HDV在路1的量]
#   ⇒ x1 + mC - γ*(uI+uH) = a
# HDV: c1 = c2 → x1 = a   (自私: 路1=路2成本)

# 三个内部解方程(若都内部):
eqI = sp.Eq(x1 + mI, a)
eqC = sp.Eq(x1 + mC - gamma*(uI+uH), a)
eqH = sp.Eq(x1, a)

# 从 eqH: x1=a. 代入 eqI: a+mI=a → mI=0 矛盾(除非 I 在角点).
# ⇒ 三类不能同时内部解. 分 regime. 关键 regime: HDV 内部解钉 x1=a, I/C 在角点.
# 这正是 plateau 铁律! HDV 把 x1 钉在 a. 
# 但 C 的 spite 项能打破: 若 C 走路1"过量"抬高 x1>a, HDV 会撤( uH 减到 0 角点).
# regime 划分:
print("=== Regime 分析(两路) ===")
print("HDV 内部解要求 x1=a. 但若 I+C 的路1总量已 > a, HDV 全撤路1(uH=0,角点).")
print()

# Regime 2: HDV 全撤路1 (uH=0). I,C 内部化决定 x1.
# x1 = uI+uC. I内部解: x1+mI=a? 不,HDV撤了 c2仍=a, I: x1+mI = a → 若成立
# 更清晰: uH=0 时, x1=uI+uC.
# I: x1 + mI = a (I 内部解, 走路1到边际=a)
# C: x1 + mC - γ*uI = a (C spite: 只有 I 在路1可被抬,HDV已撤)
uH0 = 0
x1_r2 = uI + uC
eqI2 = sp.Eq(x1_r2 + mI, a)
eqC2 = sp.Eq(x1_r2 + mC - gamma*uI, a)
sol2 = sp.solve([eqI2, eqC2], [uI, uC], dict=True)
if sol2:
    s2 = sol2[0]
    uI2 = sp.simplify(s2[uI]); uC2 = sp.simplify(s2[uC])
    x1_2 = sp.simplify(uI2 + uC2)
    print("Regime2 (HDV撤路1):")
    print("  uI =", uI2)
    print("  uC =", uC2)
    print("  x1 =", x1_2)
    # 系统成本 J = x1*c1 + x2*c2 = x1^2 + (1-x1)*a
    J2 = sp.simplify(x1_2**2 + (1-x1_2)*a)
    print("  J =", J2)
    # 代入 mI=αs, mC=α(1-s)
    J2_as = sp.simplify(J2.subs({mI: alpha*s, mC: alpha*(1-s)}))
    dJ2 = sp.simplify(sp.diff(J2_as, s))
    print("  J(α,s) =", J2_as)
    print("  dJ/ds =", dJ2)
    crit = sp.solve(sp.Eq(dJ2, 0), s)
    print("  dJ/ds=0 临界 s* =", crit)
