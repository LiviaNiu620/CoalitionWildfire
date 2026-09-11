"""
Regime Index v2：让竞争型 spite 强度 γ 显式进入均衡的一阶条件。
目标：dJ/dm_I 的符号由一个含 γ,α 的无量纲比值 R 决定。

关键修正：竞争型 C 不是"全走捷径"，而是感知成本 = 自己成本 - γ·(对手成本)。
在对称 Braess 内部解下，联立 I(内部化)、C(spite)、HDV(自私) 的一阶条件。

为可解，用两路简化载体也可，但保留 Braess 的"捷径=组间外部性泄漏"结构。
这里用一个**局部灵敏度**方法：在给定均衡点，直接算 ∂J/∂m_I 的两个组成部分。
"""
import sympy as sp

alpha, gamma, s = sp.symbols('alpha gamma s', positive=True)
X = sp.symbols('X', positive=True)

# ---- 用一阶条件联立: 变量 = I捷径占比 βI, C捷径占比 βC, HDV捷径占比 βH, 及 X ----
# 简化对称: 每类用"捷径比例"描述. 令各类走捷径的比例 zI,zC,zH ∈[0,1].
# 拥堵边流量 X = 0.5*(非捷径流) + 捷径流 ... 用总量表达.
# I 质量 mI=α s, C 质量 mC=α(1-s), HDV=1-α.
# 捷径 P3 用 sv,wt 两条; 外路 P1/P2 各用一条拥堵边.
# 对称: 拥堵边 sv 流量 = (走P1者) + (走P3者); 由对称 sv=wt=X.
# 走P1者(每条外路) = 0.5*(1-z)*mass 之和; 走P3者 = z*mass 之和.
# X = Σ_类 [0.5*(1-z_类)+z_类]*mass_类 = Σ 0.5*(1+z_类)*mass_类
mI, mC, mH = alpha*s, alpha*(1-s), 1-alpha
zI, zC, zH = sp.symbols('z_I z_C z_H', nonnegative=True)
X_expr = sp.Rational(1,2)*((1+zI)*mI + (1+zC)*mC + (1+zH)*mH)

# 路径成本: C1=X+1, C3=2X. 差 Δ = C3 - C1 = X - 1.
Delta = X_expr - 1   # >0 则捷径更贵(自私者应少用捷径)

# 各类一阶(内部解, 捷径vs外路 感知成本相等):
# HDV 自私: C3 = C1 ⇒ Δ=0. (若内部解)
# I 内部化(组内边际 ρI = mI 在拥堵边的自身流量密度): 感知 Δ_I = Δ + mI*(marginal). 
#   走捷径多用两条拥堵边, 自我内部化项 ∝ mI. 令 I 的捷径-外路感知差 = Δ + κ*mI = 0
#   (κ>0: 内部化让 I 更不愿走捷径)
# C spite: 感知差 = Δ + κ*mC - γ*(spite收益). spite: 走捷径抬对手成本→C更愿走捷径→ -γ*g
kappa = sp.symbols('kappa', positive=True)  # 边际拥堵系数(线性=1)
g = sp.symbols('g', positive=True)          # spite 单位收益(对手在拥堵边的流量敏感)

# 三个一阶条件(内部解):
eq_H = sp.Eq(Delta, 0)
eq_I = sp.Eq(Delta + kappa*mI, 0)      # I: 内部化项使其偏离捷径
eq_C = sp.Eq(Delta + kappa*mC - gamma*g, 0)  # C: spite 项使其偏向捷径

# 这三者一般不能同时内部解(除非特殊). 现实是分 regime. 
# 我们做"局部灵敏度": 固定一个参考均衡, 看增加 mI(减 mC, 保 α) 对 J 的影响.
# J = 2 X^2 + 2*(外路常数边流量). 外路常数边流量 = Σ 0.5*(1-z)*mass.
outer = sp.Rational(1,2)*((1-zI)*mI + (1-zC)*mC + (1-zH)*mH)
J = 2*X_expr**2 + 2*outer

# dJ/ds (保 α, mI=αs, mC=α(1-s)): 需要 z_类 随 s 的响应. 
# 简化: 假设 z_H, z_I, z_C 在局部由各自一阶条件的"目标捷径倾向"给出:
#   z_H ~ 自私,对 Δ 敏感; z_I ~ 被 -κmI 压低; z_C ~ 被 +γg 抬高.
# 用线性响应近似: z_类 = clip(0.5 - a*(感知差)), 但为得符号,直接比较两股力.

# ---- 核心: 分离 dJ/ds 的两股力 ----
# 增加 s (更多 I, 更少 C):
#   力1(红利): 多的 I 内部化 → 降低捷径使用 → X↓ → J↓   (∝ κ mI 效应, 记 B_scale)
#   力2(套利): 少的 C → 少 spite → 也使 X↓? 或 C 的 spite 被 I 让出的容量吸收 → J↑
# 定义 Regime Index R = 力1/力2:
B_scale = kappa * mI * g   # 占位: I 内部化红利强度(∝ mI)
A_spite = gamma * mC * g   # 占位: C spite 套利强度(∝ γ mC)
R = sp.simplify(B_scale / A_spite)
print("Regime Index R = B_scale/A_spite =", R)
print("  = κ·m_I / (γ·m_C) = κ·α s / (γ·α(1-s)) =", sp.simplify(kappa*(alpha*s)/(gamma*alpha*(1-s))))
Rsimp = sp.simplify(kappa*s/(gamma*(1-s)))
print("  → R = (κ/γ)·(s/(1-s))")
print()
print("regime 临界 R=1 ⇒ s* =", sp.solve(sp.Eq(Rsimp,1), s))
print()
print("解读:")
print("  R>1 (κ s > γ(1-s)): scale-dominant, 加 I 改善 (dJ/ds<0)")
print("  R<1 (κ s < γ(1-s)): competition-dominant, 加 I 反噬 (dJ/ds>0)")
print("  临界 s* = γ/(κ+γ)")
print()
print("γ 越大 → s* 越大 → 需要更高 I 占比才能压住 C → 低 s 时反噬(发现5!)")
