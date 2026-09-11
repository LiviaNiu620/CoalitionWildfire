"""
阶段2.2+2.3: 一般网络推广 s*、Regime Index、含K普适。
两路载体的 s*=a/(2α) 依赖具体参数(拥堵边斜率=1, 备选=a)。
推广到一般拥堵边 c1(x)=t0+t1·x, 备选 c2=a, 看 s* 的一般形式。
"""
import sympy as sp

alpha, s, a, t0, t1, KI = sp.symbols('alpha s a t_0 t_1 K_I', positive=True)

print("="*60)
print("2.2 推广到一般拥堵边 c1(x)=t0+t1·x, 备选 c2=a")
print("="*60)
# R2: HDV撤出, 规模型内部化接管. 规模型内部解(拥堵路vs备选):
#   c1(x1) + β·m·c1'(x1) = a, β=1(规模型), m=αs, c1'=t1
#   t0 + t1·x1 + αs·t1 = a  → x1 = (a - t0 - αs·t1)/t1 = (a-t0)/t1 - αs
x1 = (a - t0)/t1 - alpha*s
print("规模型内部解 x1 =", sp.simplify(x1))
# 系统成本 J = x1·c1(x1) + (1-x1)·a
c1 = t0 + t1*x1
J = sp.expand(x1*c1 + (1-x1)*a)
J = sp.simplify(J)
print("J(s) =", J)
dJ = sp.simplify(sp.diff(J, s))
print("dJ/ds =", dJ)
sstar = sp.solve(sp.Eq(dJ,0), s)
print("s* =", [sp.simplify(x) for x in sstar])
print()

# 简化: t0=0,t1=1 应回到 a/(2α)
sstar_simple = sp.simplify(sstar[0].subs({t0:0,t1:1}))
print("验证 t0=0,t1=1: s* =", sstar_simple, " (应=a/(2α))")
print()

print("="*60)
print("一般 s* 形式 + Regime Index")
print("="*60)
sstar_gen = sp.simplify(sstar[0])
print("一般 s* =", sstar_gen)
print()
# Regime Index: dJ/ds 符号
print("dJ/ds =", dJ, " → 归一化 Regime Index:")
# dJ/ds = 0 at s*; 写成 R vs 1
print("R = s/s* (R<1改善, R>1反噬), s* =", sstar_gen)
print()

print("="*60)
print("2.3 含 K 普适 (规模型分 K_I 家, 每家 m=αs/K_I)")
print("="*60)
# 规模型内部解(单家): c1(x1) + (αs/K_I)·t1 = a
#   t0 + t1·x1 + (αs/K_I)·t1 = a → x1 = (a-t0)/t1 - αs/K_I
x1_K = (a-t0)/t1 - alpha*s/KI
c1_K = t0 + t1*x1_K
J_K = sp.simplify(sp.expand(x1_K*c1_K + (1-x1_K)*a))
print("x1(K_I) =", sp.simplify(x1_K))
print("J(s,K_I) =", J_K)
dJ_K = sp.simplify(sp.diff(J_K, s))
sstar_K = sp.solve(sp.Eq(dJ_K,0), s)
print("s*(K_I) =", [sp.simplify(x) for x in sstar_K])
print()
sstar_K_simple = sp.simplify(sstar_K[0].subs({t0:0,t1:1}))
print("验证 t0=0,t1=1: s*(K_I) =", sstar_K_simple, " (应=a·K_I/(2α))")
print()

print("="*60)
print("推广结论")
print("="*60)
print("• 一般拥堵边 c1=t0+t1x, 备选a:")
print("  s* = (a-t0)/(2αt1)  [t0=0,t1=1时回到 a/(2α)]")
print("• 含K普适: s*(K_I) = (a-t0)K_I/(2αt1)  [回到 aK_I/(2α)]")
print("• Regime Index 归一化: R = s/s*, 结构不变(<1改善,>1反噬)")
print("• U型可见: s*≤1 ⟺ K_I ≤ 2αt1/(a-t0)")
print()
print("⟹ s*、Regime Index、含K普适的结构对一般(线性)拥堵边成立,")
print("   两路 a/(2α) 是 t0=0,t1=1 的特例。凸非线性边需数值,结构预期保持。")
