"""
缝合 γ：Regime 全图。推 J(s,γ,α) 的分段表达 + regime 边界。
两路载体：拥堵路 c1=x1；备选路 c2=a(常数)。三类：I(内部化),C(spite γ),HDV(自私)。

三类"走拥堵路"量 uI,uC,uH。x1=uI+uC+uH。
感知边际(拥堵路 vs 备选路 a)：
  I: x1 + mI  vs a        (内部化组内边际, c1'=1)
  C: x1 + mC - γ·(uI+uH)  vs a   (spite: 抬高对手(I,HDV)在拥堵路的成本)
  HDV: x1  vs a           (自私)

Regime 由"谁在拥堵路内部解、谁在角点"决定。随 γ 增大：
  R1 Plateau: HDV内部解 x1=a (钉住), I/C 角点
  R2 最优构成: HDV撤出(uH=0), I/C内部解
  R3 竞争过载: C 把 x1 抬过 a, ...
"""
import sympy as sp

a, gamma, alpha, s = sp.symbols('a gamma alpha s', positive=True)
uI, uC, uH = sp.symbols('u_I u_C u_H', nonnegative=True)
mI = alpha*s
mC = alpha*(1-s)

def J_of_x1(x1):
    return sp.simplify(x1**2 + (1-x1)*a)   # x1*c1 + x2*c2, c1=x1,c2=a

print("="*60)
print("Regime 1 (Plateau): HDV 内部解 x1=a")
print("="*60)
# HDV: x1=a. I 想内部化(x1+mI>a → I 不走拥堵路, uI=0角点? 或走). 
# 若 x1=a 被HDV钉住, I感知 = a+mI > a → I 全走备选路 uI=0.
# C感知 = a+mC-γ(uI+uH). uI=0. HDV在拥堵路 uH=a(近似). 
#   C感知 = a+mC-γ·a. 若 < a (即 mC<γa) → C 想走拥堵路(spite划算).
# Plateau 成立条件: C 也不足以打破 → C 温和. J = a (x1=a: J=a²+(1-a)a=a)
J_R1 = J_of_x1(a)
print("  x1=a, J_R1 =", sp.simplify(J_R1))
print("  成立条件(C温和不打破): mC >= γ·a  即 γ <= mC/a = α(1-s)/a")

print()
print("="*60)
print("Regime 2 (最优构成): HDV撤出 uH=0, I/C 内部解")
print("="*60)
# uH=0. I: x1+mI=a. C: x1+mC-γ·uI=a. x1=uI+uC.
eqI = sp.Eq(uI+uC+mI, a)   # x1+mI=a
eqC = sp.Eq(uI+uC+mC-gamma*uI, a)
sol = sp.solve([eqI, eqC],[uI,uC],dict=True)[0]
uI2=sp.simplify(sol[uI]); uC2=sp.simplify(sol[uC]); x1_2=sp.simplify(uI2+uC2)
J_R2=J_of_x1(x1_2)
print("  uI =", uI2, "  uC =", uC2)
print("  x1 =", x1_2)
print("  J_R2 =", sp.simplify(J_R2))
print("  dJ/ds =", sp.simplify(sp.diff(J_R2,s)))
print("  s* =", sp.solve(sp.Eq(sp.diff(J_R2,s),0),s))
print("  R2成立需 uI>=0,uC>=0:")
print("    uI>=0:", sp.simplify(uI2), ">=0")
print("    uC>=0:", sp.simplify(uC2), ">=0")

print()
print("="*60)
print("Regime 3 (竞争过载): C spite 强, 把 x1 抬过 a")
print("="*60)
# γ 很大: C 大量涌入拥堵路抬对手, uI->0(I被逼走), x1>a.
# I全走备选(uI=0). C内部解? 或C也角点全拥堵路 uC=mC.
# x1 = uC + uH; 若HDV也撤 uH=0, x1=uC. C: 若全走拥堵路 uC=mC → x1=mC=α(1-s)
x1_3 = mC  # C全走拥堵路
J_R3 = J_of_x1(x1_3)
print("  x1=mC=α(1-s), J_R3 =", sp.simplify(J_R3))
print("  dJ/ds =", sp.simplify(sp.diff(J_R3.subs(mC,alpha*(1-s)),s)))

print()
print("="*60)
print("Regime 边界(γ 阈值)")
print("="*60)
print("  R1→R2 边界: γ1 = mC/a = α(1-s)/a  (γ超过则C打破plateau)")
print("  R2→R3 边界: uI=0 即", sp.simplify(uI2), "=0 → γ2 =", sp.solve(sp.Eq(uI2,0),gamma))
print()
print("解读: γ 增大 → 依次穿过 R1(plateau)→R2(最优构成)→R3(竞争过载)")
print("  γ 是 regime 切换变量; J 在各 regime 段不同, γ 通过边界进入效率")

# 数值示例: 固定 α=0.5,s=0.3, 扫 γ 看 J (s=0.3 避免 1-2s=0 退化)
print()
print("数值(α=0.5,s=0.3,a=0.6): γ 扫描下的 regime 与 J")
av,sv,aval=0.5,0.3,0.6
mC_v = av*(1-sv); mI_v=av*sv
g1 = mC_v/aval
J_R2_f = sp.lambdify(gamma, J_R2.subs({alpha:av,s:sv}), 'math')
uI2_f = sp.lambdify(gamma, uI2.subs({alpha:av,s:sv}), 'math')
J_R3_v = float(J_R3.subs({alpha:av,s:sv}))
J_R2_v = float(J_R2.subs({alpha:av,s:sv}))  # 不含γ
for gv in [0.2,0.5,1.0,2.0,4.0]:
    if gv <= g1:
        reg="R1 plateau"; J=aval
    else:
        uIv = uI2_f(gv)
        if uIv >= -1e-9:
            reg="R2 最优构成"; J=J_R2_v
        else:
            reg="R3 竞争过载"; J=J_R3_v
    print(f"  γ={gv:.1f}: {reg}, J={J:.4f}  (γ1边界={g1:.3f})")

print()
print("对比: R2 的 s* 与 U 型")
print(f"  s*=a/(2α)={aval/(2*av):.3f}, R2内 J(s)=α²s²-aαs+a 是U型")
print(f"  当前 s={sv}: {'在s*左侧(加规模型改善)' if sv<aval/(2*av) else '在s*右侧(反噬)'}")
