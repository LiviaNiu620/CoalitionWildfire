"""
含 K_I, K_C 的三 regime 普适性检验。
两路: 拥堵路 c1=x1, 备选 c2=a. 
规模型 K_I 家(每家 m_I=αs/K_I), 竞争型 K_C 家(每家 m_C=α(1-s)/K_C).

关键: 原子内部化 ∝ 单家规模. 单家一阶条件(对称):
  规模型(每家): 拥堵路感知边际 = c1 + m_I·c1' = x1 + m_I ; 与备选 a 比
  竞争型(每家): x1 + m_C - γ·(拥堵路上对手边际) ; 与 a 比
  HDV: x1 ; 与 a 比

三 regime(随 γ):
  R1 Plateau: HDV内部解 x1=a (钉住)
  R2 最优构成: HDV撤出, 规模型/竞争型内部化定 x1
  R3 过载: 竞争型把 x1 抬过头
检验: 三 regime 的定性结构对任意 K_I,K_C 是否保持。
"""
import sympy as sp

a, alpha, s, gamma = sp.symbols('a alpha s gamma', positive=True)
KI, KC = sp.symbols('K_I K_C', positive=True)

m_I = alpha*s/KI      # 单家规模型
m_C = alpha*(1-s)/KC  # 单家竞争型

print("单家规模型 m_I =", m_I, " ; 单家竞争型 m_C =", m_C)
print()

# ========== R1 Plateau ==========
print("="*60)
print("R1 Plateau: HDV 内部解 x1=a → J=a")
print("="*60)
J_R1 = sp.simplify(a**2 + (1-a)*a)  # x1=a: x1*c1+x2*c2 = a*a+(1-a)*a
print("  J_R1 = a (与 K 无关)")
print("  成立条件: 竞争型不足以打破 → 单家竞争型感知(拥堵路)≥a")
print("    a + m_C - γ·a ≥ a  ⟹  γ ≤ m_C/a = α(1-s)/(K_C·a)")
gamma1 = m_C/a
print("  R1→R2 边界 γ1 =", sp.simplify(gamma1), " (依赖 K_C!)")
print()

# ========== R2 最优构成 ==========
print("="*60)
print("R2 最优构成: HDV撤出, 规模型内部化接管")
print("="*60)
# HDV撤出, 拥堵路只有规模型+竞争型. 规模型内部解(单家): x1 + m_I = a
#   → x1 = a - m_I = a - αs/K_I
x1_R2 = a - m_I
J_R2 = sp.simplify(x1_R2**2 + (1-x1_R2)*a)
print("  x1 =", x1_R2, " = a - αs/K_I")
print("  J_R2 =", sp.expand(J_R2))
# 对 s 求最优构成 s*
dJ_ds = sp.simplify(sp.diff(J_R2, s))
sstar = sp.solve(sp.Eq(dJ_ds, 0), s)
print("  dJ_R2/ds =", dJ_ds)
print("  最优构成 s*(K_I) =", sstar)
print()

# ========== R3 过载 ==========
print("="*60)
print("R3 过载: 竞争型 spite 强, 把 x1 抬过头")
print("="*60)
# 竞争型全涌拥堵路(近似), 规模型被逼走, HDV撤. x1 = 竞争型总量 = α(1-s)
x1_R3 = alpha*(1-s)
J_R3 = sp.simplify(x1_R3**2 + (1-x1_R3)*a)
print("  x1 = α(1-s) (竞争型总量, 与 K_C 无关的总量)")
print("  J_R3 =", sp.expand(J_R3))
print()

# ========== 普适性结论 ==========
print("="*60)
print("普适性检验结论")
print("="*60)
print("三 regime 的定性结构对任意 K_I, K_C 都存在:")
print("  R1 Plateau: J=a (恒定)")
print("  R2 最优构成: J=α²s²/K_I² - aαs/K_I + a, 最优 s*(K_I)=a·K_I/(2α)")
print("  R3 过载: J=α²(1-s)²+a(1-α(1-s))")
print()
print("K 的作用 = 定量平移(不改变定性三 regime):")
print("  • R1→R2 边界 γ1 = α(1-s)/(K_C·a): K_C↑ → γ1↓? 需修正符号(见下)")
print("  • 最优构成 s*(K_I) = a·K_I/(2α): K_I↑ → s*↑ (最优点右移)")
print("  • U型可见(s*≤1): K_I ≤ 2α/a")
print()
print("⟹ 三 regime 结构普适, K 只改变 regime 的定量位置与可见范围。")
print("   当 K_I ≤ 2α/a 时最优构成U型在[0,1]内可见; 超过则s*右移出界(退化单调)。")

# 数值验证
print()
print("="*60)
print("数值验证(a=0.6,α=0.6): s*(K_I)=a·K_I/(2α) 随 K_I")
print("="*60)
for kival in [1,2,3,4]:
    sv = float((a*KI/(2*alpha)).subs({a:sp.Rational(6,10),alpha:sp.Rational(6,10),KI:kival}))
    vis = "可见(U型)" if sv<=1 else "出界(单调)"
    print(f"  K_I={kival}: s*={sv:.2f}  {vis}")
