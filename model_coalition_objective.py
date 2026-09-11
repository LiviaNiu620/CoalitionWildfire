"""
根本修正后的模型: 两类都是协同coalition, 区别在协同目标(减堵vs抢占)。
- HDV: 非原子, 无协同项(散点)。
- 规模型: 协同减堵, 协同项 +Σf_e c'_e。
- 竞争型: 协同但目标偏抢占, 协同项 (1-2γ)Σf_e c'_e。
  γ=0: +1(纯减堵,=规模型); γ=0.5: 0(中性协同); γ=1: -1(纯抢占)。
  竞争强度 γ = 协同优势用于抢占而非减堵的程度。
关键: 竞争型始终有协同项(是coalition,享协同优势,≠HDV), 只是方向随γ变。

两路载体: 拥堵路c1=x1, 备选c2=a.
"""
import sympy as sp

a, alpha, s, gamma, KI = sp.symbols('a alpha s gamma K_I', positive=True)

print("="*60)
print("修正框架: 协同项方向 = 协同目标")
print("="*60)
print("HDV(非原子): 感知 = c1 (无协同项)")
print("规模型(协同减堵): 感知 = c1 + m·c1'")
print("竞争型(协同,目标偏抢占): 感知 = c1 + (1-2γ)·m·c1'")
print("  γ=0→+m·c1'(减堵=规模型); γ=0.5→0(中性协同); γ=1→-m·c1'(抢占)")
print("  ★竞争型始终有协同项(是coalition,≠HDV), 方向随γ变")
print()

# 规模型质量 αs, 竞争型质量 α(1-s), HDV 1-α
# 竞争型协同系数 (1-2γ). 为得J(K_I,s,γ), 分析regime.
# R2: HDV撤出, coalition们内部化定x1.
# 规模型内部解(每家m_I=αs/K_I): c1 + m_I·1 = a → 贡献
# 竞争型内部解(每家m_C=α(1-s)/K_C, 协同系数(1-2γ)): c1 + (1-2γ)m_C = a
# 简化: 单规模型coalition(K_I=1,质量αs) + 单竞争型coalition(质量α(1-s))
# 两个coalition + 自私HDV. 
# 对称假设下, 拥堵路总流量 x1 由 coalition 内部解决定。
# 规模型: x1 + αs = a  (若规模型内部解主导)
# 竞争型: x1 + (1-2γ)α(1-s) = a

# 场景分析: 两类coalition都在拥堵路上, HDV撤出。
# 规模型想 x1 小(减堵), 竞争型 γ大时想 x1 大(抢占)。
# 均衡 x1 由两类协同项加权。设都内部解:
# 规模型: c1 + αs = a ; 竞争型: c1 + (1-2γ)α(1-s) = a
# 这两个一般不同时成立(除非系数配平). 取"有效协同" = 规模型减堵 - 竞争型抢占净效应。
# 有效内部化 = αs·1 + α(1-s)·(1-2γ)  [规模型全内部化 + 竞争型(1-2γ)]
eff_intern = alpha*s*1 + alpha*(1-s)*(1-2*gamma)
eff_intern = sp.simplify(eff_intern)
print("有效协同(内部化)总量 = αs·(+1) + α(1-s)·(1-2γ) =", eff_intern)
print("  = α[s + (1-s)(1-2γ)] = α[1 - 2γ(1-s)]")
eff = sp.simplify(alpha*(1 - 2*gamma*(1-s)))
print("  =", eff)
print()

# x1 = a - 有效协同 (拥堵路流量被有效内部化压低)
x1 = a - eff
J = sp.simplify(x1**2 + (1-x1)*a)
print("x1 = a - α[1-2γ(1-s)] =", x1)
print("J(s,γ) =", sp.expand(J))
print()

print("="*60)
print("三个方向的效应")
print("="*60)
# 沿 s (规模型占比), 固定γ
dJs = sp.simplify(sp.diff(J, s))
print("dJ/ds =", dJs)
# 沿 γ (竞争强度), 固定s
dJg = sp.simplify(sp.diff(J, gamma))
print("dJ/dγ =", dJg)
print()
print("解读:")
print("• γ=0: 竞争型也减堵(=规模型), 有效协同=α, x1=a-α(最低堵)")
print("• γ=0.5: 竞争型中性, 有效协同=α·s(只规模型减堵)")
print("• γ=1: 竞争型抢占, 有效协同=α(2s-1), s<0.5时为负(净抢占→加堵)")
print()

# 数值
print("数值(a=0.6,α=0.6): 有效协同 α[1-2γ(1-s)] 与 J")
print(f"{'γ\\s':>6}", end="")
for sv in [0.2,0.5,0.8]:
    print(f"  s={sv:<8}", end="")
print()
for gv in [0.0,0.5,1.0]:
    print(f"γ={gv:<4}", end="")
    for sv in [0.2,0.5,0.8]:
        jv=float(J.subs({a:0.6,alpha:0.6,s:sv,gamma:gv}))
        print(f"  J={jv:<8.3f}", end="")
    print()

print()
print("="*60)
print("结论(修正框架)")
print("="*60)
print("有效协同 = α[1 - 2γ(1-s)]:")
print("  s↑(规模型多) → 有效协同↑ → J↓ (更多减堵型coalition)")
print("  γ↑(竞争型偏抢占) → 有效协同↓ → J↑ (协同用于抢占)")
print("  γ大且s小: 有效协同为负 → 竞争型coalition协同抢占加剧拥堵")
print("★ 两类都是coalition(都有协同项), 区别在协同目标(方向±)")
print("★ 竞争强度γ = 协同用于抢占而非减堵的程度, 竞争型≠HDV(HDV无协同项)")
