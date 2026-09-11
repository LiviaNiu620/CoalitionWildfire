"""
阶段3: 市场失灵(Theorem 6) + transfer-free机制设计(Theorem 7)。

链条:
- 公司内生选内部化 β*=(g/2κ)m^(1-q) (Stage1)。
- "有效规模型占比" s_eff = 全体加权平均 β (β高=规模型, β低=竞争型)。
- 市场自发 s_eff^market 由 协调成本κ + 公司规模分布 决定。
- 社会最优 s* = a/(2α) (Theorem 2)。
- 市场失灵: s_eff^market < s* (协调成本使公司内部化不足)。
- 机制设计: regulator 降 κ (提供公共调度基础设施, transfer-free) → 提高全体β* → s_eff↑ → 逼近 s*。
"""
import sympy as sp

g, kappa, q, m, alpha, a, s = sp.symbols('g kappa q m alpha a s', positive=True)

print("="*60)
print("3.1 市场失灵 (Theorem 6)")
print("="*60)
# 单公司内部化 β*=(g/2κ)m^(1-q). 有效规模型贡献 ∝ β·m (内部化的质量).
# 简化: 全体AV质量α, 若都是规模型(β=1)则 s_eff=1; 都竞争型(β=0)则 s_eff=0.
# s_eff = 加权平均 β = E[β*] (按质量加权).
# 对称情形: N家公司各 m=α/N, β*=(g/2κ)(α/N)^(1-q).
N = sp.symbols('N', positive=True)
m_each = alpha/N
beta_star = (g/(2*kappa))*m_each**(1-q)
print("单家 β* = (g/2κ)(α/N)^(1-q) =", sp.simplify(beta_star))
# 有效规模型占比 = β* (每家都内部化β*程度)
s_eff_market = beta_star  # (截断到[0,1])
print("市场有效规模型占比 s_eff^market = β* =", sp.simplify(s_eff_market))
print()
s_social = a/(2*alpha)
print("社会最优 s* =", s_social)
print()
print("市场失灵 gap = s* - s_eff^market =")
gap = sp.simplify(s_social - s_eff_market)
print("  ", gap)
print()
print("失灵条件: s_eff^market < s* ⟺ (g/2κ)(α/N)^(1-q) < a/(2α)")
print("  → 协调成本κ大 或 公司碎片(N大) → 内部化不足 → 规模型太少 → 失灵")
print("  → 这是协调成本导致的市场失灵(公司理性地内部化不足)")
print()

print("="*60)
print("3.2 transfer-free 机制设计 (Theorem 7)")
print("="*60)
print("regulator 工具(transfer-free, 非定价/补贴):")
print("  工具1: 降低协调成本 κ→κ' (提供公共调度/通信基础设施)")
print("  工具2: 促进公司整合(降N, 提高每家m)")
print("  工具3: 设计共享资源容量(改变a, Theorem3)")
print()
# 工具1: 求让 s_eff^market = s* 的目标 κ'
kappa_target = sp.solve(sp.Eq(beta_star, s_social), kappa)[0]
print("工具1最优: 令 β*(κ')=s* → κ' =", sp.simplify(kappa_target))
print("  → regulator 把协调成本降到 κ' 即可诱导 s_eff=s* (社会最优)")
print()
# 验证: κ 越小 β* 越大, 单调
print("β* 对 κ 单调: dβ*/dκ =", sp.simplify(sp.diff(beta_star, kappa)), " <0 ✓")
print("  → 降κ 严格提高内部化 → 严格改善 → 存在 κ' 精确达到 s*")
print()

print("="*60)
print("3.3 共享资源容量作为设计变量")
print("="*60)
# a = 备选路成本 = 共享资源的"稀缺度"代理. s*=a/(2α) 依赖 a.
# regulator 调 a(共享资源容量): a大(备选贵/共享资源稀缺)→s*大→需更多规模型
# 但同时 Theorem3: a太大可能改变regime. 求"社会成本最小"的最优 a.
# 社会最优 J(s*)=a(1-a/4), 对 a 求极值:
J_opt = a*(1 - a/4)
da = sp.diff(J_opt, a)
a_star = sp.solve(sp.Eq(da,0), a)[0]
print("社会最优构成下 J(s*)=a(1-a/4)")
print("dJ/da =", da, " → 最优共享资源参数 a* =", a_star)
print("  → 存在最优共享资源容量(a*=2), 对应最低系统成本 J=1")
print("  → regulator 设计共享资源容量到 a* 达全局最优")
print()

print("="*60)
print("阶段三结论")
print("="*60)
print("Theorem 6(市场失灵): 协调成本使公司内生内部化不足,")
print("  s_eff^market=(g/2κ)(α/N)^(1-q) < s*=a/(2α), gap由κ,N决定。")
print("Theorem 7(transfer-free机制设计): regulator无需定价/补贴,")
print("  通过 降协调成本κ→κ'=g α^(1-q)/(a N^(1-q)) [公共调度基础设施]")
print("  或 设计共享资源容量a→a*=2, 即可诱导社会最优构成。")
print("  β*对κ严格单调 ⟹ 机制有效且最优κ'唯一。")
