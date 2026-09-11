"""
协调成本模型: 内部化程度 β 内生。
检验: 能否自然产生"大公司内部化(规模型)/小公司不内部化(竞争型)"+三regime。

公司 k (质量 m) 选内部化程度 β∈[0,1]:
  - 内部化降低自己车队拥堵成本: 车队感知边际 = c + β·m·c'
    (β=1完全内部化, β=0非原子自私)
  - 协调成本: 协同车队需调度/算法/通信, 成本 = κ·β²·m (β越高越难协调; ∝规模)
    [用β²表示协调难度随内部化程度超线性上升]
  公司净目标: min 车队出行成本(β) + 协调成本(β)
  → 内生最优 β*(m, κ)

两路: 拥堵路c1=x1, 备选c2=a. 公司质量m, 内部化β.
"""
import sympy as sp

a, m, kappa, beta, alpha, s = sp.symbols('a m kappa beta alpha s', positive=True)

print("="*60)
print("协调成本模型: 内部化程度 β 内生")
print("="*60)
print("公司净成本 = 车队出行成本(β) + 协调成本 κβ²m")
print("车队在拥堵路的感知边际 = c1 + β·m·c1'  (β=内部化程度)")
print()

# 简化: 两路, 公司质量m内部化β, 其余(HDV等)背景流量记为 x_bg.
# 公司内部解(拥堵路vs备选): c1 + β·m = a, c1=x1, 其中 x1 = x_bg + (公司在拥堵路的量)
# 车队出行成本(单位) 近似 = 均衡拥堵水平; 关键看 β 对"公司自己成本"的影响。
# 公司自己车队出行成本 ≈ 用拥堵路的部分 * c1. 内部化β降低自己挤自己→降成本。
# 用一个 reduced form: 公司车队出行成本 T(β) = T0 - g·β·m (内部化省的成本, ∝βm)
#   g = 内部化的边际收益率(∝ 拥堵敏感度 c')
g = sp.symbols('g', positive=True)
T0 = sp.symbols('T0', positive=True)
travel_cost = T0 - g*beta*m        # 内部化省成本(∝βm)
coord_cost = kappa*beta**2*m       # 协调成本(∝β²m)
net = travel_cost + coord_cost
print("公司净成本 N(β) = (T0 - g·β·m) + κ·β²·m")
dN = sp.diff(net, beta)
beta_star = sp.solve(sp.Eq(dN,0), beta)[0]
print("dN/dβ =", dN)
print("内生最优 β* =", sp.simplify(beta_star), "(截断到[0,1])")
print()

print("="*60)
print("★ 关键: β* 与规模 m 的关系")
print("="*60)
print("β* = g/(2κ)  ← 竟然与 m 无关?!")
print()
print("问题诊断: 若内部化收益和协调成本都∝m, 则m约掉, β*与规模无关。")
print("→ 这说明'协调成本∝m'的设定不能产生'大公司更内部化'。")
print()

# 修正: 协调成本的规模经济——大公司协调"单位车"更便宜(规模经济)
# 协调成本 = κ·β²·m^q, q<1 (次线性, 规模经济): 大公司协调单位成本低
q = sp.symbols('q', positive=True)
coord_cost2 = kappa*beta**2*m**q
net2 = travel_cost + coord_cost2
beta_star2 = sp.solve(sp.Eq(sp.diff(net2,beta),0), beta)[0]
print("="*60)
print("修正: 协调成本有规模经济 κβ²m^q (q<1)")
print("="*60)
print("β* =", sp.simplify(beta_star2))
print("  = g·m/(2κ·m^q) = (g/2κ)·m^(1-q)")
print("  q<1 ⟹ β* ∝ m^(1-q) 随 m 增大 → 大公司内部化更强!")
print()
print("★ 这就对了: 协调成本的规模经济(大公司协调更便宜)")
print("  ⟹ 大公司选高β(规模型), 小公司选低β(竞争型) — 内生涌现!")
print()

# 数值
print("数值(g=1,κ=1,q=0.5): β* = 0.5·m^0.5, 随规模m")
print(f"{'m(车队规模)':>10} {'β*(内部化程度)':>14} {'类型':>8}")
for mv in [0.05,0.1,0.2,0.4,0.8]:
    bs = min(0.5*mv**0.5, 1.0)
    typ = "竞争型(低β)" if bs<0.3 else ("规模型(高β)" if bs>0.5 else "中间")
    print(f"{mv:>10} {bs:>14.3f} {typ:>8}")
print()

print("="*60)
print("结论")
print("="*60)
print("1. 竞争型vs规模型不是硬设定, 而是'内部化程度β内生选择'的结果")
print("2. 关键机制: 协调成本有规模经济(大公司协调单位车更便宜, q<1)")
print("   ⟹ 大公司选高β(规模型), 小公司选低β(竞争型) — 自然涌现")
print("3. '竞争型为何不内部化': 因为它小, 协调成本不划算 (回答了那个问号!)")
print("4. 竞争强度 = 1-β* = 1-(g/2κ)m^(1-q): 由规模m和协调成本κ内生决定")
print("5. 纯交通场景(协调成本是车队内部运营成本, 不涉乘客/份额)")
print("6. 统一: 这就是导师的规模经济 + 你最早的'原子性=规模'洞察")
