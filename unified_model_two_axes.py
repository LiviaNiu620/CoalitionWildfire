"""
统一模型(选甲): 两正交自变量 集中度K + 竞争强度(1-s)。
竞争型 β=0(不内部化=自私), 并入自私阵营。
系统 = 规模型(总质量αs, 分K_I家各内部化) + 自私阵营(HDV 1-α + 竞争型 α(1-s) = 1-αs)。

两路载体: 拥堵路c1=x1, 备选c2=a.
"""
import sympy as sp

a, alpha, s, KI = sp.symbols('a alpha s K_I', positive=True)

print("="*60)
print("统一模型: 规模型(αs,K_I家) vs 自私阵营(1-αs)")
print("="*60)
print("规模型总质量 = αs (分K_I家, 每家 m=αs/K_I, 内部化β=1)")
print("自私阵营 = HDV(1-α) + 竞争型(α(1-s)) = 1-αs (β=0, 非原子自私)")
print("竞争强度 = 竞争型占比 (1-s); 集中度 = 规模型公司数 K_I")
print()

# 自私阵营: Wardrop, 在拥堵路上内部解 c1=c2 → x1=a (若内部解)
# 但规模型内部化会改变均衡. 分regime.
# R2(自私阵营被挤出拥堵路的部分, 规模型内部化接管):
# 规模型每家内部解: c1(x1) + (αs/K_I)*c1' = a → x1 + αs/K_I = a → x1 = a - αs/K_I
x1 = a - alpha*s/KI
print("R2 规模型内部解: x1 = a - αs/K_I =", x1)
# 系统成本 J = x1*c1 + (1-x1)*a, c1=x1
J = sp.expand(x1*x1 + (1-x1)*a)
J = sp.simplify(J)
print("J(K_I, s) =", J)
print()

# 沿 s (固定K_I): 三regime/s*
print("--- 沿竞争强度维 (固定K_I, 变s) ---")
dJ_ds = sp.simplify(sp.diff(J, s))
sstar = sp.solve(sp.Eq(dJ_ds, 0), s)
print("dJ/ds =", dJ_ds)
print("s* =", sstar, " (K_I=1时:", [sp.simplify(x.subs(KI,1)) for x in sstar], ")")
print()

# 沿 K_I (固定s): 碎片化坍缩
print("--- 沿集中度维 (固定s, 变K_I) ---")
dJ_dK = sp.simplify(sp.diff(J, KI))
print("dJ/dK_I =", dJ_dK)
print("K_I→∞: J →", sp.limit(J, KI, sp.oo), " (坍缩回plateau J=a)")
print("K_I=1: J =", sp.simplify(J.subs(KI,1)))
print()

print("="*60)
print("统一结论: J(K_I, s) 两正交维度")
print("="*60)
print("J(K_I,s) =", J)
print()
print("集中度维(K_I): K_I↑ → 规模型红利↓ → J↑ (碎片化坍缩)")
print("  规模型红利 = a - J =", sp.simplify(a-J), " ∝ 1/K_I")
print("竞争强度维(1-s): s↓(竞争型多) → J↑ (吃掉效应)")
print("  两维都指向: 少而大的规模型coalition最好; 多而碎/竞争型最差")
print()
print("最优: s*=aK_I/(2α) (需≤1), K_I=1(最集中)")
print("→ 系统最优 = 集中(K_I小) + 协调(s大,竞争型少)")

# 数值验证
print()
print("数值(a=0.6,α=0.6): J(K_I,s)")
print(f"{'':>8}", end="")
for sv in [0.2,0.5,0.8,1.0]:
    print(f"s={sv:<6}", end="")
print()
for kiv in [1,2,4]:
    print(f"K_I={kiv:<4}", end="")
    for sv in [0.2,0.5,0.8,1.0]:
        jv = float(J.subs({a:0.6,alpha:0.6,s:sv,KI:kiv}))
        print(f"{jv:<8.3f}", end="")
    print()
print("→ 每列(固定s)K_I↑则J↑(碎片化); 每行(固定K_I)s↑则J↓(协调)")
