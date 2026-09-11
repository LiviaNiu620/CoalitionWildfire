"""
推导 + 测试: U型(沿s,中间最优)与 K_I, K_C 的关系。
两路载体, 最新框架(竞争型协同项1-2γ)。

有效协同 E(s) = (αs)²/K_I + (1-2γ)(α(1-s))²/K_C   [atomic内部化∝单家质量²]
J随E反向: E大→x1小(减堵)→J小。 J(s) = g(E(s)), g递减。
简化: 用 J = -E + const 的单调关系看形状(或直接 J=x1²+(1-x1)a, x1=a-E)。

U型(J中间最优, 即E中间最大): 需 E(s) 是凹的(E''<0) 且 顶点s*∈(0,1)。
倒U(J中间最差, E中间最小): E凸(E''>0)。
"""
import sympy as sp

a, alpha, s, gamma, KI, KC = sp.symbols('a alpha s gamma K_I K_C', positive=True)

# 有效协同
E = (alpha*s)**2/KI + (1-2*gamma)*(alpha*(1-s))**2/KC
E = sp.expand(E)
print("有效协同 E(s) =", E)
print()

# E 的一阶、二阶
dE = sp.simplify(sp.diff(E, s))
d2E = sp.simplify(sp.diff(E, s, 2))
print("E'(s) =", dE)
print("E''(s) =", d2E)
print()

# E 的极值点(dE=0)
scrit = sp.solve(sp.Eq(dE, 0), s)
print("E 极值点 s_c =", [sp.simplify(x) for x in scrit])
sc = sp.simplify(scrit[0])
print("  s_c =", sc)
print()

print("="*64)
print("形状判据")
print("="*64)
print("J随E反向。E凹(E''<0)→E中间大→J中间小→U型(J有内部最优)")
print("        E凸(E''>0)→E中间小→J中间大→倒U(J中间最差)")
print()
print("E''=", d2E)
print("E''<0 (U型) ⟺ K_C + K_I(1-2γ) < 0 ⟺ γ > 1/2 + K_C/(2K_I) = γ*")
gstar = sp.Rational(1,2)+KC/(2*KI)
print("γ* =", gstar)
print()

print("="*64)
print("U型的完整条件(两个都要满足)")
print("="*64)
print("条件1(凹): γ > γ* = 1/2 + K_C/(2K_I)")
print("条件2(顶点在内部): s_c ∈ (0,1)")
print("  s_c =", sc)
# s_c 在(0,1)的条件
print("  s_c∈(0,1) 分析:")
# s_c = (1-2γ)KI ... 化简看
sc_s = sp.simplify(sc)
print("  s_c =", sc_s)
print()

# 数值: 扫 K_I, K_C, γ 看 U型是否成立
print("="*64)
print("数值测试: U型(J中间最优)与 K_I,K_C 关系 (α=0.6,a=0.6)")
print("="*64)
import numpy as np
def J_of_s(sv, av, aval, gv, ki, kc):
    Ev = (av*sv)**2/ki + (1-2*gv)*(av*(1-sv))**2/kc
    x1 = aval - Ev
    return x1**2 + (1-x1)*aval

def shape(av, aval, gv, ki, kc):
    ss = np.linspace(0.01,0.99,50)
    Js = np.array([J_of_s(s_,av,aval,gv,ki,kc) for s_ in ss])
    imin=Js.argmin(); imax=Js.argmax()
    # U型: 最小在内部; 倒U: 最大在内部
    if 0<imin<49 and Js[imin]<Js[0]-1e-6 and Js[imin]<Js[-1]-1e-6:
        return f"U型(s*={ss[imin]:.2f})"
    if 0<imax<49 and Js[imax]>Js[0]+1e-6 and Js[imax]>Js[-1]+1e-6:
        return "倒U(中间差)"
    return "单调降" if Js[-1]<Js[0] else "单调升"

print(f"{'K_I':>4}{'K_C':>4} | " + " ".join(f"γ={g}".ljust(14) for g in [0.5,1.0,1.5,2.0]))
for ki in [1,2,4]:
    for kc in [1,2,4]:
        row=f"{ki:>4}{kc:>4} | "
        for gv in [0.5,1.0,1.5,2.0]:
            gc=0.5+kc/(2*ki)
            row += f"{shape(0.6,0.6,gv,ki,kc):<14} "
        print(row)
print()
print("理论 γ* = 1/2 + K_C/(2K_I):")
for ki in [1,2,4]:
    for kc in [1,2,4]:
        print(f"  K_I={ki},K_C={kc}: γ*={0.5+kc/(2*ki):.2f} (γ>γ*才U型)")
