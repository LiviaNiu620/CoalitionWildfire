"""
彻底查清: U型(沿s,中间最优)到底依赖哪些变量?
严格推导 + 全变量数值扫描 (α, γ, K_I, K_C, a)。
"""
import sympy as sp
import numpy as np

# ============ 严格推导 ============
a, alpha, s, gamma, KI, KC = sp.symbols('a alpha s gamma K_I K_C', positive=True)
E = (alpha*s)**2/KI + (1-2*gamma)*(alpha*(1-s))**2/KC
# J = x1²+(1-x1)a, x1=a-E
x1 = a - E
J = sp.expand(x1**2 + (1-x1)*a)
dJ = sp.simplify(sp.diff(J, s))
d2J = sp.simplify(sp.diff(J, s, 2))
print("="*66)
print("严格推导: J(s) 的二阶导 (决定凹凸/U型)")
print("="*66)
print("J''(s) =", d2J)
# 内部最优点(dJ=0)
scrit = sp.solve(sp.Eq(dJ,0), s)
print("内部最优 s* =", [sp.simplify(x) for x in scrit])
sc = sp.simplify(scrit[0])
print("  s* =", sc)
print()
print("★ s* 表达式里出现哪些变量:", sc.free_symbols)
print("★ J''  表达式里出现哪些变量:", d2J.free_symbols)
print()

# ============ 全变量数值扫描 ============
def J_of_s(sv, av, aval, gv, ki, kc):
    Ev=(av*sv)**2/ki + (1-2*gv)*(av*(1-sv))**2/kc
    x1=aval-Ev
    return x1**2+(1-x1)*aval

def is_Ushape(av, aval, gv, ki, kc):
    ss=np.linspace(0.01,0.99,80)
    Js=np.array([J_of_s(x,av,aval,gv,ki,kc) for x in ss])
    imin=Js.argmin()
    return 0<imin<79 and Js[imin]<Js[0]-1e-7 and Js[imin]<Js[-1]-1e-7, (ss[imin] if 0<imin<79 else -1)

print("="*66)
print("数值扫描: 逐个变量看U型是否出现")
print("="*66)

print("\n[1] 固定 K_I=1,K_C=1,a=0.6, 扫 α 和 γ:")
print(f"{'α\\γ':>6}", end="")
for gv in [0.3,0.6,1.0,1.5]: print(f" γ={gv:<5}", end="")
print()
for av in [0.2,0.4,0.6,0.8,1.0]:
    print(f"{av:>6}", end="")
    for gv in [0.3,0.6,1.0,1.5]:
        u,ss=is_Ushape(av,0.6,gv,1,1)
        print(f" {'U'+f'{ss:.2f}' if u else '单调':<7}", end="")
    print()

print("\n[2] 固定 α=0.6,a=0.6,γ=1, 扫 K_I 和 K_C:")
print(f"{'K_I\\K_C':>8}", end="")
for kc in [1,2,4,8]: print(f" K_C={kc:<4}", end="")
print()
for ki in [1,2,3,4,8]:
    print(f"{ki:>8}", end="")
    for kc in [1,2,4,8]:
        u,ss=is_Ushape(0.6,0.6,1.0,ki,kc)
        print(f" {'U'+f'{ss:.2f}' if u else '单调':<6}", end="")
    print()

print("\n[3] 固定 K_C=1,a=0.6,γ=1, 扫 K_I 和 α:")
print(f"{'K_I\\α':>7}", end="")
for av in [0.3,0.5,0.7,0.9]: print(f" α={av:<5}", end="")
print()
for ki in [1,2,3,4]:
    print(f"{ki:>7}", end="")
    for av in [0.3,0.5,0.7,0.9]:
        u,ss=is_Ushape(av,0.6,1.0,ki,1)
        print(f" {'U'+f'{ss:.2f}' if u else '单调':<6}", end="")
    print()

print("\n" + "="*66)
print("结论")
print("="*66)
print("看 s* 和 J'' 的变量依赖 + 数值扫描, 判断U型依赖哪些变量。")
