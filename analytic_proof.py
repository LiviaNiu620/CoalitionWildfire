"""
解析证明 v2：严格解规模型原子玩家的一阶条件（不预设 regime）。
Braess 对称, 拥堵边 sv,wt 成本=流量, 常数边=1, 捷径=0.

三类流量在对称结构下的变量：
  规模型公司(原子, 内部化组内边际): 走捷径 a3, 走外路 a1(=a2 对称), 2*a1 + a3 = mI
  自私(HDV+竞争型当自私): 走捷径 h3, 走外路 h1(=h2), 2*h1 + h3 = M  (M=1-alpha+mC)
对称拥堵边流量 X = a1 + h1 + a3 + h3  (sv上: 走P1的a1+h1, 走P3的a3+h3)
路径成本: C1 = X + 1 ; C3 = 2X
系统成本 J = 2*X^2 + 2*(a1+h1) + 0

规模型原子内部化: 感知边际 走P1 = C1 + (自己在sv,wt的流量)*c'  ; 但 P1 只用 sv(拥堵)+vt(常数)
  规模型在 sv 的流量 = a1 + a3 ; 在 wt 的流量 = a1 + a3 (对称)
  走 P1 感知边际 = C1 + (a1+a3)*1   [c'=1, P1 只经一条拥堵边 sv]
  走 P3 感知边际 = C3 + (a1+a3)*1 + (a1+a3)*1 = C3 + 2(a1+a3)  [P3 经 sv,wt 两条拥堵边]
自私 Wardrop: C1 vs C3 相等(内部解)。
"""
import sympy as sp

alpha, s = sp.symbols('alpha s', positive=True)
a1, a3, h1, h3, X = sp.symbols('a1 a3 h1 h3 X', nonnegative=True)

mI = alpha*s
M  = 1 - alpha + alpha*(1-s)   # 自私总量 = HDV + 竞争型

# 对称边流量
Xexpr = a1 + h1 + a3 + h3

C1 = X + 1
C3 = 2*X

# 规模型自身拥堵边流量(sv 或 wt) = a1 + a3
avload = a1 + a3
# 规模型感知边际
mc_P1 = C1 + avload*1
mc_P3 = C3 + avload*2

# --- 假设内部解: 规模型 P1,P3 都用(感知边际相等); 自私 P1,P3 都用(实际成本相等) ---
# 但先试: 规模型只走 P1 (a3=0), 自私内部解.
# 情形A: 规模型全外路 a3=0, 规模型感知边际 P1<=P3 (不愿走捷径)
# 自私内部解: C1=C3 => X+1=2X => X=1
# 守恒: 2a1 = mI => a1=mI/2 ; 2h1+h3=M ; X=a1+h1+h3=1 => h1+h3=1-mI/2
# 又需 h3>=0, h1>=0
solA = sp.solve([sp.Eq(2*a1, mI), sp.Eq(2*h1+h3, M),
                 sp.Eq(a1+h1+a3+h3, 1), sp.Eq(a3,0)], [a1,a3,h1,h3], dict=True)
print("情形A(规模型全外路, 自私内部解 X=1):")
if solA:
    sA = solA[0]
    h3A = sp.simplify(sA[h3]); h1A = sp.simplify(sA[h1]); a1A=sp.simplify(sA[a1])
    print("  a1=",a1A," h1=",h1A," h3=",h3A)
    # 检查规模型不愿走捷径: mc_P1 <= mc_P3 at this point
    avl = a1A + 0
    mcP1 = 1+1 + avl   # C1=X+1=2
    mcP3 = 2*1 + 2*avl # C3=2
    print("  规模型 mc_P1=",sp.simplify(mcP1)," mc_P3=",sp.simplify(mcP3),
          " 不走捷径条件 mc_P1<=mc_P3:", sp.simplify(mcP3-mcP1),">=0?")
    # 可行性 h3>=0:
    print("  h3>=0 要求:", sp.solve(h3A>=0, s))
    J_A = 2*(1)**2 + 2*(a1A+h1A)
    print("  J_A(s)=", sp.simplify(J_A))
    print("  dJ_A/ds=", sp.simplify(sp.diff(J_A,s)))
    for sv in [0, sp.Rational(1,2), 1]:
        print(f"    s={float(sv):.1f}: J={float(J_A.subs({alpha:sp.Rational(1,2),s:sv})):.4f}, h3={float(h3A.subs({alpha:sp.Rational(1,2),s:sv})):.4f}")

print()
# 情形B: 当 h3 触零(自私不再用捷径), 规模型内部化"接管", X<1
# 此时自私全走外路 h3=0, 2h1=M => h1=M/2. 规模型内部解 mc_P1=mc_P3.
# X = a1+h1+a3 ; mc_P1=mc_P3: (X+1)+ (a1+a3) = 2X + 2(a1+a3) => X+1+avl = 2X+2avl
#   => 1 = X + avl => X + (a1+a3) = 1
# 守恒规模型: 2a1+a3=mI
print("情形B(自私全外路 h3=0, 规模型内部解):")
h1B = M/2
# X = a1+h1B+a3 ; 条件 X + (a1+a3) = 1
Xdef = a1 + h1B + a3
cond1 = sp.Eq(Xdef + (a1+a3), 1)
cond2 = sp.Eq(2*a1+a3, mI)
solB = sp.solve([cond1, cond2], [a1, a3], dict=True)
if solB:
    sB=solB[0]; a1B=sp.simplify(sB[a1]); a3B=sp.simplify(sB[a3])
    XB = sp.simplify(a1B+h1B+a3B)
    print("  a1=",a1B," a3=",a3B," X=",XB)
    J_B = sp.simplify(2*XB**2 + 2*(a1B+h1B))
    print("  J_B(s)=", J_B)
    print("  dJ_B/ds=", sp.simplify(sp.diff(J_B,s)))
    for sv in [sp.Rational(1,2), sp.Rational(8,10), 1]:
        print(f"    s={float(sv):.1f}: J={float(J_B.subs({alpha:sp.Rational(1,2),s:sv})):.4f}, X={float(XB.subs({alpha:sp.Rational(1,2),s:sv})):.4f}, a3={float(a3B.subs({alpha:sp.Rational(1,2),s:sv})):.4f}")
