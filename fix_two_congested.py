"""
修复: 双拥堵路模型。两条路都用拥堵函数。
路A(共享瓶颈): c_A(x_A)=A0+A1·x_A
路B(备选):     c_B(x_B)=B0+B1·x_B   (B1>0, 不再常数)
需求: 总量N, OD从O到D, 两路并联. AV总量αN(K个对称coalition), HDV=(1-α)N.
决策: coalition每家上A量 a_k; HDV上A量 y. x_A=Σa_k+y, x_B=N-x_A.
"""
import sympy as sp

N,K,gamma,alpha=sp.symbols('N K gamma alpha',positive=True)
A0,A1,B0,B1=sp.symbols('A0 A1 B0 B1',positive=True)
ak,y=sp.symbols('a_k y',nonnegative=True)

print("="*74)
print("双拥堵路模型 — 均衡条件")
print("="*74)
# 对称: K家各上A量ak, x_A=K*ak+y, x_B=N-x_A
xA=K*ak+y
xB=N-xA
cA=A0+A1*xA
cB=B0+B1*xB
# coalition目标 J_k=ak*cA+(g-ak)*cB - γ ak² A1  (g=αN/K)
g=alpha*N/K
# coalition FOC: ∂/∂ak [ak*cA+(g-ak)*cB-γak²A1]
Jk=ak*cA+(g-ak)*cB-gamma*ak**2*A1
foc_co=sp.expand(sp.diff(Jk,ak))
print("coalition FOC (∂J_k/∂a_k=0):")
print("  ",foc_co,"= 0")
# HDV FOC: cA=cB (Wardrop)
foc_hdv=sp.Eq(cA,cB)
print("HDV FOC: c_A = c_B  ⟹", sp.expand(cA-cB),"= 0")
print()
print("="*74)
print("联立求解 (对称内部解 a_k=a*, y*)")
print("="*74)
a,yy=sp.symbols('a y',nonnegative=True)
xA_=K*a+yy; xB_=N-xA_
cA_=A0+A1*xA_; cB_=B0+B1*xB_
# coalition FOC (用a,yy)
Jk_=a*cA_+(g-a)*cB_-gamma*a**2*A1
eq1=sp.Eq(sp.diff(Jk_,a),0)
# 注意: ∂cB/∂a: x_B=N-Ka-y, ∂/∂a=-K. 已被diff自动处理
eq2=sp.Eq(cA_-cB_,0)  # HDV Wardrop
sol=sp.solve([eq1,eq2],[a,yy],dict=True)
print("解:",)
for s_ in sol:
    a_star=sp.simplify(s_[a]); y_star=sp.simplify(s_[yy])
    print("  a* =",a_star)
    print("  y* =",y_star)
    # AV总上A量
    S_av=sp.simplify(K*a_star)
    print("  S_av=K·a* =",S_av)
    # 系统总成本
    xA_eq=sp.simplify(K*a_star+y_star)
    J=sp.simplify(xA_eq*(A0+A1*xA_eq)+(N-xA_eq)*(B0+B1*(N-xA_eq)))
    print("  x_A* =",xA_eq)
    print()
    # 存档解
    globals()['A_STAR']=a_star; globals()['Y_STAR']=y_star; globals()['SAV']=S_av

print("="*74)
print("验证: 减堵型(γ<1/2)是否恢复正影响 + 社会最优对比")
print("="*74)
a_star=A_STAR; S_av=SAV
# AV份额 s_eq
s_eq=sp.simplify(S_av/(alpha*N))
print("均衡AV份额 s_eq = S_av/(αN) =",s_eq)
print("  = B1·K/(A1(K-2γ)+B1·K) = B1·K/((A1+B1)K - 2A1γ)")
print()
# 社会最优: 规划者选x_A最小化系统总成本 J=x_A cA + x_B cB
xAv=sp.symbols('x_A',positive=True)
Jsys=xAv*(A0+A1*xAv)+(N-xAv)*(B0+B1*(N-xAv))
xA_soc=sp.solve(sp.diff(Jsys,xAv),xAv)[0]
print("社会最优总流量 x_A^SO =",sp.simplify(xA_soc))
xA_we=sp.simplify((-A0+B0+B1*N)/(A1+B1))
print("用户均衡总流量 x_A^UE =",xA_we)
print()
# AV在共享路的社会最优份额: 规划者希望AV承担多少? 
# 关键: coalition内部化使x_A偏离. 均衡x_A^eq vs x_A^UE
xA_eq=sp.simplify(K*a_star+Y_STAR)
print("均衡总流量 x_A^eq =",xA_eq,"(=UE, 因HDV拉平)")
print()
print("★ coalition影响体现在 AV/HDV 的构成, 及系统成本 vs 纯HDV基准:")
# 基准: 全是HDV(α=0)的Wardrop成本 J_UE
J_UE=sp.simplify(xA_we*(A0+A1*xA_we)+(N-xA_we)*(B0+B1*(N-xA_we)))
J_SO=sp.simplify(xA_soc*(A0+A1*xA_soc)+(N-xA_soc)*(B0+B1*(N-xA_soc)))
print("J_UE(纯Wardrop)=",J_UE)
print("J_SO(社会最优)=",J_SO)
print("PoA gap J_UE-J_SO =",sp.simplify(J_UE-J_SO))
