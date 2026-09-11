"""Symbolic checks used in the review (sympy):
1) Theorem 2 bottleneck derivative dJ/dSigma (residual should be 0);
2) frozen-slope coalition model: VI Jacobian is symmetric (=> exact potential game);
3) true atomic splittable internalisation (state-dependent slope c'(x)): Jacobian is NOT symmetric."""
import sympy as sp
H, S, x, t1, N, a = sp.symbols('H Sigma x t1 N a', positive=True)
x1 = (H + S*x)/(1 + S)
J = N*a + t1*(x1**2 - x*x1)
print('1) bottleneck residual:', sp.simplify(sp.diff(J, S) - t1*(x - H)*(2*H + (S - 1)*x)/(1 + S)**3))
f1, f2, h, bb, tau1, tau2, g = sp.symbols('f1 f2 h bbar tau1 tau2 gamma', positive=True)
X = f1 + f2 + h
c = 1 + sp.Rational(15, 100)*X**4                        # BPR with t0=kappa=1
D = sp.diag(sp.sqrt(tau1), sp.sqrt(tau2))
Om = bb*D*((1 - g)*sp.eye(2) + g*sp.ones(2, 2))*D        # governance matrix, frozen slope bbar
F = sp.Matrix([c, c, c]) + sp.Matrix.vstack(Om*sp.Matrix([f1, f2]), sp.Matrix([0]))
Jac = F.jacobian([f1, f2, h])
print('2) frozen-slope Jacobian asymmetry:', sp.simplify(Jac - Jac.T))
Ft = sp.Matrix([c + sp.diff(c, f1)*f1, c + sp.diff(c, f2)*f2, c])
Jt = Ft.jacobian([f1, f2, h])
print('3) true-atomic Jacobian asymmetry (0,1):', sp.factor(sp.simplify((Jt - Jt.T)[0, 1])))
