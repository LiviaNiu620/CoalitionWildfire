"""
Numerical verification of the Regime Index R for the scale-vs-competitive game.

Theory (Regime_Index_推导.md):
  R_0(m_I, gamma, alpha) = m_I / (gamma * [2(1-alpha) - m_I/2])
  Threshold: R_0 = 1  ->  m_I* ≈ 2*gamma*(1-alpha)

Prediction:
  R_0 < 1 (i.e. m_I < 2*gamma*(1-alpha))  =>  dJ/dm_I > 0  (adding I hurts, 反噬)
  R_0 > 1                                  =>  dJ/dm_I < 0  (adding I helps)

Validation strategy:
  1. Scan (m_I, gamma) grid at fixed alpha (0.6, 0.4).
  2. Compute J(m_I) and dJ/dm_I (finite difference) at each point.
  3. Compute R_0 at each point.
  4. Plot heatmap of sign(dJ/dm_I) and overlay R_0 = 1 contour.
  5. Check overlap.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt

from braess_scale_vs_competitive import equilibrium


# ── Theoretical regime index ─────────────────────────────────────────────────

def R_0(m_I, gamma, alpha):
    """Simplified regime index (a priori, no equilibrium needed)."""
    denom = gamma * (2 * (1 - alpha) - m_I / 2)
    # guard against tiny/negative denom
    if denom <= 0:
        return np.inf
    return m_I / denom


def m_I_star(gamma, alpha):
    """Predicted threshold: m_I* = 2*gamma*(1-alpha) (leading-order)."""
    return 2 * gamma * (1 - alpha)


# ── Compute J and dJ/dm_I on a grid ──────────────────────────────────────────

def compute_J_grid(alpha, gammas, m_I_frac_grid,
                   n_runs=5, seed=123, iters=25_000, lam=0.015):
    """
    Returns:
      J[i, j]      = mean equilibrium J at gammas[i], m_I = m_I_frac_grid[j]*alpha
      dJ[i, j]     = numerical dJ/dm_I via central finite difference
      R0[i, j]     = simplified regime index R_0
    """
    rng = np.random.default_rng(seed)
    n_g = len(gammas)
    n_m = len(m_I_frac_grid)
    J = np.zeros((n_g, n_m))

    for i, gamma in enumerate(gammas):
        for j, frac in enumerate(m_I_frac_grid):
            m_I = frac * alpha
            Js = []
            for _ in range(n_runs):
                eq = equilibrium(alpha, m_I, gamma, iters=iters, lam=lam, rng=rng)
                Js.append(eq['J'])
            J[i, j] = np.mean(Js)

    # central finite difference for dJ/dm_I
    dJ = np.zeros_like(J)
    m_I_grid = m_I_frac_grid * alpha
    for j in range(n_m):
        if j == 0:
            dJ[:, j] = (J[:, j+1] - J[:, j]) / (m_I_grid[j+1] - m_I_grid[j])
        elif j == n_m - 1:
            dJ[:, j] = (J[:, j] - J[:, j-1]) / (m_I_grid[j] - m_I_grid[j-1])
        else:
            dJ[:, j] = (J[:, j+1] - J[:, j-1]) / (m_I_grid[j+1] - m_I_grid[j-1])

    # regime index
    R0 = np.zeros_like(J)
    for i, gamma in enumerate(gammas):
        for j, frac in enumerate(m_I_frac_grid):
            m_I = frac * alpha
            R0[i, j] = R_0(m_I, gamma, alpha)

    return J, dJ, R0


# ── Plot heatmap ─────────────────────────────────────────────────────────────

def plot_regime_diagram(alpha, gammas, m_I_frac_grid, J, dJ, R0, filename):
    """
    Two-panel figure:
      (a) heatmap of dJ/dm_I with R_0 = 1 contour and predicted threshold curve.
      (b) heatmap of J with same overlays.
    """
    m_I_grid = m_I_frac_grid * alpha
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ---- Panel (a): dJ/dm_I ------------------------------------------------
    ax = axes[0]
    # symmetric colormap around 0
    vmax = np.percentile(np.abs(dJ), 95)
    im = ax.pcolormesh(m_I_frac_grid, gammas, dJ,
                       cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    plt.colorbar(im, ax=ax, label=r'$\partial J / \partial m_I$')

    # dJ/dm_I = 0 contour (empirical regime boundary)
    cs_emp = ax.contour(m_I_frac_grid, gammas, dJ, levels=[0.0],
                        colors='black', linewidths=2.5, linestyles='-')
    ax.clabel(cs_emp, fmt={0.0: r'$\partial J/\partial m_I=0$ (empirical)'}, fontsize=9)

    # R_0 = 1 contour (analytical prediction)
    # R_0 uses log for stability
    with np.errstate(divide='ignore', invalid='ignore'):
        logR = np.log(np.maximum(R0, 1e-6))
    cs_ana = ax.contour(m_I_frac_grid, gammas, logR, levels=[0.0],
                        colors='lime', linewidths=2.5, linestyles='--')
    ax.clabel(cs_ana, fmt={0.0: r'$\mathcal{R}_0=1$ (theory)'}, fontsize=9)

    # predicted analytical curve m_I = 2*gamma*(1-alpha)
    gamma_dense = np.linspace(gammas[0], gammas[-1], 200)
    m_star_frac = np.array([m_I_star(g, alpha) / alpha for g in gamma_dense])
    m_star_frac = np.clip(m_star_frac, 0, 1)
    ax.plot(m_star_frac, gamma_dense, color='yellow', linewidth=2.5, linestyle=':',
            label=r'$m_I^*=2\gamma(1-\alpha)$')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlabel(r'$m_I / \alpha$   (I 公司占比)')
    ax.set_ylabel(r'$\gamma$   (spite 强度)')
    ax.set_title(rf'(a) $\partial J/\partial m_I$   (α={alpha})' + '\n'
                 r'Blue = 加 I 有益; Red = 加 I 反噬')

    # ---- Panel (b): J itself ----------------------------------------------
    ax = axes[1]
    im2 = ax.pcolormesh(m_I_frac_grid, gammas, J,
                        cmap='viridis', shading='auto')
    plt.colorbar(im2, ax=ax, label=r'$J$ (equilibrium system cost)')
    # overlay empirical and analytical boundaries
    ax.contour(m_I_frac_grid, gammas, dJ, levels=[0.0],
               colors='white', linewidths=2.5, linestyles='-')
    ax.plot(m_star_frac, gamma_dense, color='red', linewidth=2.5, linestyle=':',
            label=r'$m_I^*=2\gamma(1-\alpha)$')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlabel(r'$m_I / \alpha$')
    ax.set_ylabel(r'$\gamma$')
    ax.set_title(rf'(b) System cost $J$   (α={alpha})' + '\n'
                 r'White line = empirical regime boundary')

    plt.tight_layout()
    plt.savefig(filename, dpi=140, bbox_inches='tight')
    plt.close()
    print(f'  saved: {filename}')


def plot_slice(alpha, gammas_slice, m_I_frac_grid, filename,
               n_runs=5, iters=25_000, lam=0.015, seed=456):
    """
    1D slices: for each gamma, plot J vs m_I/alpha with predicted threshold m_I*.
    """
    rng = np.random.default_rng(seed)
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(gammas_slice)))

    for c, gamma in zip(colors, gammas_slice):
        Js = []
        for frac in m_I_frac_grid:
            m_I = frac * alpha
            Js_reps = [equilibrium(alpha, m_I, gamma, iters=iters, lam=lam, rng=rng)['J']
                       for _ in range(n_runs)]
            Js.append(np.mean(Js_reps))
        ax.plot(m_I_frac_grid, Js, 'o-', color=c, label=rf'$\gamma={gamma:.2f}$', ms=4)

        # predicted threshold
        frac_star = min(m_I_star(gamma, alpha) / alpha, 1.0)
        if frac_star > 0:
            ax.axvline(frac_star, color=c, linestyle=':', alpha=0.7)

    ax.set_xlabel(r'$m_I / \alpha$')
    ax.set_ylabel(r'$J$')
    ax.set_title(rf'J vs I-fraction   (α={alpha})' + '\n'
                 r'虚线 = 预测阈值 $m_I^* / \alpha = 2\gamma(1-\alpha)/\alpha$')
    ax.legend(fontsize=9, loc='best')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=140, bbox_inches='tight')
    plt.close()
    print(f'  saved: {filename}')


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print('='*70)
    print('Regime Index R validation on Braess (alpha=0.6)')
    print('='*70)

    # --- 2D scan for heatmap ---
    alpha = 0.6
    gammas = np.linspace(0.15, 1.4, 10)          # 10 gamma values
    m_I_frac_grid = np.linspace(0.05, 1.0, 15)   # 15 m_I fractions

    print(f'\nScanning ({len(gammas)} gammas) x ({len(m_I_frac_grid)} m_I fractions) = {len(gammas)*len(m_I_frac_grid)} points...')
    J, dJ, R0 = compute_J_grid(alpha, gammas, m_I_frac_grid, n_runs=2, iters=12_000, lam=0.02)
    print(f'  J range: [{J.min():.3f}, {J.max():.3f}]')
    print(f'  dJ range: [{dJ.min():.3f}, {dJ.max():.3f}]')

    plot_regime_diagram(alpha, gammas, m_I_frac_grid, J, dJ, R0,
                        filename='regime_diagram_alpha0.6.png')

    # --- 1D slices for direct comparison ---
    print('\nProducing 1D slices for J vs m_I/alpha ...')
    plot_slice(alpha, gammas_slice=[0.25, 0.5, 0.75, 1.0, 1.25],
               m_I_frac_grid=np.linspace(0.02, 1.0, 20),
               filename='J_slices_alpha0.6.png', n_runs=3, iters=15_000)

    # --- Repeat at alpha = 0.4 for comparison ---
    print('\n' + '='*70)
    print('Repeat at alpha=0.4')
    print('='*70)
    alpha = 0.4
    J2, dJ2, R02 = compute_J_grid(alpha, gammas, m_I_frac_grid, n_runs=2, iters=12_000, lam=0.02)
    plot_regime_diagram(alpha, gammas, m_I_frac_grid, J2, dJ2, R02,
                        filename='regime_diagram_alpha0.4.png')
    plot_slice(alpha, gammas_slice=[0.25, 0.5, 0.75, 1.0, 1.25],
               m_I_frac_grid=np.linspace(0.02, 1.0, 20),
               filename='J_slices_alpha0.4.png', n_runs=3, iters=15_000)

    # --- Numerical summary: how well do empirical and analytical boundaries match? ---
    print('\n' + '='*70)
    print('Match summary: predicted vs empirical threshold')
    print('='*70)
    for alpha_i, J_i, dJ_i in [(0.6, J, dJ), (0.4, J2, dJ2)]:
        print(f'\nalpha = {alpha_i}')
        print(f'  {"gamma":>6} {"predicted m_I*/α":>18} {"empirical m_I*/α (sign flip)":>32}')
        for i, gamma in enumerate(gammas):
            pred = min(m_I_star(gamma, alpha_i) / alpha_i, 1.0)
            # find empirical sign-flip point
            sign = np.sign(dJ_i[i, :])
            emp = None
            for j in range(len(m_I_frac_grid)-1):
                if sign[j] > 0 and sign[j+1] <= 0:
                    # linear interp for zero crossing
                    emp = m_I_frac_grid[j] + (m_I_frac_grid[j+1] - m_I_frac_grid[j]) * \
                          dJ_i[i, j] / (dJ_i[i, j] - dJ_i[i, j+1])
                    break
            emp_str = f'{emp:.3f}' if emp is not None else '(no flip)'
            print(f'  {gamma:6.2f} {pred:18.3f} {emp_str:>32}')
