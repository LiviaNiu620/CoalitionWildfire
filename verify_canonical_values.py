"""High-precision source of all canonical values reported in the paper."""

import mpmath as mp


mp.mp.dps = 50


def beta(m):
    return 1 - mp.exp(-2 * mp.mpf(m))


def equilibrium(masses, interior, h=0):
    sigma = sum(1 / beta(masses[k]) for k in interior)
    x1 = (mp.mpf(h) + sigma) / (1 + sigma)
    mu = 1 - x1
    flows = [mu / beta(masses[k]) for k in interior]
    rho = 4 * x1 * (1 - x1)
    return sigma, x1, mu, flows, rho


def print_case(name, masses, interior, h=0):
    sigma, x1, mu, flows, rho = equilibrium(masses, interior, h)
    print(name)
    print("  beta  =", [mp.nstr(beta(m), 20) for m in masses])
    print("  Sigma =", mp.nstr(sigma, 25))
    print("  x1    =", mp.nstr(x1, 25))
    print("  mu    =", mp.nstr(mu, 25))
    print("  flows =", [mp.nstr(f, 25) for f in flows])
    print("  check =", mp.nstr(sum(flows) + h, 25))
    print("  rho   =", mp.nstr(rho, 25))


if __name__ == "__main__":
    print_case("(0.5, 0.5)", [0.5, 0.5], [0, 1])
    print_case("(0.55, 0.45)", [0.55, 0.45], [0, 1])
    print_case("(0.9, 0.1)", [0.9, 0.1], [0], mp.mpf("0.1"))
    print_case("monopoly", [1.0], [0])
    print_case("4 equal", [0.25] * 4, [0, 1, 2, 3])
