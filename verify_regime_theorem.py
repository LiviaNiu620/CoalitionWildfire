"""Dependency-free symbolic check for the canonical regime theorem."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

# A polynomial is indexed by powers of (sigma, h, x_ue).
Poly = dict[tuple[int, int, int], Fraction]


def const(value: int) -> Poly:
    return {(0, 0, 0): Fraction(value)}


def var(axis: int) -> Poly:
    exponent = [0, 0, 0]
    exponent[axis] = 1
    return {tuple(exponent): Fraction(1)}


def add(left: Poly, right: Poly, scale: int = 1) -> Poly:
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, Fraction(0)) + scale * coefficient
        if result[exponent] == 0:
            del result[exponent]
    return result


def mul(left: Poly, right: Poly) -> Poly:
    result: Poly = {}
    for a, ca in left.items():
        for b, cb in right.items():
            exponent = tuple(a[i] + b[i] for i in range(3))
            result[exponent] = result.get(exponent, Fraction(0)) + ca * cb
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def derivative_sigma(poly: Poly) -> Poly:
    result: Poly = {}
    for (power_sigma, power_h, power_xue), coefficient in poly.items():
        if power_sigma:
            exponent = (power_sigma - 1, power_h, power_xue)
            result[exponent] = coefficient * power_sigma
    return result


def main() -> None:
    sigma, h, x_ue = var(0), var(1), var(2)
    denominator = add(const(1), sigma)
    numerator_x = add(h, mul(sigma, x_ue))

    # Up to the positive factor t1 and constant N*a,
    # J = x^2 - x_ue*x = P / denominator^2.
    p = add(mul(numerator_x, numerator_x), mul(mul(x_ue, numerator_x), denominator), scale=-1)
    derivative_numerator = add(mul(derivative_sigma(p), denominator), p, scale=-2)

    expected = mul(
        add(x_ue, h, scale=-1),
        add(add(mul(const(2), h), mul(sigma, x_ue)), x_ue, scale=-1),
    )
    if derivative_numerator != expected:
        raise AssertionError(
            f"regime derivative identity failed: {derivative_numerator} != {expected}"
        )

    output = {
        "derivative": "t1*(x_ue-h)*(2*h+(sigma-1)*x_ue)/(1+sigma)^3",
        "threshold": "1-2*h/x_ue",
        "verification": "exact coefficient equality after clearing (1+sigma)^3",
        "assumptions": ["0 <= h < x_ue", "sigma >= 0", "t1 > 0", "fixed active support"],
        "regimes": {
            "screened": "HDVs remain interior and pin x=x_ue; organization is neutral",
            "centralization": "h >= x_ue/2, so dJ/dsigma >= 0 and lower effective response is preferred",
            "interior_target": "h < x_ue/2, so J decreases below sigma_star and increases above sigma_star",
        },
    }
    Path("regime_theorem_verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
