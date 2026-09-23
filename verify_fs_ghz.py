#!/usr/bin/env python3
"""Exact GHZ checks for n = 2,...,6 and symbolic scalar identities."""
import itertools
import json

import sympy as sp

def zero(matrix):
    assert all(sp.cancel(e) == 0 for e in matrix)


def universal_checks():
    D, p, r = sp.symbols("D p r", real=True)
    beta, alpha = D / 4, D / (2 * (D - 2))
    identities = {
        "lower_product_coefficient": alpha * (D - 2) - 2 * beta,
        "upper_product_bound": alpha + 2 * (beta - alpha) / D - 1,
        "standard_identity_trace": alpha * (D - 2) / D - sp.Rational(1, 2),
        "generalized_identity_trace": alpha / beta * (D - 2) / D - 2 / D,
        "white_noise_ratio": beta * (1 - p - 2 * p / D) - (beta * (1 - p) - p / 2),
        "threshold": (1 - p - 2 * p / D).subs(p, D / (D + 2)),
        "upper_margin": 1 - alpha - 2 * (beta - alpha) * r
                       - 2 * (beta - alpha) * (1 / D - r),
    }
    for name, value in identities.items():
        assert sp.cancel(value) == 0, name
    return {name: "exact_zero" for name in identities}


def ghz_checks(n):
    D = 2**n
    eye, E, C = sp.eye(D), sp.zeros(D), sp.zeros(D)
    E[0, 0] = E[-1, -1] = C[0, -1] = C[-1, 0] = 1
    Pi = eye - E
    G, Gminus = (E + C) / 2, (E - C) / 2
    plus, minus = (eye + C) / D, (eye - C) / D
    beta, alpha = sp.Rational(D, 4), sp.Rational(D, 2 * (D - 2))
    Zs = alpha * Pi - beta * C
    Zg = Zs / beta
    assert sp.trace(G) == sp.trace(plus) == sp.trace(minus) == 1
    zero(G + Gminus - E)
    zero(G + beta * minus - beta * plus - E / 2)
    assert -sp.trace(Zs * G) == beta and -sp.trace(Zg * G) == 1
    assert (eye - Zg).is_positive_semidefinite is True
    for x, y in itertools.product(range(D), repeat=2):
        xb = [(x >> (n - 1 - j)) & 1 for j in range(n)]
        yb = [(y >> (n - 1 - j)) & 1 for j in range(n)]
        survives = all((xb[j] - yb[j] - xb[-1] + yb[-1]) % 3 == 0
                       for j in range(n - 1))
        expected = x == y or {x, y} == {0, D - 1}
        assert survives == expected
        # Endpoint phase at chi=pi.
        coefficient = (sp.Integer(-1)**(xb[-1] - yb[-1])) if survives else 0
        assert coefficient == D * minus[x, y]
    p = sp.symbols("p", real=True)
    rho = (1 - p) * G + p * eye / D
    s, t = 1 - p - 2 * p / D, beta * (1 - p) - p / 2
    zero(rho + s * Gminus - p * plus - (1 - p - p / D) * E)
    zero(rho + t * minus - (p + t) * plus - (1 - p) * E / 2)
    assert sp.cancel(-sp.trace(Zs * rho) - t) == 0
    assert sp.cancel(-sp.trace(Zg * rho) - s) == 0
    zero(rho - D * (1 - p) / 2 * plus - p / D * E
         - (p / D - (1 - p) / 2) * Pi)
    return {"n": n, "dimension": D, "Rg": "1", "Rs": str(beta),
            "threshold": str(sp.Rational(D, D + 2)),
            "phase_characters_and_primal_dual_identities": "PASS"}


def main():
    report = {"status": "PASS", "sympy_version": sp.__version__,
              "scope": "Exact regression; universal inequalities and domain signs are analytic SM arguments.",
              "universal_scalar_identities": universal_checks(),
              "ghz": [ghz_checks(n) for n in range(2, 7)]}
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O; assertions perform the checks.")
    main()
