#!/usr/bin/env python3
"""Exact GHZ certificate checks for Supplemental Material Sec. X.

Symbolic regression for n = 2,...,6 with sympy. Finite character checks do
not replace the all-n product-vector proof. No optimizer, no floating-point
arithmetic, no network access.
"""
import hashlib
import itertools
import json
import re
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent


def find_source(name):
    for candidate in (ROOT / name, ROOT.parent / name, Path.cwd() / name):
        if candidate.is_file():
            return candidate
    raise SystemExit("Cannot find " + name + "; run the script next to the manuscript source.")


SUPP_NAME = "minimal-dimension-entangled-noise-advantage-supplement.tex"
COMBINED_NAME = "minimal-dimension-entangled-noise-advantage-complete.tex"


def load_supplement():
    for candidate in (ROOT / SUPP_NAME, ROOT.parent / SUPP_NAME, Path.cwd() / SUPP_NAME):
        if candidate.is_file():
            return candidate, candidate.read_text(), "supplement file"
    try:
        combined = find_source(COMBINED_NAME)
    except SystemExit:
        fixture = ROOT / "data/ghz-equations.tex"
        return fixture, fixture.read_text(), "bundled equation fixtures (not full manuscript)"
    text = combined.read_text()
    if text.count("\\onecolumngrid") != 1 or text.count("\\end{document}") != 1:
        raise SystemExit("Unexpected document layout in " + COMBINED_NAME)
    tail = text.split("\\onecolumngrid", 1)[1].split("\\end{document}", 1)[0]
    tail = re.sub(r"\\label\{sm:", r"\\label{", tail)
    tail = re.sub(r"\\(eqref|ref|pageref|autoref)\{sm:", r"\\\1{", tail)
    return combined, tail, "single-file source"


SOURCE, SOURCE_TEXT, SOURCE_LAYOUT = load_supplement()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
        # Only endpoint coherences carry a minus sign for chi=pi.
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
    source = SOURCE_TEXT
    labels = ["eq:fs-ghz-witnesses", "eq:fs-ghz-values",
              "eq:fs-noisy-ghz-values", "eq:fs-noisy-primal", "eq:fs-noisy-free"]
    for label in labels:
        assert r"\label{" + label + "}" in source
    report = {"status": "PASS", "sympy_version": sp.__version__,
              "source_layout": SOURCE_LAYOUT,
              "source_sha256": sha(SOURCE), "script_sha256": sha(Path(__file__)),
              "source_labels": labels,
              "scope": "Exact regression; universal inequalities and domain signs are analytic SM arguments.",
              "universal_scalar_identities": universal_checks(),
              "ghz": [ghz_checks(n) for n in range(2, 7)]}
    if (ROOT / "audit").is_dir():
        output = ROOT / "audit/results/fs_ghz_exact.json"
    else:
        output = Path.cwd() / "fs_ghz_exact.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
