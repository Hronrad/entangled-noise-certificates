# Entangled-noise certificates

Exact-arithmetic verification scripts and minimal equation inputs for entanglement robustness certificates. This repository contains no full manuscript, manuscript PDF, figures, or private review files.

## Run

Python 3.11 or later is recommended. The rational-certificate checker uses only the standard library; the GHZ checker uses SymPy.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify_certificates.py --json
python verify_fs_ghz.py
```

Both scripts print JSON results. The GHZ script also writes `fs_ghz_exact.json` in the working directory. Do not run with Python's `-O` option: assertions are part of these checks.

## Inputs and scope

- `data/certificate-main.tex` and `data/certificate-supplement.tex` are exact, minimal equation excerpts. The rational checker parses the labelled matrices/vectors and checks primal feasibility, rational gaps, a counterexample to an earlier expression, and finite X-state grids.
- `data/ghz-equations.tex` supplies the referenced GHZ equations. The GHZ checker checks that its expected labels exist; its mathematical expressions are implemented in Python, **not parsed from the TeX formulas**. It checks symbolic scalar identities and exact finite cases n = 2,...,6.
- The snippets are data, not independently compilable TeX documents. They omit the surrounding proofs. A successful run is not a proof of all-state or all-n claims, witness bounds over every product state, or the manuscript's completeness.

For private manuscript-source checks, run either script with the current directory containing `minimal-dimension-entangled-noise-advantage.tex` and `minimal-dimension-entangled-noise-advantage-supplement.tex`, or their merged `minimal-dimension-entangled-noise-advantage-complete.tex`. Those sources are preferred when found; otherwise the bundled fixtures are used. The JSON records the selected source layout and SHA-256 hashes. Keep full manuscripts out of this public repository.

## Version and rights

The accompanying paper cites the fixed tag `v2026.09.23` and its commit. The mathematical checking logic is unchanged from the source scripts; this distribution only adds a fallback to bundled equation inputs. Public availability does not by itself grant an open-source license. No license has been selected by the rights holders.
