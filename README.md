# Entangled-noise certificates

Exact checks with hardcoded data; no TeX or input files required.

```sh
python3 verify_certificates.py --json
python3 -m pip install sympy==1.14.0
python3 verify_fs_ghz.py
```

The first script uses only the standard library. Both print results without writing files. Run without `-O`. Finite checks do not replace the general proofs.
