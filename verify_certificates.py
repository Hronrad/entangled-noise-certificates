#!/usr/bin/env python3
"""Exact rational certificate checks; standard library only."""
import argparse
import itertools
import json
from fractions import Fraction as F


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def pt(a, dim_b):
    dim_a = len(a)//dim_b
    return [[a[dim_b*(i//dim_b)+(j % dim_b)]
               [dim_b*(j//dim_b)+(i % dim_b)]
             for j in range(dim_a*dim_b)]
            for i in range(dim_a*dim_b)]


def determinant(a):
    a = [row[:] for row in a]
    value = F(1)
    for i in range(len(a)):
        p = next((j for j in range(i, len(a)) if a[j][i]), None)
        if p is None:
            return F(0)
        if p != i:
            a[i], a[p] = a[p], a[i]
            value = -value
        pivot = a[i][i]
        value *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return value


def psd(a):
    n = len(a)
    assert a == transpose(a)
    return all(determinant([[a[i][j] for j in ids] for i in ids]) >= 0
               for k in range(1, n+1)
               for ids in itertools.combinations(range(n), k))


def ldlt(a):
    n = len(a)
    assert a == transpose(a)
    lower = eye(n)
    pivots = []
    for i in range(n):
        pivot = a[i][i]-sum((lower[i][k]**2*pivots[k]
                            for k in range(i)), F(0))
        assert pivot > 0
        pivots.append(pivot)
        for j in range(i+1, n):
            lower[j][i] = (a[j][i]-sum(
                (lower[j][k]*lower[i][k]*pivots[k]
                 for k in range(i)), F(0)))/pivot
    diagonal = [[pivots[i] if i == j else F(0) for j in range(n)]
                for i in range(n)]
    assert multiply(multiply(lower, diagonal), transpose(lower)) == a
    return pivots


def outer(v):
    return [[x*y for y in v] for x in v]


def x_witness(a, b, c, d, w, z):
    rho = [[a,0,0,w],[0,b,z,0],[0,z,c,0],[w,0,0,d]]
    rho = [[F(x) for x in row] for row in rho]
    if w*w > b*c:
        i, j, x, y, t = 1, 2, b, c, w
    elif z*z > a*d:
        i, j, x, y, t = 0, 3, a, d, z
    else:
        return rho, scale(F(0), eye(4)), scale(F(0), eye(4))
    if t >= max(x, y):
        u, v = t-x, t-y
    elif x > t:
        u, v = F(0), t*t/x-y
    else:
        u, v = t*t/y-x, F(0)
    noise = scale(F(0), eye(4))
    noise[i][i], noise[j][j] = u, v
    if x > y:
        i, j, x, y = j, i, y, x
    mu = min(t/y, F(1)) if y else F(1)
    vector = [F(0)]*4
    vector[i], vector[j] = F(1), -mu
    return rho, noise, outer(vector)


def verify():
    result = {"arithmetic": "exact fractions; no floating-point positivity"}
    rho = scale(F(1,40), [
        [2,0,0,0,0,0], [0,1,0,-4,0,0], [0,0,13,0,-5,0],
        [0,-4,0,17,0,0], [0,0,-5,0,2,0], [0,0,0,0,0,5]])
    assert trace(rho) == 1
    rho_pivots = ldlt(rho)

    u = list(map(F, [0,2,0,1,0,0]))
    v1 = list(map(F, [1,0,0,0,1,0]))
    v2 = list(map(F, [0,1,0,0,0,1]))
    yg = scale(F(1,40), outer(u))
    w = scale(F(2,3), add(outer(v1), outer(v2)))
    z = pt(w, 3)
    target = add(rho, yg)
    assert psd(yg) and psd(w)
    assert min(ldlt(target)) > 0
    assert psd(pt(target,3))
    assert not psd(pt(yg,3)), "the displayed noise is entangled"
    assert multiply(w, pt(target,3)) == scale(F(0), eye(6))
    # Independent block check.
    order = [0,4,1,5,2,3]
    blocks = [[40*pt(target,3)[i][j] for j in order] for i in order]
    assert blocks == [[2,-2,0,0,0,0],[-2,2,0,0,0,0],
                      [0,0,5,-5,0,0],[0,0,-5,5,0,0],
                      [0,0,0,0,13,0],[0,0,0,0,0,18]]
    upper = trace(yg)
    lower = -trace(multiply(z,rho))
    gap = lower-upper
    assert upper == F(1,8)
    assert lower == F(2,15)
    assert gap == F(1,120)
    assert trace(multiply(z,outer(u)))/sum(x*x for x in u) == F(16,15)
    assert trace(multiply(add(z,scale(-1,eye(6))),yg)) == gap
    assert not psd(add(eye(6),scale(-1,z)))
    result["qubit_qutrit"] = {
        "Rg_upper": str(upper), "Rs_lower": str(lower), "gap_lower": str(gap),
        "rho_LDLT_pivots": list(map(str,rho_pivots)), "rho_rank": 6,
        "target_pt_integer_spectrum": [0,0,4,10,13,18],
        "witness_noise_expectation": "16/15",
        "universal_witness_bound": "Analytical complex-product norm proof in the supplement, not a finite-sampling assertion.",
        "optimality": "The two robustness optima and noise optimality are not asserted."}

    epsilon = F(1,200)
    fullrank_rho = add(scale(1-epsilon,rho),scale(epsilon/6,eye(6)))
    fullrank_noise = scale(1-epsilon,yg)
    assert min(ldlt(fullrank_rho)) > 0
    assert psd(fullrank_noise)
    assert psd(pt(add(fullrank_rho,fullrank_noise),3))
    witness_trace = trace(z)
    assert witness_trace == trace(w) == F(8,3)
    fullrank_gap = -trace(multiply(z,fullrank_rho))-trace(fullrank_noise)
    assert fullrank_gap == (1-epsilon)*gap-epsilon*witness_trace/6
    assert fullrank_gap == (3-163*epsilon)/360
    assert fullrank_gap == F(437,72000) > F(606,100000)
    result["fullrank_stability"] = {"epsilon":str(epsilon),
        "state_eigenvalue_lower_bound":"1/1200",
        "witness_trace":str(witness_trace),"gap_lower":str(fullrank_gap),
        "certified_positive_epsilon_range":"0 <= epsilon < 3/163"}

    ex = [[F(3,8),0,0,F(1,4)], [0,F(1,5),F(1,100),0],
          [0,F(1,100),F(1,20),0], [F(1,4),0,0,F(3,8)]]
    assert trace(ex) == 1 and min(ldlt(ex)) > 0
    noise = scale(F(0), eye(4))
    noise[1][1], noise[2][2] = F(1,20), F(1,5)
    witness = outer([F(0),F(1),F(-1),F(0)])
    assert psd(noise) and psd(pt(noise,2))
    assert psd(add(ex,noise)) and psd(pt(add(ex,noise),2))
    assert psd(witness) and psd(add(eye(4),scale(-1,pt(witness,2))))
    value = trace(noise)
    assert value == -trace(multiply(witness,pt(ex,2))) == F(1,4)
    # Squared Wootters numbers.
    spin = [[F(x) for x in row] for row in
            [[0,0,0,-1],[0,0,1,0],[0,1,0,0],[-1,0,0,0]]]
    product = multiply(ex, multiply(multiply(spin,ex),spin))
    lambdas = [F(5,8),F(1,8),F(11,100),F(9,100)]
    assert len(set(lambdas)) == 4 and min(lambdas) > 0
    for lam in lambdas:
        assert determinant(add(product,scale(-lam*lam,eye(4)))) == 0
    norms = [F(1),F(1),F(5,4),F(5,4)]
    concurrence = lambdas[0]-sum(lambdas[1:])
    predicted = concurrence*min(sum(pair)
                 for pair in itertools.combinations(norms[1:],2))/2
    assert predicted == F(27,80) > value
    result["JMR_counterexample"] = {"old_expression": str(predicted),
                                   "certified_value": str(value)}

    old_noise = scale(F(1,80), [[6,0,0,-6],[0,12,6,0],
                               [0,6,3,0],[-6,0,0,6]])
    assert trace(old_noise) == predicted
    assert psd(old_noise) and psd(pt(old_noise,2))
    assert psd(pt(add(ex,old_noise),2))
    # Endpoint feasibility extends by convexity.
    for tau in [F(0),F(1,1000000),F(1,3),F(1,2),F(1)]:
        interpolated = add(scale(1-tau,old_noise),scale(tau,noise))
        assert psd(interpolated) and psd(pt(interpolated,2))
        assert psd(pt(add(ex,interpolated),2))
        assert trace(interpolated) == F(27,80)-F(7,80)*tau
    result["JMR_feasible_descent"] = {
        "endpoints_PPT":"PASS for noise and target",
        "cost":"27/80 - (7/80)*tau", "strict_cost_derivative":"-7/80",
        "whole_interval":"feasible by exact endpoints and separable-cone convexity"}

    cases = branches = 0
    rank_histogram = {}
    for aa,bb,cc,dd in itertools.product(range(4),repeat=4):
        norm = aa+bb+cc+dd
        if not norm:
            continue
        for ww,zz in itertools.product(range(4),repeat=2):
            if ww*ww > aa*dd or zz*zz > bb*cc:
                continue
            params = [F(v,norm) for v in (aa,bb,cc,dd,ww,zz)]
            rho,noise,witness = x_witness(*params)
            assert psd(pt(add(rho,noise),2))
            assert min(noise[i][i] for i in range(4)) >= 0
            assert psd(add(eye(4),scale(-1,pt(witness,2))))
            assert trace(noise) == -trace(multiply(witness,pt(rho,2)))
            assert psd(witness)
            cases += 1
            branches += trace(noise) > 0
            rank = sum(2 if x*y-t*t > 0 else int(x+y > 0)
                       for x,y,t in [(aa,dd,ww),(bb,cc,zz)])
            rank_histogram[rank] = rank_histogram.get(rank,0)+1
    assert set(rank_histogram) == {1,2,3,4}
    result["X_state_grid"] = {"cases": cases, "entangled": branches,
                              "ranks": rank_histogram,
                              "checks": "exact primal/dual matching; includes zero and branch boundaries"}

    figure_examples = []
    for coherence, weights in [(F(1,5),(F(1,14),F(0))),
                               (F(3,10),(F(1,5),F(1,15)))]:
        rho, noise, witness = x_witness(F(1,3),F(1,10),F(7,30),
                                         F(1,3),coherence,F(0))
        assert psd(rho) and trace(rho) == 1
        assert (noise[1][1],noise[2][2]) == weights
        assert psd(pt(add(rho,noise),2))
        assert psd(add(eye(4),scale(-1,pt(witness,2))))
        assert trace(noise) == -trace(multiply(witness,pt(rho,2)))
        figure_examples.append({"w":str(coherence),
                                "weights":list(map(str,weights))})
    result["BSA_supplement_figure_mixed_inputs"] = figure_examples

    phase_cases = phase_entangled = 0
    for delta, kappa in itertools.product([F(-1),F(-1,2),F(0),F(1,2),F(1)],
                                         [F(0),F(1,4),F(1,2),F(3,4),F(1)]):
        rho, noise, witness = x_witness(F(1,3),(1-delta)/6,(1+delta)/6,
                                        F(1,3),kappa/3,F(0))
        assert psd(rho) and trace(rho)==1
        assert psd(pt(add(rho,noise),2))
        assert trace(noise)==-trace(multiply(witness,pt(rho,2)))
        entangled=4*kappa*kappa>1-delta*delta
        assert (trace(noise)>0)==entangled
        if entangled:
            terms=sum(noise[i][i]>0 for i in range(4))
            assert terms==(2 if 2*kappa>1+abs(delta) else 1)
            q=noise[1][1]/trace(noise)
            if delta==0: assert q==F(1,2)
            assert (q-F(1,2))*delta>=0
            phase_entangled+=1
        phase_cases+=1
    result["supplement_X_preparation_figure"]={"exact_slice_cases":phase_cases,
        "entangled_cases":phase_entangled,"checks":"PPT and one/two-product boundaries, representative weights and exact certificate matching; minimality itself is proved analytically in SI"}

    result["status"] = "PASS"
    return result


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O; assertions perform the checks.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",action="store_true")
    args = parser.parse_args()
    report = verify()
    if args.json:
        print(json.dumps(report,ensure_ascii=False,indent=2))
    else:
        print("PASS: exact manuscript certificates and X-state grid")
        print(json.dumps({k:v for k,v in report.items()
                          if k != "qubit_qutrit"},
                         ensure_ascii=False,indent=2))
        print("2x3 certified gap =", report["qubit_qutrit"]["gap_lower"])
