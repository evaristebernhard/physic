#!/usr/bin/env python3
"""Finite-dimensional numerical test of OpenAI Prop. 4.1 low-J balancing.

The target equation is
   |s.gamma| Prod_{j=1}^J |s.W_j| = tau |b_s|^2 (up to exp(+-kappa)).
The scalar tau is eliminated analytically: the best log-error kappa for
given vectors is one-half the range of
   r_s=Sum_j log|s.W_j| - 2 log|b_s|.

A least-squares solver explores candidate directions. This is NOT a global
certificate. Unlike the heuristic d=3 fits, the d=2 cubic factorization
provides an exact algebraic construction at J=2.
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def original_fourier(depth):
    p = {}
    for level in range(depth):
        if not p:
            p = {(-1,): Fraction(1,4), (1,): Fraction(1,4)}
            continue
        q = {}
        for s,b in p.items():
            for t,c in p.items():
                u = tuple(x+y for x,y in zip(s,t))
                q[u] = q.get(u, Fraction(0))+b*c
        z = (0,)*len(next(iter(p)))
        q[z] = q.get(z, Fraction(0))-1
        new = {s+(0,): b for s,b in p.items()}
        for s,b in q.items():
            if b:
                new[s+(-1,)] = -b/4
                new[s+(1,)] = -b/4
        p = {s:b for s,b in new.items() if b}
    return p


def balanced_error(S, target, vectors):
    v = S @ vectors.T
    if min(abs(v.flat))<1e-12:
        return None
    r = np.log(abs(v)).sum(axis=1)-target
    return (float((r.max()-r.min())/2),
            float(np.std(r)),
            float(np.exp((r.max()+r.min())/2)))


def exact_level_two():
    # A scalar multiple of the complete product of its three real factors:
    # Q(x,y)= x^3+(49/256)x^2 y-(2047/512)x y^2-(49/64)y^3.
    # Exact checks: Q(1,0)=1, Q(0,1)=-49/64,
    # Q(2,1)=Q(2,-1)=1/256.
    coeffs = [-49/64, -2047/512, 49/256, 1.]
    roots = np.roots(coeffs)
    if max(abs(roots.imag))>1e-8:
        raise ArithmeticError("The source cubic unexpectedly has non-real roots")
    W = np.array([[1.,-1/float(r)] for r in roots.real])
    p = original_fourier(2)
    items = [(s,b) for s,b in p.items() if s>tuple(-x for x in s)]
    S = np.array([s for s,b in items], dtype=float)
    values = np.prod(S@W.T,axis=1)
    ideal = np.array([16*float(b*b) for s,b in items])
    return {"J":2,"D":len(items),"frequency_vectors":S.astype(int).tolist(),
            "root_y_over_x":roots.real.tolist(),
            "factor_vectors":W.tolist(),
            "max_relative_error_float":float(max(abs(abs(values)/ideal-1))),
            "exact_statement":"the rational cubic matches 16*|b_s|^2 at every original d=2 frequency, and has three distinct real roots"}


def fit(depth, J, restarts, max_nfev, seed):
    p = original_fourier(depth)
    items = [(s,b) for s,b in p.items() if s>tuple(-x for x in s)]
    S = np.array([s for s,b in items], dtype=float)
    target = np.array([2*math.log(abs(float(b))) for s,b in items])
    d = depth
    rng = np.random.default_rng(seed)
    best = None
    for k in range(restarts):
        v = rng.normal(size=(J+1,d))
        v/=np.linalg.norm(v,axis=1,keepdims=True)
        def residual_and_jac(z):
            W = z.reshape(J+1,d)
            dots = S@W.T
            safe = np.where(abs(dots)>1e-12,dots,
                            np.where(dots>=0,1e-12,-1e-12))
            r = np.log(abs(safe)).sum(axis=1)-target
            jac = (S[:,None,:]/safe[:,:,None]).reshape(len(S),-1)
            return r-r.mean(),jac-jac.mean(axis=0)
        solve = least_squares(lambda z: residual_and_jac(z)[0],v.ravel(),
             jac=lambda z: residual_and_jac(z)[1],max_nfev=max_nfev,
             method="trf")
        W = solve.x.reshape(J+1,d)
        stats = balanced_error(S,target,W)
        if stats is None: continue
        kappa,rms,tau = stats
        record = {"J":J,"kappa":kappa,"rms_log_residual":rms,
                  "tau_fit":tau,"vectors":W.tolist(),
                  "restart":k,"solver_calls":int(solve.nfev)}
        if best is None or kappa<best["kappa"]:
            best=record
    return {"d":depth,"D":len(items),"J":J,"best":best,
            "qualification":"local least-squares candidate, not a proof of optimality or feasibility at target delta"}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--depth",type=int,default=3,choices=[2,3])
    ap.add_argument("--directions",nargs="+",type=int,default=[1,2,3,5,8])
    ap.add_argument("--restarts",type=int,default=5)
    ap.add_argument("--max-nfev",type=int,default=150)
    ap.add_argument("--seed",type=int,default=20261010)
    ap.add_argument("--output",default="")
    a=ap.parse_args()
    data={"method":"low-J balancing of original OpenAI Fourier coefficients",
          "equation":"kappa=1/2(range_s(Sum_j log|s.W_j|-2log|b_s|))",
          "analytic_d2":exact_level_two(),
          "fits":[fit(a.depth,J,a.restarts,a.max_nfev,a.seed+J) for J in a.directions],
          "status":"analytic J=2,d=2 feasible; numerical higher-d tests only; Sections 5-6 not run"}
    blob=json.dumps(data,indent=2)
    print(blob)
    if a.output:
        dest=Path(a.output)
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(blob+"\n")
if __name__=="__main__":
    main()
