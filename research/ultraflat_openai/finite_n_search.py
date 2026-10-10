#!/usr/bin/env python3
"""Finite-N numerical Littlewood-polynomial feasibility benchmark.

Constructs actual sign sequences by stochastic single-flip annealing of
a sampled SMOOTH MINIMAX AMPLITUDE objective. NOT the OpenAI Sections 4-6
construction. This is a separate finite-N comparator, useful to establish
measured N / flatness levels without inventing asymptotic N0.

For a sign flip at j, the exact (on a fixed frequency grid) squared
spectrum update is:
 q'_m-q_m = (4/N)*(1-a_j Re(conj(P_m)*exp(2*pi*i*j*m/M)))
 where q_m=|P_m|^2/N-1. No FFT needed per candidate flip.

For the final grid extrema U,V of p(t)=|P(exp(2*pi*i*t))|^2/N,
Bernstein's inequality certifies, in exact arithmetic,
   max p <= U/(1-alpha), min p >= V-alpha*U/(1-alpha)
 where alpha=pi^2*(N-1)^2/(2*M^2)<1.
Returned numbers use floating-point grids; full rigorous numerical
certification would require outward-rounded FFT / interval arithmetic.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


def sampled_polynomial(a, M):
    n=len(a)
    if M<2*n: raise ValueError("Require M>=2N")
    x=np.zeros(M,dtype=float)
    x[:n]=a
    return np.fft.ifft(x)*M  # positive Fourier exponent


def soft_amplitude_objective(q, kappa):
    residual=np.sqrt(np.maximum(1+q,0))-1
    z=np.concatenate((kappa*residual,-kappa*residual))
    return float((logsumexp(z)-math.log(len(z)))/kappa)


def bernstein_grid_bounds(a, grid_factor=32):
    n=len(a)
    M=max(2*n,int(grid_factor*n))
    P=sampled_polynomial(a,M)
    values=(abs(P)**2)/n
    V=float(min(values))
    U=float(max(values))
    alpha=math.pi**2*(n-1)**2/(2*M*M)
    if alpha>=1: raise ValueError("Grid too coarse for Bernstein bound")
    upper=U/(1-alpha)
    lower=max(0,V-alpha*upper)
    return {"N":n,"grid_M":M,"alpha":alpha,"sample_min_square":V,
      "sample_max_square":U,"analytic_min_square_lower":lower,
      "analytic_max_square_upper":upper,
      "sample_relative_error":max(1-math.sqrt(V),math.sqrt(U)-1),
      "analytic_full_circle_relative_error_bound":
       max(1-math.sqrt(lower),math.sqrt(upper)-1),
      "numeric_qualification":"Bernstein inequality exact; FFT/grid extrema use float64 and are NOT interval certified"}


def optimize(n, steps=20000, restarts=4, seed=20261010,
             sharpness=8., temperature=.07, grid_factor=8):
    rng=np.random.default_rng(seed)
    M=grid_factor*n
    if M<2*n: raise ValueError("grid factor must be >=2")
    theta=2*math.pi*np.arange(M)/M
    exponentials=np.exp(1j*np.outer(np.arange(n),theta))
    best_score=float("inf")
    best_a=None
    for rep in range(restarts):
        a=rng.choice([-1,1],size=n).astype(np.int8)
        P=a@exponentials
        q=(abs(P)**2)/n-1
        current=soft_amplitude_objective(q,sharpness)
        for it in range(steps):
            j=int(rng.integers(n))
            E=exponentials[j]
            q_new=q+(4/n)*(1-int(a[j])*np.real(np.conj(P)*E))
            score=soft_amplitude_objective(q_new,sharpness)
            t=temperature*(1-it/steps)**1.5+.004
            if score<=current or rng.random()<math.exp(min(0,(current-score)/t)):
                P-=2*int(a[j])*E
                a[j]*=-1
                q=q_new
                current=score
            if current<best_score:
                best_score=current
                best_a=a.copy()
    certificate=bernstein_grid_bounds(best_a)
    return {"N":n,"seed":seed,"steps_per_restart":steps,"restarts":restarts,
        "soft_max_objective":best_score,
        "signs":"".join('+' if int(x)>0 else '-' for x in best_a),
        "grid_certificate":certificate,
        "status":"finished finite-N sign search; NOT Proposition 4.1 /5.1, NOT theorem N0"}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lengths",type=int,nargs="+",default=[32,64,128])
    ap.add_argument("--steps",type=int,default=20000)
    ap.add_argument("--restarts",type=int,default=4)
    ap.add_argument("--seed",type=int,default=20261010)
    ap.add_argument("--sharpness",type=float,default=8.)
    ap.add_argument("--temperature",type=float,default=.07)
    ap.add_argument("--output",default="")
    a=ap.parse_args()
    data={"method":"actual finite-length sign-sequence search with analytic grid-to-circle envelope",
       "scope":"Independent numerical comparator only, NOT the full original OpenAI proof",
       "results":[optimize(n,steps=a.steps,restarts=a.restarts,
                seed=a.seed+i,sharpness=a.sharpness,temperature=a.temperature)
                for i,n in enumerate(a.lengths)]}
    payload=json.dumps(data,indent=2)
    print(payload)
    if a.output:
        dest=Path(a.output)
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(payload+"\n")
if __name__=="__main__":
    main()
