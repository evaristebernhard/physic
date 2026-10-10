#!/usr/bin/env python3
"""FFT acceleration of the UNMODIFIED OpenAI p_d Fourier recursion.

At depth d, 1-p_d(y)^2 is obtained from the multidimensional
convolution of the Fourier array with itself. Construct the next layer:
  p_{d+1}=p_d + (1-p_d^2)cos(2*pi*y_{d+1})/2.
This script computes depth 6 from the EXACT dyadic depth-5 coefficients,
then optionally prunes ONLY depth-6 coefficients with an explicit L1
loss budget. It does not attempt Section 4 coefficient balancing or N0.

The convolution itself uses float64 FFT and is *numerical*, not interval
certified. Algebraic truncation bounds are valid in exact arithmetic;
their numerical application needs a separate FFT roundoff bound.
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

try:
    from .exact_pruning import step
except ImportError:
    from exact_pruning import step


def original_poly_exact(d):
    poly={}
    exponent=0
    for i in range(d):
        poly,exponent=step(poly,exponent,max_pairs=2_000_000)
    return poly,exponent


def convolve_dense(parent, exponent, limit=3_000_000):
    d=len(next(iter(parent)))
    x=np.asarray(list(parent),dtype=int)
    lo=x.min(axis=0)
    hi=x.max(axis=0)
    shape=tuple(int(h-l+1) for l,h in zip(lo,hi))
    full=tuple(2*s-1 for s in shape)
    count=math.prod(full)
    if count>limit:
        raise RuntimeError(
            f"FFT_CELL_BUDGET: dense convolution requires {count:,} cells > {limit:,}")
    old=np.zeros(shape,dtype=np.float64)
    denom=float(1<<exponent)
    for s,b in parent.items():
        old[tuple(int(v-l) for v,l in zip(s,lo))]=float(b)/denom
    t0=time.perf_counter()
    convolution=fftconvolve(old,old,mode="full")
    fft_seconds=time.perf_counter()-t0
    # Real Fourier symmetry is invariant under the true convolution.
    convolution=(convolution+np.flip(convolution))/2
    # delta_0 - p_d^2 in the new-frequency +/-1 planes.
    zero=tuple(int(-2*l) for l in lo)
    convolution[zero]-=1
    new=-convolution/4
    M_old=float(sum((b/denom)**2 for b in parent.values()))
    M_new=M_old+2*float(np.sum(new*new))
    return new,{"d":d+1,"old_fourier_support":len(parent),
                "dense_parent_shape":list(shape),
                "dense_convolution_shape":list(full),
                "dense_cells":int(count),
                "fft_convolution_seconds":fft_seconds,
                "original_energy_estimate":M_new,
                "old_energy":M_old}


def prune_new_symmetric(new, old_energy, old_terms, budget):
    """Prune both new-frequency slabs, preserving the previous slab.

    Removal mask is shared by +/-1 Fourier planes. Floating sums are
    estimated, not interval-certified. An exact version obeys
       ||unpruned - pruned||_infty <= deleted_L1 <= budget
       ||pruned/(1+deleted_L1)||_infty <= 1
    when the parent is the unmodified original p_d.
    """
    magnitudes=np.abs(new)
    if budget<=0:
        cutoff=0.
        keep=magnitudes>0
    else:
        ordered=np.sort(magnitudes.ravel())
        prefix=np.cumsum(ordered)
        k=int(np.searchsorted(prefix,budget/2,side="right"))
        if k>=len(ordered):
            cutoff=float("inf")
            keep=np.zeros_like(new,dtype=bool)
        else:
            cutoff=float(ordered[k])
            keep=magnitudes>=cutoff
    lost=float(2*np.sum(magnitudes[~keep]))
    # The threshold chosen via prefix makes the loss <= budget
    # within floating-point accumulation noise.
    if lost>budget*(1+1e-8)+1e-12:
        raise ArithmeticError("L1 budget exceeded")
    trimmed=float(old_energy+2*np.sum(new[keep]**2))
    normalized=trimmed/(1+lost)**2
    return {"budget":float(budget),
            "cutoff":cutoff,"new_frequency_sign_pairs_retained":int(np.count_nonzero(keep)),
            "total_support_estimate":int(old_terms+2*np.count_nonzero(keep)),
            "deleted_L1_estimate":lost,
            "unscaled_energy_estimate":trimmed,
            "norm_one_rescaled_energy_estimate":normalized,
            "warning":"FFT float64 estimates, not rigorous interval-certified numerical bounds"}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parent-depth",type=int,default=5,choices=[2,3,4,5])
    ap.add_argument("--budgets",type=float,nargs="+",default=[0.0001,0.0005,0.001,0.005])
    ap.add_argument("--max-cells",type=int,default=3_000_000)
    ap.add_argument("--output",default="")
    a=ap.parse_args()
    parent,exponent=original_poly_exact(a.parent_depth)
    squares,meta=convolve_dense(parent,exponent,a.max_cells)
    result={
      "algorithm":"original OpenAI p_d -> p_{d+1} by multidimensional FFT convolution",
      "stage":"Original Section 4 Fourier recurrence, BEFORE spreading",
      "precision":"float64, not verified interval arithmetic",
      "meta":meta,
      "pruning":[prune_new_symmetric(squares,meta["old_energy"],len(parent),t)
                 for t in a.budgets],
      "final_N":"Not computed; Sections 4 spreading, 5 and 6 have not run"
    }
    blob=json.dumps(result,indent=2)
    print(blob)
    if a.output:
        path=Path(a.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(blob+"\n")
if __name__=="__main__":
    main()
