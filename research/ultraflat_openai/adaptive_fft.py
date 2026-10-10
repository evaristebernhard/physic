#!/usr/bin/env python3
"""Adaptive FFT + Fourier L1 pruning for OpenAI's ORIGINAL p_d recursion.

Start with exactly represented depth-5 coefficients from exact_pruning.py;
perform p -> p+(1-p^2)cos(2*pi*y_new)/2 by a multidimensional FFT; prune
only new +/-1-frequency slabs with an explicit L1 loss. Iterate until
a conservative dense-cell guard is reached.

Error propagation, IN EXACT ARITHMETIC:
  E_{j+1} <= (2+E_j/2)*E_j + eta_{j+1}.
Proof: T_y(u)-T_y(v)=(u-v)[1-(u+v)cos(2*pi*y)/2].
Original |p_j|<=1 and approximate |q_j|<=1+E_j.
Deleting Fourier coefficients with L1 mass eta adds at most eta in L-inf.

Numerical FFTs are performed in float64. Thus the reported error
recurrence only bounds analytical pruning, NOT roundoff. Certified
interval FFT and tracking roundoff are still needed for a fully rigorous
Proposition 4.1 input. No coefficient spreading or final N performed.
"""
import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

try:
    from .exact_pruning import step
except ImportError:
    from exact_pruning import step


def initial_exact_pruned(depth=5, eta=Fraction(5,1000)):
    p={}
    exponent=0
    for j in range(depth):
        p,exponent=step(p,exponent,max_pairs=2_000_000)
    R=len(p)
    cutoff=eta*(1<<exponent)/R
    kept={s:b for s,b in p.items() if abs(b)>=cutoff}
    discarded=sum(abs(b) for s,b in p.items() if abs(b)<cutoff)
    loss=Fraction(discarded,1<<exponent)
    assert loss<=eta
    p_float={s:float(b)/float(1<<exponent) for s,b in kept.items()}
    return p_float, float(loss), {
        "depth":depth,"raw_terms":R,"retained_terms":len(kept),
        "deleted_L1_exact":str(loss),"deleted_L1_float":float(loss),
        "energy_before_pruning":float(Fraction(sum(b*b for b in p.values()),1<<(2*exponent))),
        "energy_after_pruning_float":sum(v*v for v in p_float.values())
    }


def next_step_fft(p, budget, max_cells=35_000_000, max_modes=2_000_000):
    """One approximate original recurrence; stops before resource overflow."""
    t0=time.perf_counter()
    coordinates=np.array(list(p),dtype=int)
    lo=coordinates.min(axis=0)
    hi=coordinates.max(axis=0)
    original_shape=tuple(int(x) for x in (hi-lo+1))
    convolution_shape=tuple(2*x-1 for x in original_shape)
    cells=math.prod(convolution_shape)
    if cells>max_cells:
        return None,{"status":"FFT_CELL_BUDGET","dense_cells_required":cells,
                     "dense_shape":list(convolution_shape),
                     "max_cells":max_cells}
    arr=np.zeros(original_shape,dtype=np.float64)
    for s,v in p.items():
        arr[tuple(np.array(s)-lo)]=v
    fft_start=time.perf_counter()
    c=fftconvolve(arr,arr,mode="full")
    c=(c+np.flip(c))/2
    fft_seconds=time.perf_counter()-fft_start
    c[tuple(-2*lo)]-=1.0
    slab=-c/4.0
    magnitudes=np.abs(slab)
    if budget>0:
        ordered=np.sort(magnitudes.ravel())
        prefix=np.cumsum(ordered)
        k=int(np.searchsorted(prefix,budget/2,side="right"))
        threshold=float(ordered[k]) if k<len(ordered) else float("inf")
        keep=magnitudes>=threshold
    else:
        threshold=0.
        keep=magnitudes>0
    lost=2*float(np.sum(magnitudes[~keep]))
    if lost>budget*(1+1e-8)+1e-12:
        raise ArithmeticError("Pruning exceeded desired numerical L1 budget")
    indices=np.array(np.where(keep)).T
    if len(p)+2*len(indices)>max_modes:
        return None,{"status":"OUTPUT_MODE_BUDGET",
                     "proposed_output_terms":len(p)+2*len(indices),
                     "max_modes":max_modes,"dense_cells":cells}
    out={s+(0,):v for s,v in p.items()}
    start=2*lo
    for idx in indices:
        s=tuple(int(v) for v in idx+start)
        value=float(slab[tuple(idx)])
        out[s+(-1,)]=value
        out[s+(1,)]=value
    old_energy=sum(v*v for v in p.values())
    estimated_energy=old_energy+2*float(np.sum((slab[keep])**2))
    return out,{
        "status":"COMPUTED_FLOAT64",
        "parent_terms":len(p),
        "parent_shape":list(original_shape),
        "dense_shape":list(convolution_shape),
        "dense_cells":cells,
        "new_terms":2*len(indices),
        "output_terms":len(out),
        "dropped_L1_float":lost,
        "threshold_float":threshold,
        "energy_float":estimated_energy,
        "fft_seconds":fft_seconds,
        "stage_seconds":time.perf_counter()-t0}


def run(target_depth=9, eta_start="0.005", stage_budget=0.0005,
        max_cells=35_000_000, max_modes=2_000_000):
    if target_depth<5:
        raise ValueError("target depth >=5 required")
    p,err,first=initial_exact_pruned(5,Fraction(eta_start))
    stage=[dict(first,analytic_pruning_error_upper=err)]
    for d in range(6,target_depth+1):
        nxt,record=next_step_fft(p,stage_budget,max_cells,max_modes)
        record["depth"]=d
        if nxt is None:
            stage.append(record)
            break
        err=(2+err/2)*err+record["dropped_L1_float"]
        record["analytic_pruning_error_upper"]=err
        record["energy_interval_conditional_on_exact_fft"]=[
            max(0,math.sqrt(record["energy_float"])-err)**2,
            (math.sqrt(record["energy_float"])+err)**2]
        stage.append(record)
        p=nxt
    return {
      "method":"OpenAI unmodified p_d + adaptive multidimensional FFT and Fourier L1 pruning",
      "eta_start":eta_start,"eta_each_following_stage":stage_budget,
      "pruning_error_bound_rule":"E_next <= (2+E/2)E + new_deleted_L1, assuming exact convolution",
      "numeric_status":"Float64 FFT used, so values and error estimates are not interval-certified",
      "no_claim":"No Proposition 4.1 balanced F, no Sections 5-6, no final N or N0",
      "stages":stage}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target-depth",type=int,default=9)
    ap.add_argument("--eta-start",default="0.005")
    ap.add_argument("--stage-budget",type=float,default=0.0005)
    ap.add_argument("--max-cells",type=int,default=35_000_000)
    ap.add_argument("--max-modes",type=int,default=2_000_000)
    ap.add_argument("--output",default="")
    a=ap.parse_args()
    result=run(a.target_depth,a.eta_start,a.stage_budget,
               a.max_cells,a.max_modes)
    txt=json.dumps(result,indent=2)
    print(txt)
    if a.output:
        path=Path(a.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(txt+"\n")
if __name__=="__main__":
    main()
