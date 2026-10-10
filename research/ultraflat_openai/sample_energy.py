#!/usr/bin/env python3
"""Monte Carlo evaluation of the UNMODIFIED OpenAI recurrence p_d.

This computes M_d=E[p_d(Y_1,...,Y_d)^2] with independent uniform torus
samples without constructing its astronomically large Fourier support.
It makes no claim about Section-4 spreading or the final Littlewood
polynomial length N.

With n samples and a preselected max depth D, Hoeffding + a union bound
gives with probability at least 1-confidence_failure:
  |estimated (1-M_d) - true (1-M_d)| <=
  sqrt(log(2 D / confidence_failure)/(2 n))
simultaneously for all 1<=d<=D.
"""
import argparse
import json
import math
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples",type=int,default=1_000_000)
    ap.add_argument("--max-depth",type=int,default=90)
    ap.add_argument("--seed",type=int,default=20261010)
    ap.add_argument("--delta",type=float,default=.005)
    ap.add_argument("--failure-probability",type=float,default=.01)
    ap.add_argument("--output",default="")
    args=ap.parse_args()
    if args.samples<2 or args.max_depth<1:
        ap.error("samples>=2 and max-depth>=1 required")
    if not 0<args.delta<1 or not 0<args.failure_probability<1:
        ap.error("delta and failure-probability must lie in (0,1)")
    try:
        import numpy as np
    except ImportError as ex:
        raise SystemExit("Install numpy to run this separate Monte Carlo experiment") from ex
    rng=np.random.default_rng(args.seed)
    x=np.zeros(args.samples,dtype=np.float64)
    radius=math.sqrt(
        math.log(2*args.max_depth/args.failure_probability)/
        (2*args.samples))
    selected=set([2,3,4,5,10,20,30,40,45,50,55,56,57,58,59,60,
                  61,62,63,64,65,66,67,70,75,80,90,args.max_depth])
    first_estimate=None
    first_high_confidence=None
    checkpoints=[]
    for d in range(1,args.max_depth+1):
        x+=.5*(1-x*x)*np.cos(2*math.pi*rng.random(args.samples))
        residual=1-x*x
        observed=float(np.mean(residual))
        if first_estimate is None and observed<=args.delta:
            first_estimate=d
        if first_high_confidence is None and observed+radius<=args.delta:
            first_high_confidence=d
        if d in selected:
            checkpoints.append({
                "d":d,
                "estimated_deficit":observed,
                "standard_error_descriptive":float(np.std(residual,ddof=1)/math.sqrt(args.samples)),
                "hoeffding_upper_simultaneous":min(1.,observed+radius),
                "energy_requirement_supported_at_confidence":observed+radius<=args.delta
            })
    result={
      "scope":"Monte Carlo of original p_d energy ONLY (not full Proposition 4.1 or final N)",
      "samples":args.samples,"depth":args.max_depth,"seed":args.seed,
      "delta":args.delta,"failure_probability":args.failure_probability,
      "uniform_hoeffding_radius":radius,
      "estimated_first_depth":first_estimate,
      "first_depth_supported_with_simultaneous_confidence":first_high_confidence,
      "important":"Probabilistic sampling guarantee, NOT deterministic mathematical proof",
      "checkpoints":checkpoints}
    blob=json.dumps(result,indent=2,ensure_ascii=False)
    print(blob)
    if args.output:
        dest=Path(args.output)
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(blob+"\n",encoding="utf8")
if __name__=="__main__":
    main()
