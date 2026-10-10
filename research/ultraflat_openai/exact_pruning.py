#!/usr/bin/env python3
"""Exact Fourier audit and norm-preserving spectral pruning of OpenAI p_d.

No nonlinear modification to the paper's recursive polynomial is made.
All Fourier coefficients are represented by integers with a common power
of two denominator; no floats are used for recurrence or pruning.
This is an input optimization BEFORE the original Proposition 4.1 spreading.
It does not construct a Littlewood polynomial or establish final N.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

def step(poly: dict[tuple[int,...],int], exponent: int,
         max_pairs: int = 2_000_000):
    """p_d coefficients are numerators / 2**exponent.

    From exponent q, next exponent is 2q+2.  The old p_{d-1}
    coefficients persist with zero last coordinate; the new +/-1
    last-coordinate coefficients equal 1/4 times (1-p^2).
    """
    if not poly:
        return {(-1,):1,(1,):1},2
    if len(poly)**2 > max_pairs:
        raise RuntimeError(
            f"CONVOLUTION_BUDGET: {len(poly)}**2 > {max_pairs}")
    q2=2*exponent+2
    convolution=defaultdict(int)
    for s,b in poly.items():
        for t,c in poly.items():
            convolution[tuple(x+y for x,y in zip(s,t))]+=b*c
    convolution[(0,)*len(next(iter(poly)))] -= 1<<(2*exponent)
    scale=1<<(q2-exponent)
    nxt={s+(0,):b*scale for s,b in poly.items()}
    for s,b in convolution.items():
        if b:
            nxt[s+(-1,)]=-b
            nxt[s+(1,)]=-b
    return {s:b for s,b in nxt.items() if b},q2

def prune(poly, exponent, eta):
    """Delete |b_s|<eta/|S|; renormalize by 1+actual L1 loss.

    eta is a Fraction.  Guaranteed ||q||infty<=1 if ||p||infty<=1.
    This exact norm control makes the original paper's section-4
    supnorm argument reusable without an extra triangle-error budget.
    """
    R=len(poly)
    bound=eta/R
    cutoff=bound*(1<<exponent)
    kept={s:b for s,b in poly.items() if abs(b)>=cutoff}
    deleted={s:b for s,b in poly.items() if abs(b)<cutoff}
    lost_l1=Fraction(sum(abs(b) for b in deleted.values()),1<<exponent)
    lost_l2=Fraction(sum(b*b for b in deleted.values()),1<<(2*exponent))
    normalization=1+lost_l1
    assert lost_l1<=eta
    assert lost_l2<=eta*eta/R
    assert all(tuple(-x for x in s) in kept for s in kept)
    assert all(tuple(-x for x in s) in deleted for s in deleted)
    raw_energy=Fraction(sum(b*b for b in poly.values()),1<<(2*exponent))
    remaining_energy=(raw_energy-lost_l2)/normalization**2
    return kept,lost_l1,lost_l2,normalization,raw_energy,remaining_energy

def projection_gamma(keys):
    d=len(next(iter(keys)))
    K=max(abs(x) for s in keys for x in s)
    base=2*K+1
    return tuple(base**i for i in range(d))

def log10_required_support_lower_bound(keys, exponent, gamma):
    """For a fixed gamma and original uniform-per-box construction.

    From 1<=|c_a|/sqrt(w_a)<=2 and n_s>=1,
    |A| >= max_s x_s /(4 min_s x_s), x_s=|b_s|^2/|s.gamma|.
    Only a lower bound for the given gamma, NOT best-possible gamma or N.
    """
    values=[]
    for s,b in keys.items():
        rho=abs(sum(x*g for x,g in zip(s,gamma)))
        assert rho>0
        values.append(2*math.log10(abs(b))-math.log10(rho))
    return max(0.,max(values)-min(values)-math.log10(4))

def audit(depths, eta=Fraction(5,1000), max_pairs=2_000_000):
    poly={}
    q=0
    results=[]
    levels=set(depths)
    for d in range(1,max(levels)+1):
        try:
            poly,q=step(poly,q,max_pairs)
        except RuntimeError as e:
            for j in sorted(levels):
                if j>=d: results.append({"d":j,"status":"BUDGET","reason":str(e)})
            break
        if d not in levels: continue
        kept,l1,l2,norm,raw_en,clean_en=prune(poly,q,eta)
        gamma=projection_gamma(poly)
        pairs=len(poly)//2
        kept_pairs=len(kept)//2
        info={
            "d":d,
            "status":"COMPLETE_ORIGINAL_RECURRENCE_AND_PRUNING",
            "dyadic_exponent":q,
            "raw_support":len(poly),
            "raw_sign_pairs":pairs,
            "raw_energy":float(raw_en),
            "pruned_support":len(kept),
            "pruned_sign_pairs":kept_pairs,
            "deleted_L1":float(l1),
            "deleted_L2_squared":float(l2),
            "normalization":float(norm),
            "norm_one_energy_after_pruning":float(clean_en),
            "raw_min_coefficient":float(Fraction(min(abs(b) for b in poly.values()),1<<q)),
            "pruned_min_coefficient":float(Fraction(min(abs(b) for b in kept.values()),1<<q)/norm),
            "original_equal_box_log10_support_lower_bound_fixed_gamma":
                log10_required_support_lower_bound(poly,q,gamma),
            "pruned_equal_box_log10_support_lower_bound_fixed_gamma":
                log10_required_support_lower_bound(kept,q,gamma),
            "interpretation":"Lower bounds depend on fixed gamma and equal-box scheme; NOT final N"
        }
        results.append(info)
    return {
      "status":"section4_original_p_recursion_with_optional_pre_spreading_pruning",
      "eta":float(eta),
      "guarantee":"pruned polynomial scaled by 1+L1 loss has supremum <=1; energy is exact",
      "not_done":"Lemma 4.2 full chirp rounding, Section 5 packing, Section 6 final N",
      "runs":results
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depths",type=int,nargs="+",default=[2,3,4,5,6])
    parser.add_argument("--eta",default="0.005")
    parser.add_argument("--max-pairs",type=int,default=2_000_000)
    parser.add_argument("--output",default="")
    a=parser.parse_args()
    if any(d<1 for d in a.depths): parser.error("depths must be positive")
    eta=Fraction(a.eta)
    if not 0<eta<1: parser.error("eta must be in (0,1)")
    result=audit(a.depths,eta,a.max_pairs)
    blob=json.dumps(result,indent=2,ensure_ascii=False)
    print(blob)
    if a.output:
        path=Path(a.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(blob+"\n",encoding="utf8")

if __name__=="__main__":
    main()
