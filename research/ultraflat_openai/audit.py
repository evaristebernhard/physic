#!/usr/bin/env python3
"""Faithful Section-4 parameter audit for OpenAI, Ultraflat real Littlewood polynomials.

Implements the EXACT original p_j Fourier recursion, signed frequency pairs,
the original one-direction-per-pair W_j construction, exact dyadic tau search,
and a transparent (non-certified) chirp T/support proxy.

It does NOT construct the Section-5 waves, invoke hypergraph packing, or find N0.
"""
import argparse
import json
import math
import time
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def next_poly(p, cap=12000):
    """Fourier coefficients as exact rational numbers; p_0 is {}."""
    if not p:
        return {(-1,): Q(1,4), (1,): Q(1,4)}
    d=len(next(iter(p)))
    sq=defaultdict(Q)
    for s,b in p.items():
        for t,c in p.items():
            sq[tuple(x+y for x,y in zip(s,t))]+=b*c
        if len(sq)>cap:
            raise RuntimeError("MAX_TERMS_SQUARED")
    sq[(0,)*d]-=1
    out={s+(0,):b for s,b in p.items()}
    for s,b in sq.items():
        if b:
            out[s+(-1,)]=-b/4
            out[s+(1,)]=-b/4
    out={s:b for s,b in out.items() if b}
    if len(out)>cap:
        raise RuntimeError("MAX_TERMS")
    return out

def energy(p):
    return sum((b*b for b in p.values()),Q())

def reps_of(p):
    return sorted((s,b) for s,b in p.items()
                  if s>tuple(-x for x in s))

def gamma_of(reps):
    d=len(reps[0][0])
    K=max(abs(x) for s,b in reps for x in s)
    B=2*K+1
    g=tuple(B**i for i in range(d))
    assert all(dot(s,g)!=0 for s,b in reps)
    return g

def perpendiculars(reps, limit=30):
    freq=[s for s,b in reps]
    d=len(freq[0])
    if d==2:
        return [(Q(-s[1]),Q(s[0])) for s in freq]
    ans=[]
    for i,s in enumerate(freq):
        candidate=None
        for radius in range(1,limit+1):
            for vals in product(range(-radius,radius+1),repeat=d):
                if max(abs(x) for x in vals)!=radius:
                    continue
                v=tuple(Q(int(x)) for x in vals)
                proj=tuple(v[k]-Q(dot(s,v),dot(s,s))*s[k]
                           for k in range(d))
                if any(proj) and all(dot(other,proj)!=0
                                      for j,other in enumerate(freq) if j!=i):
                    candidate=proj
                    break
            if candidate is not None:
                break
        if candidate is None:
            raise RuntimeError("W0_SEARCH_BUDGET")
        ans.append(candidate)
    return ans

def w1_of(reps,gamma,w0):
    ans=[]
    for j,(s,b) in enumerate(reps):
        prod=Q(1)
        for k,w in enumerate(w0):
            if k!=j:
                prod*=abs(dot(s,w))
        req=b*b/(abs(dot(s,gamma))*prod)
        ans.append(tuple(req*Q(x,dot(s,s)) for x in s))
    return ans

def choose_tau(reps,gamma,w0,w1,delta,max_bits=1024):
    for bits in range(max_bits+1):
        tau=Q(1,2**bits)
        worst=Q()
        passed=True
        for s,b in reps:
            lhs=Q(1)
            for v0,v1 in zip(w0,w1):
                f=dot(s,v0)+tau*dot(s,v1)
                if f==0:
                    passed=False
                    break
                lhs*=abs(f)
            if not passed:
                break
            target=tau*b*b/abs(dot(s,gamma))
            e=abs(lhs/target-1)
            worst=max(worst,e)
            if e>delta:
                passed=False
                break
        if passed:
            return tau,bits,worst
    raise RuntimeError("TAU_BUDGET")

def analytic_obstruction(delta):
    """Necessary d; lower bound on |A| for ORIGINAL untrimmed p_d boxes.

    The |A| theorem assumes the original equal-magnitude boxes and
    coefficient/weight ratios 1<=|c_a|/sqrt(w_a)<=2, not arbitrary
    modifications to Proposition 4.1.  It is NOT a bound on final N.
    """
    d=max(2,math.ceil(math.log(delta)/math.log(7/8)-1e-12))
    exponent=(1 << (d+2))+d-11
    return {"delta":delta,"necessary_d":d,
            "A_lower_bound_log10":exponent*math.log10(2),
            "warning":"A-support lower bound, NOT final N0"}

def run_depth(d,delta,max_terms,max_D):
    t0=time.monotonic()
    p={}
    for _ in range(d):
        p=next_poly(p,max_terms)
    reps=reps_of(p)
    result={"d":d,"support_p":len(p),"D":len(reps),
            "energy":float(energy(p)),
            "energy_gate":float(1-energy(p))<=delta}
    if d<2:
        result["status"]="PROP4_REQUIRES_D>=2"
        return result
    if len(reps)>max_D:
        result["status"]="DIRECTION_BUDGET"
        result["elapsed_seconds"]=round(time.monotonic()-t0,4)
        return result
    gamma=gamma_of(reps)
    w0=perpendiculars(reps)
    w1=w1_of(reps,gamma,w0)
    tau,bits,err=choose_tau(reps,gamma,w0,w1,Q(str(delta)))
    all_beta=[abs(dot(s,x)+tau*dot(s,y))
              for s,b in reps for x,y in zip(w0,w1)]
    beta_min=min(all_beta)
    B1=float(sum(abs(b) for b in p.values()))
    alpha=min(1/len(reps),delta/(math.e*len(reps)*B1))
    # PROXY only. Real Lemma 4.2 also requires T depending on smooth g.
    # The formula is not a necessary OR a sufficient T bound.
    scale=max(4*len(reps)/delta,alpha**-2)
    logT=math.log10(scale)-math.log10(float(beta_min))
    logs=[]
    for s,b in reps:
        logbox=len(reps)*logT+sum(
            math.log10(float(abs(dot(s,x)+tau*dot(s,y))))
            for x,y in zip(w0,w1))
        logs.append(logbox)
    mx=max(logs)
    logA=mx+math.log10(2*sum(10**(v-mx) for v in logs))
    result.update({"status":"SECTION4_PROXY","gamma":list(gamma),
                   "tau_bits":bits,"tau_max_relative_error":float(err),
                   "alpha":alpha,"log10_beta_min":math.log10(float(beta_min)),
                   "proxy_log10_T":logT,"proxy_log10_A":logA,
                   "proxy_warning":"NOT certified T or A; ignores g-dependent Lemma 4.2 constants",
                   "elapsed_seconds":round(time.monotonic()-t0,4)})
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depths",type=int,nargs="+",default=[2,3,4,5])
    p.add_argument("--delta",type=float,default=.005)
    p.add_argument("--max-terms",type=int,default=12000)
    p.add_argument("--max-directions",type=int,default=64)
    p.add_argument("--output",default="")
    a=p.parse_args()
    if not 0<a.delta<.01:
        raise SystemExit("delta must be in (0, 1/100), as in Proposition 4.1")
    results=[]
    for depth in a.depths:
        try:
            results.append(run_depth(depth,a.delta,a.max_terms,a.max_directions))
        except (RuntimeError,ValueError,OverflowError) as exc:
            results.append({"d":depth,"status":"STOPPED","reason":str(exc)})
    data={"paper":"OpenAI, Ultraflat real Littlewood polynomials (2026-10-05)",
          "method":"exact original p_j and original D-direction spreading",
          "scope":"Section 4 only: full N0 and stages 5-6 NOT executed",
          "analytic_obstruction":analytic_obstruction(a.delta),
          "runs":results}
    payload=json.dumps(data,ensure_ascii=False,indent=2)
    print(payload)
    if a.output:
        from pathlib import Path
        q=Path(a.output)
        q.parent.mkdir(parents=True,exist_ok=True)
        q.write_text(payload+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
