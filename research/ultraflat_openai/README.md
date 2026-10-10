# OpenAI ultraflat real Littlewood polynomials — original-method reproduction

Isolated research module. The rest of the \`physic\` repository (atomic excess-charge paper) is unchanged.

**Source:** OpenAI, *Ultraflat real Littlewood polynomials*, October 5, 2026, Sections 3–6: https://github.com/openai/math/tree/main/preprints/Ultraflat-real-Littlewood-polynomials-October-5-2026

**Scope and scientific status.** This is the FIRST executable reproduction experiment, not a claim to reproduce the complete theorem. The Python code follows the **unchanged** Section 4 recurrence and constructs the same one-direction-per-sign-pair coefficients \(\gamma,W_j^{(0)},W_j^{(1)},\tau\). The rational Fourier coefficients and the relative errors of the Section 4 equation (4.5) are checked exactly. It additionally estimates the scale of \(T\) and the expanded support using clearly **non-certified** heuristic thresholds.

**Not implemented:** numerical Lemma 4.2's chirp/discrepancy rounding, Proposition 5.1's signed interval packing, phase gluing and integration, and Section 6's sign rounding. **We have not measured any successful final sequence length \(N\), and no finite \(N_0(\varepsilon)\) has been established.** The audit intentionally stops rather than inventing an \(N\).

## Reproduce

No third-party Python dependencies:

\`\`\`bash
python3 -m unittest discover -s research/ultraflat_openai -v
python3 research/ultraflat_openai/audit.py \
  --depths 2 3 4 5 --delta 0.005 \
  --max-directions 64 --max-terms 12000 \
  --output research/ultraflat_openai/results/delta_0p005.json
\`\`\`

\`--max-directions\` and \`--max-terms\` are deliberate resource guards, not theoretical limitations. A stopped run is **not** a result about possible final \(N\).

## Baseline measurements (Python 3, local execution)

At the admissible internal parameter \(\delta=0.005 < 0.01\):

| Original recursion depth | Support of \(p_d\) | \(D\), sign pairs | \(M_d\) | Outcome |
|---|---:|---:|---:|---|
| 2 | 8 | 4 | 0.2216796875 | Exact \(W,\tau\) checked; \(\log_{10} T\) proxy 13.54, \(\log_{10}|A|\) proxy 51.82 |
| 3 | 54 | 27 | 0.3037941754 | Exact \(W,\tau\) checked; \(\log_{10} T\) proxy 42.30, \(\log_{10}|A|\) proxy 1132.74 |
| 4 | 700 | 350 | 0.3740568532 | Stop: exceeds \(D=64\) resource guard |
| 5 | — | — | — | Stop: Fourier expansion exceeds 12,000-term guard |

For depths 2 and 3 the energy condition \(\|p_d\|_2^2\ge1-\delta=0.995\) **fails**, so their enormous expansion sizes are **toy-stage scale diagnostics**, not valid full instances of Proposition 4.1. The numbers for \(T\) and \(|A|\) are proxy estimates based on a box-count threshold and unit-constant rounding scale, **not** rigorous necessary or sufficient bounds for Lemma 4.2.

## A rigorous obstruction for the unchanged Section 4 route

For the original recursion,
\[
e_d:=1-\|p_d\|_2^2\ge(7/8)^d.
\]
Indeed
\[
M_{j+1}-M_j=\tfrac18\int(1-p_j^2)^2\le\tfrac18(1-M_j).
\]
Consequently, meeting the required \(M_d\ge1-\delta\) forces
\[
d\ge\frac{\log(1/\delta)}{\log(8/7)}.
\]

There is a persistent frequency \(e_1=(1,0,\ldots,0)\) with coefficient magnitude \(1/4\). There are also \(2^{d-1}\) frequencies
\[
s=(2^{d-1},\pm2^{d-2},\ldots,\pm1)
\]
of identical coefficient magnitude \(2^{-(2^{d+1}-2)}\). For **any** vector \(\gamma\) with \(\gamma_1\ne0\), at least one of these extreme frequencies has
\[
|s\cdot\gamma|\ge2^{d-1}|\gamma_1|.
\]
Put \(x_s=|b_s|^2/|s\cdot\gamma|\). This implies
\[
\frac{x_{e_1}}{x_{s_*}}\ge 2^{2^{d+2}+d-9}.
\]
For OpenAI's original equal-magnitude boxes, the coefficient-weight ratio \(1\le |c_a|/\sqrt{w_a}\le2\) implies (with common scale \(Q\))
\[
Qx_s/4\le n_s\le Qx_s,\qquad n_s\ge1.
\]
Hence **for this specific original construction**,
\[
\boxed{|A|\ge n_{e_1}\ge2^{2^{d+2}+d-11}.}
\]
At \(\delta=0.005\), the necessary depth is \(d\ge40\), and the lower bound on the intermediate support is
\[
\log_{10}|A|\ge 1.3239\times10^{12}.
\]

This is a **mathematical lower bound on the number of intermediate Fourier modes**, **NOT** a lower bound on final \(N\), and **NOT** an obstruction to modified/optimized algorithms satisfying the same Proposition 4.1.

## Next faithful steps

The next implementable stage is numerical Lemma 4.2: Fourier integrals of a smooth truncated quadratic chirp, complex-phase discrepancy rounding, and uniform numerical/certified errors. Only when this passes can one build the finite Fourier polynomial \(F\). The later signed interval packing and joined oscillatory function \(B_N\) are distinct stages. Measuring \(N\) can begin only when those stages produce a valid candidate \(B_N\) and real Fourier coefficients for Section 6.

The principal research question is whether one can retain the **same Proposition 4.1 interface** while avoiding its original tensor-spreading support blowup. This module keeps the original benchmark intact for that comparison.


## Second experiment: exact dyadic expansion and norm-preserving pruning (2026-10-10)

The initial \`audit.py\` intentionally stops at a 12,000-term cap. The
separate \`exact_pruning.py\` represents all coefficients as integer
numerators with common denominator \(2^{q_d}\), where
\[
q_1=2,\qquad q_{d+1}=2q_d+2,\qquad q_d=2^{d+1}-2.
\]
This uses **the original paper's recursion exactly** but computes level 5
without arbitrary-precision rational objects for each coefficient.

Reproduce:
\`\`\`bash
python3 research/ultraflat_openai/exact_pruning.py \
  --depths 2 3 4 5 6 --eta 0.005 \
  --output research/ultraflat_openai/results/pruning_delta_0p005.json
\`\`\`

At \(\eta=0.005\), exact computations (not Monte Carlo):

| d | Original terms | Pruned terms | Actual deleted Fourier \(\ell^1\) | \(M_d\) before | Energy after norm-one rescaling |
|--:|--:|--:|--:|--:|--:|
| 3 | 54 | 46 | 0.00048828125 | 0.303794175386 | 0.303497688763 |
| 4 | 700 | 452 | 0.000233292580 | 0.374056853249 | 0.373882384318 |
| 5 | 18,362 | 5,994 | 0.000230615190 | 0.435424475962 | 0.435223714387 |

At d=6 the program deliberately stops before squaring 18,362 terms
(\(>337\) million coefficient-pair operations). This resource ceiling is
a computational guard, not a mathematical impossibility.

### Norm-one pruning theorem; compatible with Proposition 4.1

Let \(p(y)=\sum_{s\in S}b_s e(s\cdot y)\) be an original OpenAI auxiliary
polynomial, \(R=|S|\), \(\|p\|_\infty\le1\). For \(\eta>0\), keep exactly the
coefficients \(|b_s|\ge\eta/R\) and define
\[
\ell=\sum_{\text{deleted}}|b_s|\le\eta,\quad
q(y)=\frac{\sum_{\text{kept}}b_se(s\cdot y)}{1+\ell}.
\]
Then
\[
\|q\|_\infty\le1,\qquad
\|q\|_2^2
=\frac{\|p\|_2^2-\sum_{\rm deleted}|b_s|^2}{(1+\ell)^2}
\ge\frac{\|p\|_2^2-\eta^2/R}{(1+\eta)^2},
\]
and every nonzero retained coefficient has magnitude at least
\(\eta/[R(1+\eta)]\). The nonzero mean and frequency-nonparallel properties
of the original construction are unchanged, because only frequency
deletion and positive scalar rescaling occur.

**Useful corollary:** If \(0<\delta<1/100\), take an original \(p_d\) with
\(M_d\ge1-\delta/2\) and choose \(\eta=\delta/8\). The rescaled \(q\)
satisfies \(\|q\|_\infty\le1\), \(\|q\|_2^2\ge1-\delta\), and has a
controlled positive smallest coefficient. Therefore the ORIGINAL
equal-magnitude spreading proof of Proposition 4.1 can be applied to
\(q\) (with its new coefficients) without changing its logical interface.
This removes one severe coefficient-conditioning obstacle but does
**not** remove the \(D\)-fold tensor spreading complexity.

### Dependence on the weight vector gamma

For a **fixed** \(\gamma\), the original spreading boxes obey
\[
n_s\in[Qx_s/4,Qx_s],\quad x_s=|b_s|^2/|s\cdot\gamma|,\quad n_s\ge1.
\]
Thus
\[
|A|\ge\frac{\max_s x_s}{4\min_s x_s}.
\]
This is a genuine lower bound **for that fixed gamma and original
equal-box method**, not for an optimized gamma or for final \(N\).
With the explicit signed-base integer \(\gamma\) used by the code, the
\(\log_{10}\) lower bounds for \(|A|\) at \(d=5\) change from about
**41.62 (raw)** to **17.38 (pruned)**. The reduction is substantial but
even the pruned instance remains far beyond explicit final-wave
computation and fails the \(\delta=0.005\) energy gate.

## Third experiment: high-dimensional energy without Fourier expansion

The paper's \(p_d(Y_1,\dots,Y_d)\) may be sampled by independently
drawing \(U_j\) uniform on \([0,1)\) and updating
\[
x_j=x_{j-1}+\tfrac12(1-x_{j-1}^2)\cos(2\pi U_j).
\]
This directly samples \(1-M_d=\mathbb E(1-x_d^2)\). It does NOT generate
the Fourier support or the final signs.

\`\`\`bash
python3 -m pip install numpy
python3 research/ultraflat_openai/sample_energy.py \
  --samples 1000000 --max-depth 90 --seed 20261010 \
  --delta 0.005 --failure-probability 0.01 \
  --output research/ultraflat_openai/results/energy_delta_0p005.json
\`\`\`

Observed deficits: at depth 58, 0.005205290; 59, 0.004775503;
60, 0.004365253; 65, 0.002797559; 66, 0.002561365.
The first **point estimate** below 0.005 is depth 59, not an exact
certificate. Hoeffding's inequality and a union bound over all 90 depths,
with the stated 1,000,000 samples, give a simultaneous radius
\[
\sqrt{\log(2\cdot90/0.01)/(2\cdot 10^6)}
=0.002213383.
\]
Accordingly depth 66 is the first depth at which this **99%-confidence
sampling bound** establishes deficit at most 0.005 (conditional on
ideal i.i.d. uniform sampling and exact arithmetic). This is a
probabilistic statement, not deterministic proof of the exact integral.

**Main conclusion:** the original recurrence's needed depth is
plausibly around 59 rather than the much larger conservative analytic
upper bound. The actual bottleneck is explicitly materializing and
spreading its enormous Fourier support. These experiments do not
yield a value or rigorous numerical bound for final \(N_0(\varepsilon)\).
