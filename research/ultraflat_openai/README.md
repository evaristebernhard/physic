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
