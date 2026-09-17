# Exact beta_3 paper draft

Research draft on the exact radial mean-field constant `beta_3` arising in the Hundertmark–Pattakos–Schulz excess-charge framework.

## Structure

- `main.tex` — master LaTeX file
- `chapter/01_introduction.tex` — context and main theorem
- `chapter/02_radial_reduction.tex` — Newton/radial reduction
- `chapter/03_logarithmic_problem.tex` — log-radius transform, compactness, obstacle condition
- `chapter/04_extremizer.tex` — explicit cosine extremizer and obstacle factorization
- `chapter/05_global_optimality.tex` — atomlessness, endpoint identities, conditional negativity, global proof
- `chapter/06_atomic_consequence.tex` — consequence for the HPS mean-field input
- `chapter/07_outlook.tex` — route toward Pauli/Hardy stability improvements
- `chapter/A_proof_audit.tex` — dedicated audit of delicate proof steps
- `references.bib` — bibliography

## Main claimed value

```text
beta_3 = 3 exp(2 a_*) / (1 + 3 exp(2 a_*))
a_*    = arctan(sqrt(2)) / sqrt(2)

beta_3      = 0.9205346822904167...
beta_3^{-1} = 1.086325175181735...
```

The draft deliberately distinguishes this exact mean-field result from a fully propagated finite-Z atomic excess-charge theorem. The finite-N HPS comparison and lower-order optimization still need to be rerun before claiming a final replacement of the published `1.1185 Z + O(Z^{1/3})` bound.

## Compile

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Branch: `paper/exact-beta3`
