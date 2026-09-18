# Exact beta_3 and finite-Z excess-charge bound

Research draft on the exact radial mean-field constant beta_3 arising in the Hundertmark–Pattakos–Schulz (HPS) excess-charge framework, together with its explicit finite-N / finite-Z propagation.

## Structure

- main.tex — master LaTeX file and abstract
- chapter/01_introduction.tex — context and the two main theorems
- chapter/02_radial_reduction.tex — Newton/radial reduction
- chapter/03_logarithmic_problem.tex — log-radius transform, compactness, obstacle condition
- chapter/04_extremizer.tex — explicit cosine extremizer and obstacle factorization
- chapter/05_global_optimality.tex — atomlessness, endpoint identities, conditional negativity, global proof
- chapter/06_atomic_consequence.tex — HPS finite-N propagation and the sharpened finite-Z theorem
- chapter/07_outlook.tex — route toward Pauli/Hardy stability improvements beyond the cubic mean-field obstruction
- chapter/A_proof_audit.tex — audit of the delicate radial and finite-N steps
- references.bib — bibliography
- verify_constants.py — numerical audit of the closed finite-Z constants

## Exact mean-field constant

    beta_3 = 3 exp(2 a_*) / (1 + 3 exp(2 a_*))
    a_*    = arctan(sqrt(2)) / sqrt(2)

    beta_3      = 0.9205346822904167...
    beta_3^{-1} = 1.086325175181735...

## Finite-Z theorem

Using the HPS cubic-weight finite-configuration estimate and weighted kinetic bound with the exact beta_3, the draft proves, for every Z >= 4,

    N_c(Z) <
        beta_3^{-1} Z
      + 3.812 Z^(1/3)
      + 0.01294
      + 0.18185 Z^(-1/3)
      + 0.01941 Z^(-2/3).

In particular,

    N_c(Z) < 1.086326 Z + 3.90 Z^(1/3),   Z >= 4.

The lower-order coefficients before safe rounding are

    a1* = 3.811470425068532...
    a2* = 0.012932442561687...
    a3* = 0.181840841982563...
    a4* = 0.019402438389967...

The a1 audit uses the sharper consequence of Lieb's bound valid for Z >= 4,

    N/Z < 2 + 1/Z <= 9/4,

and checks both endpoints of the HPS optimization function. This avoids an artificial endpoint loss that appears if the interval is unnecessarily enlarged to 5/2.

## Compile

    pdflatex main.tex
    bibtex main
    pdflatex main.tex
    pdflatex main.tex

Working branch: fix/finite-N-excess-charge
