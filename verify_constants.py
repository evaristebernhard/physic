"""Reproduce the explicit constants in the finite-Z theorem.

This script is only a numerical audit of the closed formulas in the paper;
the proof is contained in chapter/06_atomic_consequence.tex.
"""

from math import atan, exp, pi, sqrt

a_star = atan(sqrt(2.0)) / sqrt(2.0)
beta3 = 3.0 * exp(2.0 * a_star) / (1.0 + 3.0 * exp(2.0 * a_star))
beta_inv = 1.0 / beta3

c0 = (
    sqrt(5.0)
    * (2.0 * 1.456 / (9.0 * pi**2)) ** (1.0 / 3.0)
    * (4.0 * pi ** (2.0 / 3.0) / sqrt(15.0))
)

a2 = beta_inv / 84.0
a3 = (c0 / 5.0) * (5.0 / 12.0) ** (2.0 / 3.0) * beta3 ** (-1.0 / 3.0)
a4 = (c0 / 84.0) * beta3 ** (-1.0 / 3.0)

A = 3.0 * (3.0 / 10.0) ** (1.0 / 3.0) * beta3 ** (-2.0 / 3.0)
B = c0 * beta_inv


def F(x: float) -> float:
    return A * x ** (1.0 / 3.0) + B * x ** (-2.0 / 3.0)


left = beta_inv
right = 9.0 / 4.0
a1 = max(F(left), F(right))

absorption_at_Z4 = (
    0.01294 * 4.0 ** (-1.0 / 3.0)
    + 0.18185 * 4.0 ** (-2.0 / 3.0)
    + 0.01941 / 4.0
)

print(f"a_*            = {a_star:.15f}")
print(f"beta_3         = {beta3:.16f}")
print(f"beta_3^(-1)    = {beta_inv:.16f}")
print(f"c_0             = {c0:.16f}")
print(f"F(beta^-1)      = {F(left):.16f}")
print(f"F(9/4)          = {F(right):.16f}")
print(f"a_1*            = {a1:.16f}")
print(f"a_2*            = {a2:.16f}")
print(f"a_3*            = {a3:.16f}")
print(f"a_4*            = {a4:.16f}")
print(f"absorb at Z=4   = {absorption_at_Z4:.16f}")

assert beta_inv < 1.086326
assert a1 < 3.812
assert a2 < 0.01294
assert a3 < 0.18185
assert a4 < 0.01941
assert 3.812 + absorption_at_Z4 < 3.90
