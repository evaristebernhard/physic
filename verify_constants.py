"""Interval audit for the explicit constants in the finite-Z theorem.

Requires mpmath. All assertions are interval comparisons: they only pass if
the whole computed interval lies on the safe side of the rounded constant
used in the paper.
"""

from mpmath import iv

ONE = iv.mpf(1)
TWO = iv.mpf(2)
THREE = iv.mpf(3)
FIVE = iv.mpf(5)

sqrt2 = iv.sqrt(TWO)
a_star = iv.atan2(sqrt2, ONE) / sqrt2
beta3 = THREE * iv.exp(TWO * a_star) / (ONE + THREE * iv.exp(TWO * a_star))
beta_inv = ONE / beta3

c0 = (
    iv.sqrt(FIVE)
    * (TWO * iv.mpf("1.456") / (iv.mpf(9) * iv.pi**2)) ** (ONE / THREE)
    * (iv.mpf(4) * iv.pi ** (TWO / THREE) / iv.sqrt(iv.mpf(15)))
)

a2 = beta_inv / iv.mpf(84)
a3 = (
    (c0 / FIVE)
    * (FIVE / iv.mpf(12)) ** (TWO / THREE)
    * beta3 ** (-ONE / THREE)
)
a4 = (c0 / iv.mpf(84)) * beta3 ** (-ONE / THREE)

A = (
    THREE
    * (THREE / iv.mpf(10)) ** (ONE / THREE)
    * beta3 ** (-TWO / THREE)
)
B = c0 * beta_inv


def F(x):
    return A * x ** (ONE / THREE) + B * x ** (-TWO / THREE)


left = beta_inv
right = iv.mpf(9) / iv.mpf(4)
F_left = F(left)
F_right = F(right)

absorption_at_Z4 = (
    iv.mpf("0.01294") * iv.mpf(4) ** (-ONE / THREE)
    + iv.mpf("0.18185") * iv.mpf(4) ** (-TWO / THREE)
    + iv.mpf("0.01941") / iv.mpf(4)
)

print("a_*           =", a_star)
print("beta_3        =", beta3)
print("beta_3^(-1)   =", beta_inv)
print("c_0           =", c0)
print("F(beta^-1)    =", F_left)
print("F(9/4)        =", F_right)
print("a_2*          =", a2)
print("a_3*          =", a3)
print("a_4*          =", a4)
print("absorb at Z=4 =", absorption_at_Z4)

assert beta_inv < iv.mpf("1.086326")
assert F_left > F_right
assert F_left < iv.mpf("3.812")
assert a2 < iv.mpf("0.01294")
assert a3 < iv.mpf("0.18185")
assert a4 < iv.mpf("0.01941")
assert iv.mpf("3.812") + absorption_at_Z4 < iv.mpf("3.90")

print("All interval-certified margins passed.")
