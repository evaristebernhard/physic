"""Check exact Section-4 identities, not the unimplemented Sections 5-6."""
import unittest
from fractions import Fraction as Q
try:
    from .audit import (
        next_poly, energy, reps_of, gamma_of, perpendiculars,
        w1_of, choose_tau, analytic_obstruction
    )
except ImportError:
    from audit import (
        next_poly, energy, reps_of, gamma_of, perpendiculars,
        w1_of, choose_tau, analytic_obstruction
    )

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

class TestOpenAIOriginalRecurrence(unittest.TestCase):
    def test_first_three_levels(self):
        p1=next_poly({})
        self.assertEqual(energy(p1),Q(1,8))
        p2=next_poly(p1)
        self.assertEqual(energy(p2),Q(227,1024))
        self.assertEqual(len(p2),8)
        self.assertEqual(p2[(2,1)],Q(-1,64))
        p3=next_poly(p2)
        self.assertEqual(len(p3),54)
        self.assertEqual(p3[(4,2,1)],Q(-1,16384))
        for s,b in p3.items():
            self.assertEqual(b,p3[tuple(-x for x in s)])

    def test_spreading_input_exactness(self):
        p2=next_poly(next_poly({}))
        reps=reps_of(p2)
        gamma=gamma_of(reps)
        w0=perpendiculars(reps)
        w1=w1_of(reps,gamma,w0)
        for j,(s,b) in enumerate(reps):
            self.assertEqual(dot(s,w0[j]),0)
            self.assertNotEqual(dot(s,gamma),0)
            prod=Q(1)
            for k in range(len(reps)):
                if k!=j:
                    prod*=abs(dot(s,w0[k]))
            self.assertEqual(
                dot(s,w1[j])*prod,
                b*b/abs(dot(s,gamma)))
        tau,bits,worst=choose_tau(reps,gamma,w0,w1,Q(1,200))
        self.assertLessEqual(worst,Q(1,200))

    def test_raw_scheme_obstruction(self):
        obstruction=analytic_obstruction(.005)
        self.assertGreaterEqual(obstruction["necessary_d"],40)
        self.assertGreater(obstruction["A_lower_bound_log10"],1e12)

if __name__=="__main__":
    unittest.main()
