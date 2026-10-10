"""Exact-arithmetic regression tests for original p_d plus safe pruning."""
import unittest
from fractions import Fraction as Q

try:
    from .exact_pruning import step, prune, projection_gamma, audit
    from .audit import next_poly, energy
except ImportError:
    from exact_pruning import step, prune, projection_gamma, audit
    from audit import next_poly, energy

class TestExactPruning(unittest.TestCase):
    def test_dyadic_matches_original_fraction_recursion(self):
        a={}
        q=0
        old={}
        for d in range(1,5):
            a,q=step(a,q)
            old=next_poly(old)
            self.assertEqual(
                {s:Q(b,1<<q) for s,b in a.items()},
                old,
                f"Mismatch at depth {d}"
            )
            self.assertEqual(
                Q(sum(b*b for b in a.values()),1<<(2*q)),
                energy(old))

    def test_d5_exact_counts_and_energy(self):
        a={}
        q=0
        for _ in range(5): a,q=step(a,q)
        self.assertEqual(q,62)
        self.assertEqual(len(a),18362)
        exact_energy=Q(sum(v*v for v in a.values()),1<<(2*q))
        self.assertAlmostEqual(float(exact_energy),0.4354244759621249)
        retained,loss,loss_squared,scale,raw,clean=prune(a,q,Q(5,1000))
        self.assertEqual(len(retained),5994)
        self.assertLessEqual(loss,Q(5,1000))
        self.assertEqual(clean,(raw-loss_squared)/(1+loss)**2)
        self.assertEqual(scale,1+loss)
        self.assertLess(loss_squared,Q(5,1000)**2/len(a))
        self.assertTrue(all(tuple(-x for x in s) in retained for s in retained))

    def test_preprocessed_input_still_has_nonparallel_pairs(self):
        a={}
        q=0
        for _ in range(4): a,q=step(a,q)
        retained,*_=prune(a,q,Q(5,1000))
        representative=[s for s in retained if s>tuple(-v for v in s)]
        self.assertGreaterEqual(len(representative),2)
        for i,s in enumerate(representative):
            # All nonzero frequencies have final nonzero coordinate ±1,
            # so distinct sign-pairs are automatically nonparallel.
            last=next(v for v in reversed(s) if v)
            self.assertIn(last,(-1,1))
        self.assertTrue(all(sum(x*y for x,y in zip(s,projection_gamma(retained)))!=0
                            for s in retained))

    def test_budget_is_reported_not_confused_with_N(self):
        data=audit([4,5,6],Q(5,1000),max_pairs=2_000_000)
        self.assertEqual(data["runs"][2]["status"],"BUDGET")
        self.assertIn("NOT final N",data["runs"][0]["interpretation"])

if __name__=="__main__":
    unittest.main()
