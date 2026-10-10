"""Algorithm-level regression checks; do not mistake numerical fits for proofs."""
import math
import unittest
import numpy as np
try:
    from .lowrank_balance import exact_level_two
    from .finite_n_search import sampled_polynomial, bernstein_grid_bounds, optimize
except ImportError:
    from lowrank_balance import exact_level_two
    from finite_n_search import sampled_polynomial, bernstein_grid_bounds, optimize

class TestFiniteAlgorithms(unittest.TestCase):
    def test_analytic_two_direction_exact_balance(self):
        result=exact_level_two()
        self.assertEqual(result["J"],2)
        self.assertEqual(result["D"],4)
        self.assertLess(result["max_relative_error_float"],1e-9)

    def test_exact_grid_flip_identity(self):
        a=np.array([1,-1,1,1,-1,1,-1,1], dtype=np.int8)
        n=len(a); M=128; j=3
        P=sampled_polynomial(a,M)
        z=np.exp(2j*math.pi*j*np.arange(M)/M)
        q=abs(P)**2/n-1
        predicted=q+(4/n)*(1-int(a[j])*np.real(np.conj(P)*z))
        new=a.copy();new[j]*=-1
        actual=abs(sampled_polynomial(new,M))**2/n-1
        self.assertTrue(np.allclose(predicted,actual,atol=1e-12))

    def test_grid_to_circle_bernstein_envelope(self):
        signs=np.array([1,-1,1,1,1,-1,-1,1,1,-1,-1,-1],dtype=np.int8)
        r=bernstein_grid_bounds(signs,grid_factor=32)
        p=sampled_polynomial(signs,20000)
        dense=abs(p)**2/len(signs)
        self.assertLessEqual(max(dense),r["analytic_max_square_upper"]+1e-10)
        self.assertGreaterEqual(min(dense),r["analytic_min_square_lower"]-1e-10)

    def test_short_search_outputs_binary_coefficients(self):
        x=optimize(12,steps=20,restarts=1,seed=123,sharpness=5)
        self.assertEqual(len(x["signs"]),12)
        self.assertEqual(set(x["signs"])-set("+-"),set())
        self.assertTrue(math.isfinite(x["grid_certificate"]["analytic_full_circle_relative_error_bound"]))

if __name__=="__main__":
    unittest.main()
