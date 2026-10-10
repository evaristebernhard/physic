"""Adaptive convolution unit tests: exact algebra and numerical d6 smoke test."""
import unittest
import numpy as np
try:
    from .adaptive_fft import initial_exact_pruned, next_step_fft
except ImportError:
    from adaptive_fft import initial_exact_pruned, next_step_fft

class TestAdaptiveRecursion(unittest.TestCase):
    def test_recurrence_lipschitz_identity(self):
        xs=np.linspace(-1.,1.,21)
        errors=np.linspace(-0.1,0.1,9)
        for p in xs:
            for e in errors:
                q=p+e
                for cosine in (-1.,-.2,0.,.7,1.):
                    tp=p+.5*(1-p*p)*cosine
                    tq=q+.5*(1-q*q)*cosine
                    self.assertAlmostEqual(
                        tq-tp,
                        (q-p)*(1-(q+p)*cosine/2),
                        places=12)
                    self.assertLessEqual(
                        abs(tq-tp),(2+abs(e)/2)*abs(e)+1e-12)

    def test_fifth_to_sixth_layer(self):
        p,e,info=initial_exact_pruned(5)
        self.assertEqual(info["raw_terms"],18362)
        self.assertEqual(info["retained_terms"],5994)
        self.assertLess(e,.005)
        q,res=next_step_fft(p,.0005,max_cells=700000,max_modes=200000)
        self.assertEqual(res["status"],"COMPUTED_FLOAT64")
        self.assertEqual(len(q),res["output_terms"])
        self.assertLess(res["dropped_L1_float"],.0005+1e-12)
        self.assertAlmostEqual(res["energy_float"],.48947955,places=5)
        E6=(2+e/2)*e+res["dropped_L1_float"]
        self.assertLess(E6,.001)

if __name__=="__main__":
    unittest.main()
