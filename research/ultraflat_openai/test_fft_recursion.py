"""Validate multidimensional FFT recurrence against exact rational Fourier coefficients."""
import unittest
import numpy as np
try:
    from .exact_pruning import step
    from .fft_recursion import original_poly_exact, convolve_dense, prune_new_symmetric
except ImportError:
    from exact_pruning import step
    from fft_recursion import original_poly_exact, convolve_dense, prune_new_symmetric

class TestFFTRecurrence(unittest.TestCase):
    def test_d3_to_d4_matches_original_exact(self):
        parent,exp=original_poly_exact(3)
        exact,next_exp=step(parent,exp,max_pairs=2_000_000)
        slab,meta=convolve_dense(parent,exp,limit=10000)
        coords=np.asarray(list(parent),dtype=int)
        lo=2*coords.min(axis=0)
        for s,b in exact.items():
            if s[-1]==0:
                continue
            got=slab[tuple(np.asarray(s[:-1])-lo)]
            self.assertAlmostEqual(got,float(b/(1<<next_exp)),places=10)
        exact_energy=float(sum((b/(1<<next_exp))**2 for b in exact.values()))
        self.assertAlmostEqual(meta["original_energy_estimate"],exact_energy,places=10)

    def test_budget_bounds_in_float(self):
        parent,exp=original_poly_exact(3)
        slab,meta=convolve_dense(parent,exp)
        info=prune_new_symmetric(slab,meta["old_energy"],len(parent),.005)
        self.assertLessEqual(info["deleted_L1_estimate"],.005+1e-12)
        self.assertLess(info["norm_one_rescaled_energy_estimate"],1.)
        self.assertGreater(info["total_support_estimate"],len(parent))

if __name__=="__main__":
    unittest.main()
