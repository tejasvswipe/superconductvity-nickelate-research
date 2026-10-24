#!/usr/bin/env python3
"""Invariant and independent-implementation tests for the new sensitivity model."""
from __future__ import annotations
import unittest
import math
import numpy as np
import spatial_covariance_sensitivity as sp


def spans_union_find(active: np.ndarray, diagonal: bool = False) -> bool:
    """Independent union-find connectivity implementation used only as a cross-check."""
    n = active.shape[0]
    parent = np.arange(n * n + 2)
    rank = np.zeros(n * n + 2, dtype=np.int8)
    left, right = n * n, n * n + 1

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = int(parent[a])
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    offsets = ((-1, 0), (0, -1))
    if diagonal:
        offsets = ((-1, -1), (-1, 0), (-1, 1), (0, -1))
    for r in range(n):
        for c in range(n):
            if not active[r, c]:
                continue
            i = r * n + c
            if c == 0:
                union(i, left)
            if c == n - 1:
                union(i, right)
            for dr, dc in offsets:
                rr, cc = r + dr, c + dc
                if 0 <= rr < n and 0 <= cc < n and active[rr, cc]:
                    union(i, rr * n + cc)
    return find(left) == find(right)


class PercolationChecks(unittest.TestCase):
    def test_empty_and_full_limits(self):
        self.assertFalse(sp.spans(np.zeros((8, 8), bool)))
        self.assertTrue(sp.spans(np.ones((8, 8), bool)))

    def test_horizontal_path_and_disconnected_column(self):
        path = np.zeros((7, 7), bool)
        path[3, :] = True
        self.assertTrue(sp.spans(path))
        path[3, 3] = False
        self.assertFalse(sp.spans(path))

    def test_diagonal_rule_changes_connectivity(self):
        checker = np.eye(7, dtype=bool)
        self.assertFalse(sp.spans(checker, diagonal=False))
        self.assertTrue(sp.spans(checker, diagonal=True))

    def test_bfs_and_union_find_agree_on_seeded_masks(self):
        for seed in range(20):
            rng = np.random.default_rng(seed)
            mask = rng.random((17, 17)) < (0.35 + 0.01 * seed)
            for diagonal in (False, True):
                self.assertEqual(sp.spans(mask, diagonal), spans_union_find(mask, diagonal))

    def test_covariance_calibration(self):
        for family in ("gaussian_spectrum", "lorentzian_spectrum"):
            for target in (0.0, 0.35, 0.60, 0.80):
                ell = sp.calibrate_ell(family, target)
                self.assertAlmostEqual(sp.theoretical_rho1(family, ell), target, places=8)

    def test_field_normalization(self):
        rng = np.random.default_rng(92173)
        for family in ("gaussian_spectrum", "lorentzian_spectrum"):
            ell = sp.calibrate_ell(family, 0.6)
            z = sp.make_field(rng, family, ell)
            self.assertEqual(z.shape, (sp.N, sp.N))
            self.assertAlmostEqual(float(z.mean()), 0.0, places=12)
            self.assertAlmostEqual(float(z.std(ddof=0)), 1.0, places=12)
            self.assertTrue(np.isfinite(z).all())

    def test_wilson_interval_contains_sample_proportion(self):
        for k, n in ((0, 300), (1, 300), (150, 300), (300, 300)):
            lo, hi = sp.wilson(k, n)
            p = k / n
            self.assertLessEqual(lo, p + 1e-15)
            self.assertGreaterEqual(hi, p - 1e-15)
            self.assertLessEqual(0.0, lo)
            self.assertLessEqual(hi, 1.0 + 1e-15)

    def test_active_fraction_matches_marginal_order(self):
        rng = np.random.default_rng(331)
        total = 0
        reps = 100
        for _ in range(reps):
            z = sp.make_field(rng, "white", 0.0)
            x = sp.MU + sp.SIGMA * z
            total += int((np.abs(x - sp.X_OPT) <= sp.WINDOW_HALF_WIDTH).sum())
        empirical = total / (reps * sp.N * sp.N)
        expected = math.erf(sp.WINDOW_HALF_WIDTH / (math.sqrt(2) * sp.SIGMA))
        self.assertLess(abs(empirical - expected), 0.015)


if __name__ == "__main__":
    unittest.main(verbosity=2)
