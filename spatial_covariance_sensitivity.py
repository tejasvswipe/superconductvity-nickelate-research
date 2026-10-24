#!/usr/bin/env python3
"""Independent robustness test; not a reproduction of the missing original generator.

Compares two positive Fourier spectra after matching their first-neighbour
correlation, with a padded 3N field cropped to the central film-sized window.
The dimensionless local-window rule and all random seeds are recorded in CSV.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT / "spatial_config.json").read_text())
N = int(CONFIG["N"])
PAD = int(CONFIG["padding_factor"] * N)
REPLICATES = int(CONFIG["replicates_per_condition"])
TARGET_RHOS = tuple(float(x) for x in CONFIG["target_nearest_neighbor_correlations"])
FAMILIES = tuple(CONFIG["covariance_families"])
MU = float(CONFIG["mean_proxy"])
X_OPT = float(CONFIG["optimum_proxy"])
SIGMA = float(CONFIG["spatial_standard_deviation"])
WINDOW_HALF_WIDTH = float(CONFIG["active_window_half_width"])
BASE_SEED = int(CONFIG["base_seed"])


def spectral_density(family: str, ell: float, n: int = PAD) -> np.ndarray:
    """Return a nonnegative isotropic Fourier spectrum on a periodic square."""
    q = 2.0 * np.pi * np.fft.fftfreq(n)
    qx, qy = np.meshgrid(q, q, indexing="ij")
    q2 = qx * qx + qy * qy
    if family == "white":
        return np.ones((n, n), dtype=float)
    if family == "gaussian_spectrum":
        return np.exp(-0.5 * ell * ell * q2)
    if family == "lorentzian_spectrum":
        return (1.0 + ell * ell * q2) ** -2.0
    raise ValueError(f"unknown family {family}")


def theoretical_rho1(family: str, ell: float) -> float:
    s = spectral_density(family, ell)
    c = np.fft.ifft2(s).real
    return float(c[0, 1] / c[0, 0])


def calibrate_ell(family: str, target_rho: float) -> float:
    if target_rho == 0.0:
        return 0.0
    lo, hi = 0.0, 1.0
    while theoretical_rho1(family, hi) < target_rho:
        hi *= 2.0
        if hi > PAD:
            raise ValueError(f"could not reach rho={target_rho} for {family}")
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        if theoretical_rho1(family, mid) < target_rho:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def make_field(rng: np.random.Generator, family: str, ell: float) -> np.ndarray:
    white = rng.standard_normal((PAD, PAD))
    spectrum = spectral_density(family, ell)
    z = np.fft.ifft2(np.fft.fft2(white) * np.sqrt(spectrum)).real
    c0 = PAD // 2 - N // 2
    z = z[c0:c0 + N, c0:c0 + N]
    z -= z.mean()
    sd = z.std(ddof=0)
    if not np.isfinite(sd) or sd == 0.0:
        raise FloatingPointError("degenerate field after crop")
    return z / sd


def nearest_neighbor_corr(z: np.ndarray) -> float:
    a = z[:, :-1].ravel()
    b = z[:, 1:].ravel()
    return float(np.corrcoef(a, b)[0, 1])


def spans(active: np.ndarray, diagonal: bool = False) -> bool:
    n = active.shape[0]
    seen = np.zeros_like(active, dtype=bool)
    if diagonal:
        offsets = ((-1, -1), (-1, 0), (-1, 1), (0, -1),
                   (0, 1), (1, -1), (1, 0), (1, 1))
    else:
        offsets = ((-1, 0), (1, 0), (0, -1), (0, 1))
    stack: list[tuple[int, int]] = []
    for r in range(n):
        if active[r, 0]:
            seen[r, 0] = True
            stack.append((r, 0))
    while stack:
        r, c = stack.pop()
        if c == n - 1:
            return True
        for dr, dc in offsets:
            rr, cc = r + dr, c + dc
            if 0 <= rr < n and 0 <= cc < n and active[rr, cc] and not seen[rr, cc]:
                seen[rr, cc] = True
                stack.append((rr, cc))
    return False


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = k / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return center - half, center + half


def main() -> None:
    seed_seq = np.random.SeedSequence(BASE_SEED)
    condition_seeds = seed_seq.spawn(len(TARGET_RHOS) * 2)
    rows = []
    replicate_rows = []
    condition_index = 0
    for family in FAMILIES:
        for target_rho in TARGET_RHOS:
            ell = calibrate_ell(family, target_rho)
            seed = int(condition_seeds[condition_index].generate_state(1)[0])
            condition_index += 1
            rng = np.random.default_rng(seed)
            span4 = span8 = 0
            active_total = 0
            rho_total = 0.0
            for replicate in range(REPLICATES):
                z = make_field(rng, family, ell)
                x = MU + SIGMA * z
                active = np.abs(x - X_OPT) <= WINDOW_HALF_WIDTH
                active_count = int(active.sum())
                rho_emp = nearest_neighbor_corr(z)
                did_span4 = int(spans(active, diagonal=False))
                did_span8 = int(spans(active, diagonal=True))
                active_total += active_count
                rho_total += rho_emp
                span4 += did_span4
                span8 += did_span8
                replicate_rows.append({
                    "family": family,
                    "target_rho_nn": target_rho,
                    "ell_grid_units": ell,
                    "condition_seed": seed,
                    "replicate_index_zero_based": replicate,
                    "empirical_crop_rho_nn": rho_emp,
                    "active_site_count": active_count,
                    "active_fraction": active_count / (N * N),
                    "spans_four_neighbor": did_span4,
                    "spans_eight_neighbor": did_span8,
                })
            p4, p8 = span4 / REPLICATES, span8 / REPLICATES
            lo4, hi4 = wilson(span4, REPLICATES)
            lo8, hi8 = wilson(span8, REPLICATES)
            rows.append({
                "family": family,
                "target_rho_nn": target_rho,
                "spectral_scale_ell_grid_units": ell,
                "theoretical_rho_nn": theoretical_rho1(family, ell),
                "empirical_crop_rho_nn_mean": rho_total / REPLICATES,
                "N": N,
                "padded_N": PAD,
                "replicates": REPLICATES,
                "seed": seed,
                "mean_proxy": MU,
                "sigma": SIGMA,
                "window_half_width": WINDOW_HALF_WIDTH,
                "mean_active_fraction": active_total / (REPLICATES * N * N),
                "span4_count": span4,
                "span4_probability": p4,
                "span4_wilson95_low": lo4,
                "span4_wilson95_high": hi4,
                "span8_count": span8,
                "span8_probability": p8,
                "span8_wilson95_low": lo8,
                "span8_wilson95_high": hi8,
            })
    out = ROOT / "spatial_covariance_sensitivity.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    raw_out = ROOT / "spatial_covariance_sensitivity_replicates.csv"
    with raw_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=replicate_rows[0].keys())
        writer.writeheader()
        writer.writerows(replicate_rows)
    print(f"Wrote {out} and {raw_out} ({len(rows)} conditions; {len(replicate_rows)} fields)")


if __name__ == "__main__":
    main()
