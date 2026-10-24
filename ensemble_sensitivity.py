#!/usr/bin/env python3
"""Reproduce ensemble_sensitivity.csv from the supplied tight-binding model."""
from __future__ import annotations
import csv
import importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
source = ROOT / "original_reproduction" / "notebook_cell_1.py"
spec = importlib.util.spec_from_file_location("submitted_model", source)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


def density_at_fixed_mu(nlayer, mu, *, onsite, first_bond, mesh=45,
                        t_inplane=1.0, t_perp=0.2, temperature=0.04):
    energies, vectors = model.diagonalize_stack(
        nlayer, t_inplane, t_perp, onsite, first_bond, mesh
    )
    x = np.clip((energies - mu) / max(float(temperature), 1e-9), -80, 80)
    occupations = 1.0 / (1.0 + np.exp(x))
    weights = np.abs(vectors) ** 2
    return 2.0 * np.einsum("klm,km->l", weights, occupations, optimize=True) / len(energies)


def main():
    rows = []
    nlayer, mesh = 8, 45
    for filling in (0.7, 1.0, 1.3):
        reference, mu0 = model.layer_density(nlayer, filling, mesh=mesh)
        cases = []
        onsite = np.zeros(nlayer); onsite[0] = 0.30
        cases.append(("interface_onsite_plus0.30", onsite, None))
        cases.append(("interface_first_bond_0.05", np.zeros(nlayer), 0.05))
        onsite2 = np.zeros(nlayer); onsite2[:2] = 0.15
        cases.append(("two_layer_onsite_plus0.15", onsite2, None))
        for name, onsite, bond in cases:
            rho_n, mu_n = model.layer_density(
                nlayer, filling, onsite=onsite, first_bond=bond, mesh=mesh
            )
            rho_mu = density_at_fixed_mu(
                nlayer, mu0, onsite=onsite, first_bond=bond, mesh=mesh
            )
            for ensemble, rho, mu in (
                ("fixed_total_density", rho_n, mu_n),
                ("fixed_chemical_potential", rho_mu, mu0),
            ):
                delta = rho - reference
                rows.append({
                    "filling_target": filling,
                    "scenario": name,
                    "ensemble": ensemble,
                    "mu_over_t": float(mu),
                    "mean_filling_after": float(rho.mean()),
                    "sum_delta_n": float(delta.sum()),
                    "profile_rms_delta_n": float(np.sqrt(np.mean(delta**2))),
                    "interface_delta_n": float(delta[0]),
                })
    out = ROOT / "ensemble_sensitivity.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
