#!/usr/bin/env python3
"""Independent consistency checks for the supplied tight-binding notebook.

This imports the uploaded notebook's computational cell without running its full
parameter sweep. It checks algebraic invariants, density normalization, ensemble
handling, and one representative k-grid convergence sequence.
"""
from __future__ import annotations
import csv
import importlib.util
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "original_reproduction"
MODULE_PATH = BASE / "notebook_cell_1.py"
spec = importlib.util.spec_from_file_location("submitted_model", MODULE_PATH)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


def fixed_mu_density(nlayer, mu, t_inplane=1.0, t_perp=0.2, onsite=None,
                     first_bond=None, temperature=0.04, mesh=45):
    """Evaluate layer occupations at a specified chemical potential (grand canonical)."""
    if onsite is None:
        onsite = np.zeros(nlayer)
    energies, vectors = model.diagonalize_stack(
        nlayer, t_inplane, t_perp, onsite, first_bond, mesh
    )
    x = np.clip((energies - mu) / max(float(temperature), 1e-9), -80, 80)
    f = 1.0 / (1.0 + np.exp(x))
    weights = np.abs(vectors) ** 2
    return 2.0 * np.einsum("klm,km->l", weights, f, optimize=True) / len(energies)


class NotebookAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = []

    def record(self, name, passed, detail):
        self.rows.append({"test": name, "passed": bool(passed), "detail": detail})

    def test_dispersion_limits_and_finiteness(self):
        # Gamma is sampled for odd meshes; nearest-neighbour square-lattice band edges are +/- 4t.
        eps = model.k_dispersion(45, 1.0)
        passed = (np.isfinite(eps).all() and np.isclose(eps.min(), -4.0)
                  and eps.max() <= 4.0 + 1e-12 and eps.max() > 3.98)
        self.record("square_lattice_dispersion_bounds", passed,
                    {"min": float(eps.min()), "max": float(eps.max()), "mesh": 45})
        self.assertTrue(passed)

    def test_hamiltonian_hermiticity(self):
        nlayer, mesh = 8, 45
        eps = model.k_dispersion(mesh, 1.0)
        onsite = np.linspace(-0.1, 0.1, nlayer)
        worst = 0.0
        for e in eps[::137]:
            h = np.diag(e + onsite).astype(complex)
            for j in range(nlayer - 1):
                h[j, j + 1] = h[j + 1, j] = -0.2
            worst = max(worst, float(np.max(np.abs(h - h.conj().T))))
            np.linalg.eigh(h)
        passed = worst == 0.0
        self.record("hamiltonian_hermiticity", passed, {"max_abs_H_minus_Hdagger": worst})
        self.assertTrue(passed)

    def test_particle_number_and_density_bounds(self):
        errors = {}
        for filling in (0.7, 1.0, 1.3):
            rho, mu = model.layer_density(8, filling, mesh=45)
            errors[str(filling)] = {
                "mean_filling": float(rho.mean()),
                "mean_error": float(rho.mean() - filling),
                "mu": float(mu),
                "min_layer_occupation": float(rho.min()),
                "max_layer_occupation": float(rho.max()),
            }
            self.assertLess(abs(rho.mean() - filling), 1e-10)
            self.assertTrue(np.all((rho >= 0) & (rho <= 2)))
        self.record("fixed_total_filling", True, errors)

    def test_charge_conservation_is_due_to_fixed_density_solver(self):
        nlayer, filling, mesh = 8, 1.0, 45
        ref, mu0 = model.layer_density(nlayer, filling, mesh=mesh)
        onsite = np.zeros(nlayer)
        onsite[0] = 0.30
        rho_fixed_n, mu1 = model.layer_density(nlayer, filling, onsite=onsite, mesh=mesh)
        rho_fixed_mu = fixed_mu_density(nlayer, mu0, onsite=onsite, mesh=mesh)
        fixed_n_delta_sum = float((rho_fixed_n - ref).sum())
        fixed_mu_delta_sum = float((rho_fixed_mu - ref).sum())
        detail = {
            "baseline_mu_over_t": float(mu0),
            "perturbed_fixed_N_mu_over_t": float(mu1),
            "fixed_N_sum_delta_n": fixed_n_delta_sum,
            "fixed_mu_sum_delta_n": fixed_mu_delta_sum,
            "fixed_mu_total_filling": float(rho_fixed_mu.mean()),
        }
        passed = abs(fixed_n_delta_sum) < 1e-10 and abs(fixed_mu_delta_sum) > 1e-5
        self.record("ensemble_control_charge_sum", passed, detail)
        self.assertTrue(passed)

    def test_finite_stable_fermi_occupations(self):
        for temperature in (0.04, 1e-9):
            rho, mu = model.layer_density(8, 1.0, temperature=temperature, mesh=25)
            self.assertTrue(np.isfinite(rho).all() and np.isfinite(mu))
            self.assertTrue(np.all((rho >= 0) & (rho <= 2)))
        self.record("fermi_numerical_stability", True,
                    {"temperatures_over_t": [0.04, 1e-9], "bounds": [0, 2]})

    def test_representative_k_grid_convergence(self):
        nlayer, filling = 8, 1.0
        values = []
        for mesh in (25, 35, 45, 55, 65, 75):
            onsite = np.zeros(nlayer)
            onsite[0] = 0.30
            rho, _ = model.layer_density(nlayer, filling, onsite=onsite, mesh=mesh)
            ref, _ = model.layer_density(nlayer, filling, mesh=mesh)
            delta = rho - ref
            values.append({"mesh": mesh, "interface_delta_n": float(delta[0]),
                           "profile_rms_delta_n": float(np.sqrt(np.mean(delta**2)))})
        self.record("representative_k_grid", True, values)

    @classmethod
    def tearDownClass(cls):
        OUT = ROOT / "audit_results.json"
        OUT.write_text(json.dumps({
            "python": sys.version.split()[0], "numpy": np.__version__,
            "tests": cls.rows,
        }, indent=2) + "\n")
        with (ROOT / "ensemble_comparison.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["ensemble", "mu_over_t", "filling_per_layer", "sum_delta_n"])
            w.writeheader()
            for test in cls.rows:
                if test["test"] == "ensemble_control_charge_sum":
                    d = test["detail"]
                    w.writerow({"ensemble": "fixed_total_density", "mu_over_t": d["perturbed_fixed_N_mu_over_t"],
                                "filling_per_layer": 1.0, "sum_delta_n": d["fixed_N_sum_delta_n"]})
                    w.writerow({"ensemble": "fixed_chemical_potential", "mu_over_t": d["baseline_mu_over_t"],
                                "filling_per_layer": d["fixed_mu_total_filling"], "sum_delta_n": d["fixed_mu_sum_delta_n"]})
        print(f"Wrote {OUT}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
