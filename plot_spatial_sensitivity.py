#!/usr/bin/env python3
"""Plot prespecified spatial-model sensitivity data with Wilson 95% intervals."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
rows = list(csv.DictReader((ROOT / "spatial_covariance_sensitivity.csv").open()))
colors = {"gaussian_spectrum": "#1f77b4", "lorentzian_spectrum": "#d95f02"}
labels = {"gaussian_spectrum": "Gaussian spectrum", "lorentzian_spectrum": "Lorentzian spectrum"}
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), sharey=True, constrained_layout=True)
for ax, neighbors in zip(axes, (4, 8)):
    for family in colors:
        data = [r for r in rows if r["family"] == family]
        x = np.array([float(r["target_rho_nn"]) for r in data])
        p = np.array([float(r[f"span{neighbors}_probability"]) for r in data])
        low = np.array([float(r[f"span{neighbors}_wilson95_low"]) for r in data])
        high = np.array([float(r[f"span{neighbors}_wilson95_high"]) for r in data])
        yerr = np.maximum(0.0, np.vstack((p-low, high-p)))
        ax.errorbar(x, p, yerr=yerr, marker="o", lw=1.8,
                    capsize=3.5, color=colors[family], label=labels[family])
    ax.set_title(f"{neighbors}-neighbour connectivity")
    ax.set_xlabel("Target nearest-neighbour field correlation")
    ax.set_xlim(-0.04, 0.84)
    ax.set_ylim(-0.04, 1.04)
    ax.grid(alpha=0.25)
axes[0].set_ylabel("Left-to-right spanning probability")
axes[1].legend(frameon=False, loc="center right")
fig.suptitle("Connectivity depends on covariance and adjacency\n"
             r"$N=64$, 300 maps/condition, mean $=x^*=0.50$, $\sigma=0.13$, $w=0.08$")
out = ROOT / "spatial_covariance_sensitivity.png"
fig.savefig(out, dpi=200)
print(out)
