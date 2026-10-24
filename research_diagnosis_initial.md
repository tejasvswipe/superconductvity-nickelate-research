# Initial Research Diagnosis: Nickelate Percolation and Interface Charge Models

**Prepared 23 October 2026**

## Scope and status

This report covers the requested initial phases: inventory of both models, execution and red-team audit of the supplied tight-binding notebook, an evidence check of the supplied percolation manuscript/data, and the next-stage research questions. It preserves the two projects as distinct. The uploaded notebook was executed locally; Google Colab was not accessed.

The link to the Google Doc (`1nbDdd5xZNb01dwrUw7GD0v0QguuA2g3LSAg1HCxqmcE`) could not be read because the configured Google Docs CLI had no credentials. No login was initiated. The instructions describe an experimental protocol in that Doc, but its exact contents are therefore not part of this audit. Please provide its text or an accessible export for a complete protocol review.

Only three user-uploaded files were available: the percolation manuscript and pasted data in `pasted_content.txt`, research instructions in `pasted_content_2.txt`, and `nickelate_interface_model_colab.ipynb`. The pasted manuscript includes selected raw table values, but there is no standalone percolation source code, CSV, or independent padded-field implementation. The notebook states that it emits data files; those were not separate user uploads. The code below therefore distinguishes freshly reproduced notebook calculations from percolation numbers that remain reported but unreproduced.

## Baseline inventory

### Project A — correlated-field percolation

**Question and intended physical system.** The manuscript asks whether a spatial correlation length can change the probability of a film-spanning favorable region even when the mean dimensionless reduction proxy and favorable-site fraction are similar. It uses infinite-layer nickelate films only as motivation; the calculation itself is not parameterized to a material.

**Model and equations.** At lattice site *i*, the local proxy is `x_i = μ + σ z_i`, where `z_i` is a normalized Gaussian random field. A site is marked favorable when `|x_i − x*| ≤ w`. With four-neighbour site connectivity, an outcome is one if a favorable cluster touches both left and right edges. For `R` independent maps, the estimated spanning probability is `p̂ = (1/R) Σ I_r`. The report gives Wilson binomial intervals.

At the central condition `μ = x* = 0.50`, `σ = 0.13`, and `w = 0.08`, an untruncated unit-normal marginal implies an infinite-sample favorable-site fraction `erf[w/(√2 σ)] = 0.461699`. That is consistent with the reported finite-map means 0.46129 and 0.45841. It does not itself verify the random-field generator or finite-sample normalization.

**Variables and reported settings.** The manuscript varies `μ`, `σ`, correlation length, window width, and lattice size. Its primary result uses 64 × 64 maps, 500 realizations per condition, seed 247119, `μ = x* = 0.50`, `w = 0.08`, and correlation lengths 1, 2, and 3. The reported `σ = 0.13` contrast is `Pspan = 0.150` at length 1 versus 0.988 at length 3, with active fractions 0.46129 and 0.45841. I independently recalculated Wilson intervals from 75/500 and 494/500 outcomes; they match the supplied intervals 0.12138–0.18396 and 0.97407–0.99449. This verifies interval arithmetic only, not the underlying simulations.

**Boundary/numerical/statistical choices.** The text says Gaussian Fourier filtering on a 3N × 3N field followed by a central N × N crop, then map-wise recentering and rescaling. It uses open-edge left-to-right spanning and edge-sharing connectivity. The exact spectral filter, covariance function, mapping from filter parameter to stated correlation length, implementation, and per-realization samples are unavailable. The crop is described as reducing wraparound, but the separate padded implementation is also unavailable.

**What is supported and unsupported.** The pasted table and manuscript are internally consistent for the headline numerical contrast and active fraction. Conditional on the stated stochastic construction, they report a large geometry effect. They do not establish a reduction-state optimum, superconductivity at a site, a superconducting transport path, a material-specific correlation length, or improved experimental prediction. They also do not let us reproduce or independently validate any spanning probability.

### Project B — one-orbital layered tight-binding slab

**Question and intended physical system.** The notebook asks how a local interface potential, a weakened boundary hopping, a two-layer potential proxy, or uniform hopping changes layer occupations in a finite layered band model, and whether an exponential profile describes the resulting density change. It is an illustrative model, not a material-specific nickelate calculation.

**Hamiltonian and numerical solution.** For each in-plane momentum `k`, the slab matrix is

`H_ll'(k) = [ε(k) + V_l] δ_ll' − t_l δ_(l+1,l') − t_(l−1) δ_(l−1,l')`,

with `ε(k) = −2t(cos kx + cos ky)`. The interlayer direction is a finite open chain; the in-plane Brillouin zone is integrated on a uniform `K × K` mesh. The matrix is diagonalized with `numpy.linalg.eigh`. If its eigenvectors are `u_klm` and eigenvalues are `E_km`, the layer occupation including spin degeneracy is

`n_l = (2/Nk) Σ_km |u_klm|² f[(E_km − μ)/T]`,

where `f(y) = 1/(1 + exp(y))`. For the fixed-density calculation, the code bisects `μ` until `(1/L)Σ_l n_l` equals the requested filling. The reported observable is `Δn_l = n_l(perturbed) − n_l(reference)`.

**Settings and perturbations.** The code uses `t = 1`, `t_perp/t = 0.20`, `T/t = 0.04`, fillings 0.70, 1.00, and 1.30, thicknesses 4, 6, 8, and 12 layers, and a 45 × 45 main mesh. It tests `V_0 = ±0.30t`, first-bond hopping 0.05t instead of 0.20t, `V_0 = V_1 = 0.15t`, a uniform 5% increase in both in-plane and interlayer hopping, and an interface potential under that uniform hopping control. The code recalculates `μ` after each perturbation; this is a fixed-total-filling (canonical) comparison, not fixed chemical potential.

**Boundary conditions, normalization, and seeds.** The in-plane momentum sum is a periodic Brillouin-zone quadrature; the layer direction has open ends. Occupation includes a factor of two for spin and is normalized by the number of sampled momenta. The model is deterministic, so random seeds are not applicable. The notebook reports no calibration to measured band parameters.

## Notebook execution and independent checks

The principal code cell was extracted and executed locally. Its representative output matched the notebook's recorded output, including `Δn_0 = −0.126066`, profile RMS `0.047690`, and summed `Δn = −8.88 × 10⁻¹⁶` for eight layers at filling 1.0 under `V_0 = +0.30t`. It also reproduced the six-point k-grid table. The display-only code cell initially lacked IPython in the environment; after installing IPython in an isolated audit directory, it loaded the generated image successfully. Markdown cells were not executable. The run did not occur in a Google Colab account.

The six independent audit tests pass: nearest-neighbour dispersion bounds and finiteness; Hermiticity; fixed-filling normalization for all three fillings; density bounds; Fermi-function stability at `T/t = 0.04` and a near-zero-temperature limit; and representative momentum-grid sensitivity. The tests and outputs are in the linked reproducibility bundle.

### Findings that materially change interpretation

1. **The Hamiltonian is Hermitian.** The real symmetric hopping matrix passed a direct `H − H†` check with maximum absolute residual zero. The implemented square-lattice dispersion is the standard nearest-neighbour form. `eigh` supplies orthonormal eigenvectors and real eigenvalues.

2. **The chemical potential is solved afresh for every scenario.** The bisection has fixed bracketing and 70 iterations, so there is no user-supplied initial chemical potential whose value can contaminate results. At fixed density, all tested fillings were recovered to machine precision. Fermi-function arguments are clipped to ±80, avoiding overflow in the tested range.

3. **Charge conservation is a constraint, not an independent screening result.** For each perturbed calculation, solving `μ` to preserve total filling mathematically enforces `Σ_l Δn_l = 0` up to roundoff. It is not produced by post-hoc renormalization, but it is still a consequence of the ensemble constraint.

4. **Fixed chemical potential gives a different answer.** At eight layers and target filling 1.0 with `V_0 = +0.30t`, the canonical run shifts `μ/t` from `5.77 × 10⁻⁶` to 0.03841 and gives `Σ_l Δn_l ≈ 0`. Holding `μ` at the reference value instead changes the mean filling to 0.982516 and gives `Σ_l Δn_l = −0.139869`. Across the tested fillings, this distinction is material for onsite-potential scenarios. It must be stated whenever interpreting charge transfer.

5. **The reported profile is not a screened exponential.** For the eight-layer, filling-1 interface-potential case, the first layer changes by −0.1261, while the remaining seven layers each gain roughly +0.015 to +0.022 electrons per site. The profile is a depleted boundary layer plus nearly uniform compensation over the rest of this short slab, not a decay toward a bulk-like zero. The notebook's log-linear fit excludes the interface layer, takes `log |Δn|`, and returns a nominal length of about 68 layers with `R² = 0.068`. That length is not identifiable as a screening length. The fit also hides signs and can misread oscillatory or plateau profiles.

6. **The tiny first-bond response has a symmetry explanation.** At filling 1 in this bipartite nearest-neighbour model, a hopping-only change preserves particle-hole symmetry. Layer densities remain nearly one per site, so the barrier response is at numerical-noise scale. The output does not justify a general claim that interface hopping is unimportant in doped or multiorbital nickelates.

7. **The k-grid check is narrow.** For the one reported N = 8, filling-1, `V_0 = 0.30t` observable, the interface change shifts from −0.126066 at 45 × 45 to −0.125817 at 75 × 75 (about 0.20%); profile RMS shifts from 0.047690 to 0.047599 (about 0.19%). This supports local numerical stability for those observables and that condition. It does not demonstrate convergence for every filling, thickness, scenario, fit, or derived quantity. The 25 × 25 value differs by about 3.9% at the interface, so the selected mesh is not interchangeable with a coarse mesh.

8. **There is no self-consistent electrostatics.** `V_l` is an imposed one-body potential; the electron density does not feed back through Poisson/Hartree electrostatics. The model contains no long-range Coulomb interaction, dielectric boundary conditions, substrate bands, multiorbital hybridization, disorder, oxygen chemistry, or correlated electronic state. It computes response to an externally imposed perturbation, not the screening of that perturbation. Hence the current code cannot establish a material screening length.

The fit's low `R²` is correctly reported for its selected representative case, but its interpretation must be narrower: this profile is not well described by the selected log-linear `|Δn|` fit in this finite noninteracting slab. It is not evidence that exponential screening is absent in nickelates or in a model that includes self-consistent electrostatics.

## Red-team status and immediate corrections

| Issue | Why it matters | Severity | Minimum correction |
|---|---|---|---|
| Percolation source code and separate raw outputs are absent | The headline stochastic result cannot be reproduced, and the covariance/filter implementation cannot be audited | Fatal to a reproducibility claim; repairable if files exist | Supply source, configuration, raw replicate results, and the padded-field implementation; otherwise label all values as manuscript-reported |
| Google Doc protocol is inaccessible here | The full experimental design and endpoint definitions cannot be checked | Repairable access gap | Provide an export or paste the text |
| Correlation covariance is underspecified | Nominal correlation length does not uniquely determine cluster geometry | Fatal to physical inference; repairable mathematically | Define covariance and empirical estimator; test alternative covariance families while matching active fraction |
| Canonical conservation is described as a model result without ensemble contrast | A zero-sum profile can be mistaken for electrostatic screening | Repairable interpretation issue | Report the fixed-density constraint explicitly and add fixed-μ as a control, with ensemble chosen from physical boundary conditions |
| The one-orbital model has no electrostatic feedback | It cannot calculate screened potential or identify physical screening length | Fatal to a screening claim; repairable only by changing model | Add a self-consistent electrostatic kernel and boundary conditions; validate parameters and keep a no-Coulomb baseline |
| The reported decay fit is underpowered and sign-blind | A large `xi` with low fit quality can look like a length despite no identifiable exponential | Repairable analysis issue | Predefine fitting range, inspect signed profiles, compare candidate models and residuals, and report “no single length” when appropriate |
| Convergence is tested at one filling/thickness/scenario | The rest of the parameter sweep may remain grid-sensitive | Repairable | Define a numerical tolerance, then check all primary observables across meshes for representative edge and interior conditions |
| `uniform_hopping_plus5pct` changes both bandwidths but keeps absolute T fixed | It also changes reduced temperature `T/t` and is not a pure strain test | Repairable | State the control as a toy bandwidth change or scale T proportionally if fixed `T/t` is intended |
| Spatial percolation is described near superconductivity | A spanning favorable-site component is not a superconducting current path | Fatal to material-level claim; repairable by strict claim boundary | Keep it as a mathematical proxy unless calibrated conductance/Josephson physics and experimental maps are added |

No justified physical mechanism presently couples the two projects. Both remain separate research directions.

## Reproducibility and deliverables currently available

The audit bundle contains the unchanged inputs, extracted notebook cell, fresh CSV/PNG outputs, execution logs, a six-test `unittest` suite, a fixed-density/fixed-μ comparison, and environment notes. Notebook baseline execution succeeded. The uploaded percolation calculations remain unverified because their code and raw outputs were not supplied. The linked Doc remains inaccessible.

The literature review, competing-hypothesis ranking, and final research-direction recommendation will be added after source-by-source verification. A technically revised manuscript is intentionally not drafted before that review; the supplied instructions say to establish the evidence and research question first.
