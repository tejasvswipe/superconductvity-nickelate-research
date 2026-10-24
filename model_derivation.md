# Mathematical Reconstruction of the Two Nickelate Research Models

## 1. Scope

The spatial-field percolation calculation and the layered tight-binding slab are distinct mathematical models. The first maps a correlated scalar field to a binary lattice-connectivity event. The second maps one-electron eigenstates to layer occupations under either a particle-number or chemical-potential constraint. No coupling between these models is supplied, and this derivation does not merge them.

All tight-binding quantities below are dimensionless in units of the in-plane hopping `t` unless noted. All spatial-field quantities are dimensionless model inputs, not measured oxygen concentrations or material constants.

## 2. Correlated-field site percolation

### 2.1 Field and local favorable-window rule

Let `i` index sites on an `N × N` square lattice. Define

`x_i = μ + σ z_i`,

where `μ` is a mean proxy, `σ ≥ 0` is its spatial standard deviation, and `z_i` is a unit-variance, zero-mean random field. The manuscript describes generating `z` by Gaussian Fourier filtering, then recentering and rescaling each finite map. It does not supply the filter kernel or the covariance function. To state what is and is not determined, write a generic covariance as

`C_λ(r) = E[z_i z_(i+r)]`, with `C_λ(0)=1`,

where `λ` is a nominal correlation-scale parameter. The exact relationship between `λ`, the Fourier filter, and a measured correlation length cannot be reconstructed from the supplied text alone.

For assumed optimum `x*` and window half-width `w`, define the binary field

`A_i = 1{ |x_i − x*| ≤ w }`.

If `z_i` has a standard-normal marginal before finite-map recentering/rescaling, then the infinite-sample one-site favorable probability is

`p = Φ((x* + w − μ)/σ) − Φ((x* − w − μ)/σ)`,

where `Φ` is the standard-normal CDF. At `μ=x*`, this reduces to

`p = erf(w/(√2 σ))`.

For the central manuscript inputs `w=0.08` and `σ=0.13`, this predicts `p=0.461699` before the finite-map normalization operation. The reported active fractions, approximately 0.46129 and 0.45841 for two nominal correlation lengths, are consistent with this marginal value. Correlation changes cluster geometry while leaving the one-site marginal approximately fixed; finite-sample recentering and rescaling can introduce small differences in empirical active fraction.

### 2.2 Connectivity event and estimator

Use four-neighbour adjacency: sites are adjacent only if they share a lattice edge. Let `S_r` be one if realization `r` contains at least one connected active component touching both the left and right boundaries, and zero otherwise. For `R` realizations,

`p̂_span = (1/R) Σ_(r=1)^R S_r`.

A Wilson interval for binomial proportion `p̂=k/R`, normal quantile `z`, is

`center = (p̂ + z²/(2R))/(1+z²/R)`,

`halfwidth = z sqrt[p̂(1−p̂)/R + z²/(4R²)]/(1+z²/R)`.

These intervals quantify finite Monte Carlo sampling error conditional on the model and independent-replicate assumption. They do not include uncertainty in `μ`, `σ`, `w`, covariance, field normalization, or the mapping to a material.

### 2.3 Assumptions and numerical checks required

The manuscript reports `N=64`, `R=500`, `μ=x*=0.50`, `w=0.08`, `σ=0.13`, and correlation scales 1 and 3, with `p̂_span=0.150` and `0.988`, respectively. The provided manuscript includes values and Wilson limits but not the field generator or replicate maps. Therefore the derivation establishes the meaning of the statistic, not the correctness of those simulations.

A reproducible implementation must report the spectral filter or `C_λ(r)`, field normalization, crop and edge conventions, random-seed stream, treatment of invalid/nonfinite maps, exact connectivity algorithm, and per-condition replicate count. Robustness should compare alternative covariance functions, exact or carefully matched active fractions, lattice sizes, and independent seeds. A spanning component is only a geometric event; it is not by itself a superconducting state or a conducting path.

## 3. Layered one-orbital tight-binding slab

### 3.1 Hamiltonian and boundary conditions

For an in-plane wavevector `k=(kx,ky)` and `L` layers, the notebook uses the nearest-neighbour square-lattice dispersion

`ε(k) = −2t(cos kx + cos ky)`.

The `L × L` slab Hamiltonian is

`H_ll'(k) = [ε(k)+V_l] δ_ll' − t_l δ_(l+1,l') − t_(l−1) δ_(l−1,l')`,

where `V_l` is an imposed onsite energy and `t_l` is the hopping on the bond between layers `l` and `l+1`. For the standard slab, `t_l=t_perp`; for the interface-barrier scenario, only the first bond is reduced. The layer direction is an open chain: there is no bond connecting layer `L−1` back to layer 0. The in-plane direction is represented by uniform Brillouin-zone quadrature on a `K × K` mesh.

With real `V_l` and `t_l`, `H(k)=H(k)^T=H(k)†`, so it is Hermitian. Diagonalize it as

`H(k) u_(km) = E_(km) u_(km)`,

with normalized eigenvectors `Σ_l |u_(klm)|²=1`.

### 3.2 Layer density and chemical potential

Let `T` denote thermal energy `k_B T_phys` in the same units as `t`; the notebook input `TEMPERATURE=0.04` is therefore `k_B T_phys/t`. Define `f_T(E−μ)=1/[1+exp((E−μ)/T)]`. With two spin states and `N_k=K²` sampled momenta,

`n_l(μ,T) = (2/N_k) Σ_(k,m) |u_(klm)|² f_T(E_(km)−μ)`.

The code targets filling `ν` per layer/orbital including spin by solving

`(1/L) Σ_l n_l(μ,T) = ν`.

The occupation is monotone increasing with `μ`; bracketed bisection is therefore appropriate. The notebook uses an energy-derived bracket and 70 iterations. It recalculates `μ` separately for every perturbation. The occupation calculation clips `(E−μ)/T` to `[-80,80]` before exponentiation, preventing overflow for the tested finite-temperature and near-zero-temperature checks.

### 3.3 Canonical versus grand-canonical response

For a perturbed potential or hopping configuration, define the response relative to an unperturbed reference:

`Δn_l = n_l^pert − n_l^ref`.

In the fixed-total-density calculation, both configurations are solved at the same target `ν`, so

`Σ_l n_l^pert = Lν = Σ_l n_l^ref`,

and consequently

`Σ_l Δn_l = 0`.

This follows from the imposed particle-number constraint, up to numerical root-finding and summation error. It is not an independent prediction of Coulomb screening. At fixed chemical potential, `μ` is held at the reference value and generally `Σ_l Δn_l ≠ 0`, because the system exchanges particles with a reservoir.

In the supplied calculation at `L=8`, `ν=1`, and `V_0=+0.30t`, the fixed-density solver changes `μ/t` from `5.77×10⁻⁶` to `0.03841` and preserves the layer-summed density to `−8.88×10⁻¹⁶`. Holding `μ` fixed gives mean filling 0.982516 and `Σ_l Δn_l=−0.139869`. These are two different ensembles and must not be conflated.

### 3.4 Why this Hamiltonian does not calculate screened electrostatics

`V_l` is an external one-body perturbation. The potential is not updated from `Δn_l`; there is no Poisson equation, Coulomb kernel, dielectric environment, substrate, or reservoir self-energy. Thus the calculation determines a band-electron density response to a prescribed perturbation. It cannot determine the screened potential or a material screening length.

A minimal next model would add an electrostatic feedback term, for example

`V_l^eff = V_l^ext + Σ_(l') U_(l,l') [n_(l') − n_(l')^ref]`,

where `U_(l,l')` is an explicitly specified Coulomb Green's-function matrix in energy units, including lattice spacing, dielectric environment, and boundary conditions. Solve the eigenproblem and density, update `V^eff`, and iterate until both the density/number residual and potential residual satisfy preregistered tolerances. The kernel, dielectric constants, fixed-charge or reservoir boundary condition, and neutrality rule must be stated. Without material-specific values this remains a model sensitivity study, not a nickelate prediction.

### 3.5 Decay-profile fitting

An exponential model for a signed response could be `Δn_l=A exp(−z_l/ξ)` only if the profile approaches a bulk baseline and remains sign-consistent over a justified range. If signs alternate, a candidate is a damped oscillation `Δn_l=A exp(−z_l/ξ) cos(qz_l+φ)`. A multi-component response can require more than one length. The notebook instead fits a straight line to `log |Δn_l|`, excludes the interface layer, and filters small values. This can hide sign changes and does not distinguish an exponential from a plateau, oscillation, or finite-size compensation.

For the representative eight-layer, half-filled `V_0=+0.30t` result, the interface changes by −0.1261 and the other layers gain about +0.015 to +0.022 each. A reported `ξ≈68` layers with `R²≈0.068` is therefore not an identifiable screening length. A defensible model comparison should inspect signed profiles first, preregister fitting windows, compare only physically justified models using residual diagnostics/information criteria, and return “no single length identifiable” when estimates vary materially with thickness or fit window.

## 4. Verification boundary

The notebook calculations establish a reproducible mathematical result for the stated one-orbital noninteracting model. They are not a material-specific result. The percolation manuscript's headline contrast is supported by its pasted table and internally consistent Wilson arithmetic, but its numerical implementation is not independently executable from the supplied files. Neither model establishes that the other is a physical description of infinite-layer nickelate films.
