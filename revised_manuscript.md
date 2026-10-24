# Covariance and adjacency sensitivity in a correlated-field threshold-percolation model

**Working manuscript draft for scientific feedback.** Author list and affiliations should be retained from the submitted draft only after all authors approve this revised scope. This draft reports the new sensitivity computation separately from the user-supplied percolation manuscript's unreproduced numerical claims.

## Abstract

A supplied manuscript proposes that spatial correlation can alter the probability of a film-spanning favorable region even when the mean of a local scalar proxy and the favorable-site fraction are similar. Its original percolation implementation and raw fields were unavailable for independent reproduction. We therefore implemented a separately specified sensitivity experiment to examine two sources of model dependence: covariance shape and adjacency convention. A Gaussian random field was generated spectrally on a padded 192 × 192 lattice, cropped to 64 × 64, recentered and rescaled, and thresholded within a symmetric proxy window. Two nonnegative isotropic spectra were calibrated to the same full-grid nearest-neighbour correlations. Three hundred fields were sampled for each correlation target in each spectrum. Four-neighbour and eight-neighbour connectivity were measured on the same maps. At target nearest-neighbour correlation 0.80, four-neighbour spanning probabilities were 0.203 (61/300; Wilson 95% interval 0.162–0.252) for the Gaussian spectrum and 0.150 (45/300; interval 0.114–0.195) for the Lorentzian-squared spectrum; the intervals overlap. At zero target correlation, eight-neighbour spanning was 1.000 and 0.990 for the respective spectra and remained high throughout. Thus connectivity conclusions are strongly convention-dependent in this model. These computations establish neither reproduction of the original manuscript's maps nor spatially percolative superconductivity in nickelates. A material-specific next question is whether calibrated spatial reduction-state statistics improve held-out-batch prediction of transport-plus-magnetic response beyond mean proxy, precursor composition, phase fraction and strain.

## 1. Introduction

The supplied manuscript considers a binary spatial model motivated by reduction heterogeneity in infinite-layer nickelate films. It labels each lattice site favorable when a scalar local proxy lies within a finite window around an assumed optimum, then estimates the probability that favorable sites connect the left and right edges. The manuscript reports a strong difference between nominal correlation lengths at similar mean and active fraction. However, the submitted generator, exact covariance, raw realizations, and standalone exported data were not supplied. Its headline simulation cannot therefore be independently rerun from the available inputs.

The broad concept that local superconducting regions may acquire macroscopic response through a percolative process is not new. Pelc et al. modeled superconducting precursors in cuprates with inhomogeneous local transition temperatures and patch percolation [1]. Jang et al. discussed percolative paths in a thickness-driven cuprate film transition, while Tan et al. discussed superconducting islands and spanning superclusters in disordered indium-oxide films [2] [3]. Those studies do not establish the specific proxy, covariance, or mechanism in nickelates. The present numerical work therefore addresses a narrower issue: how much the geometry result depends on covariance and adjacency choices in a clearly specified threshold model.

## 2. Model and methods

### 2.1 Local window and random fields

For site `i`, define the dimensionless proxy `x_i = μ + σ z_i`, with per-map normalized field `z_i`. A site is active if

`A_i = 1{ |x_i − x*| ≤ w }`.

For an unstandardized unit-normal marginal, the one-site active probability is

`p = Φ((x* + w − μ)/σ) − Φ((x* − w − μ)/σ)`.

At `μ = x*`, this becomes `p = erf[w/(sqrt(2)σ)]`. For the chosen `σ = 0.13` and `w = 0.08`, it equals 0.461699. Each simulated crop was recentered and rescaled to exact sample mean zero and unit variance before `μ` and `σ` were applied, so the empirical active fraction fluctuates around, rather than being algebraically fixed to, this marginal value.

Fields were generated on a 192 × 192 periodic mesh using independent standard-normal white noise in Fourier space, multiplied by the square root of a nonnegative isotropic spectrum, inverse transformed, and cropped centrally to 64 × 64. The crop was then standardized. The compared spectra were `S_G(q) = exp(−ell²|q|²/2)` and `S_L(q) = [1 + ell²|q|²]⁻²`. For each family, `ell` was calibrated by bisection to target full-grid nearest-neighbour correlation `ρ1 ∈ {0, 0.35, 0.60, 0.80}`. These spectra are alternative mathematical choices; they are not calibrated nickelate covariance models. Matching one neighbour correlation does not match the full correlation function.

The proxy settings were `μ = x* = 0.50`, `σ = 0.13`, and `w = 0.08`. For each of eight spectrum-correlation conditions, 300 independent fields were generated. Seeds were deterministically derived from `SeedSequence(53821)` and recorded in the CSV. Spanning was tested horizontally under both four-neighbour (edge-sharing) and eight-neighbour (including diagonals) adjacency.

### 2.2 Statistical estimator and tests

For `R` maps, the estimated spanning probability is `p̂ = k/R`, where `k` maps span. Wilson 95% intervals are reported for binomial sampling uncertainty conditional on independent maps and the selected field generator. These intervals omit uncertainty from the covariance family, adjacency rule, finite map size, and mapping to a physical system.

The breadth-first-search spanning function was checked against an independently written union-find implementation on 20 seeded binary masks under both adjacency conventions. Additional tests checked empty and full limits, simple path and broken-path cases, covariance calibration, exact crop standardization, Wilson interval bounds, and the white-noise marginal active fraction. All eight automated tests passed. These tests verify code invariants, not physical validity or convergence with lattice size.

## 3. Results

### 3.1 Active fraction and empirical correlation

Across conditions, the simulated active fractions ranged from 0.46035 to 0.46181, close to the unstandardized marginal value 0.461699. The mean empirical nearest-neighbour correlations after crop and per-map standardization were close to the targets. Exact values, calibration scales, seeds, replicate counts, spanning counts, and Wilson intervals are reported in the machine-readable data file.

### 3.2 Four-neighbour spanning

For target correlations 0 and 0.35, neither spectrum produced any spanning maps under four-neighbour connectivity. At target 0.60, the Gaussian-spectrum family produced 0/300 spans and the Lorentzian-squared family produced 1/300. At target 0.80, the Gaussian-spectrum family produced 61/300 spans (`p̂ = 0.203`, 95% Wilson interval 0.162–0.252), and the Lorentzian-squared family produced 45/300 (`p̂ = 0.150`, interval 0.114–0.195). The intervals overlap, so this run does not establish a difference between the two covariance families. It does show that, in both selected families, the high-correlation condition lies in a partial-spanning regime whereas lower tested targets do not.

### 3.3 Adjacency convention

Eight-neighbour connectivity gave 300/300 spans for the Gaussian-spectrum case and 297/300 for the Lorentzian-squared case at zero target correlation. Spanning stayed high across the tested targets, reaching 288–300/300. The active maps are unchanged between these calculations. The difference arises only from whether diagonal contact counts as connectivity, demonstrating that the adjacency convention substantially changes the modeled conclusion.

## 4. Discussion

This sensitivity run is consistent with a general mathematical statement: spatial organization and graph connectivity can affect spanning probability while marginal site fraction is held approximately constant. It also shows that one scalar correlation length cannot fully specify spatial geometry. Two covariances matched at one nearest-neighbour statistic can differ in other spatial correlations, and eight-neighbour connectivity makes spanning much easier than edge-sharing connectivity in this parameter regime.

These are model results, not material-specific findings. The dimensionless proxy has no demonstrated calibration to oxygen concentration, nickel valence, local structure, or local superconducting response. The calculation contains no reduction kinetics, electronic states, Josephson coupling, phase stiffness, conductance, resistance, magnetic screening, or current. A spanning binary component therefore is not a superconducting transport path.

The supplied source manuscript's reported 1-to-3 correlation-length contrast and its reported padded-field sensitivity values remain manuscript-reported, not reproduced here. The current run used a newly specified spectral generator and matched nearest-neighbour targets; it must not be substituted for the missing original implementation. It also used only one map size and 300 maps per condition. Finite-size scaling, seed-batch replication, full-covariance matching, and exact-active-fraction control are still needed before claiming robust model generality.

Prior literature has already used percolation to discuss the emergence of superconductivity from local inhomogeneity in cuprates and disordered films [1] [2] [3]. The strongest distinct nickelate research direction is thus not generic geometric percolation. It is a prospective prediction test: whether independently mapped local reduction-state variance or correlation statistics add reproducible predictive value for transport-plus-diamagnetic outcomes across independent growth batches, beyond measured mean state and batch covariates. Existing nickelate work documents synthesis, stoichiometry, reduction and transport in specific platforms, but does not establish that predictive claim [4] [5].

## 5. Limitations and planned decisive tests

The next computational study should vary `σ`, `w`, target correlation, covariance family and lattice size in a preregistered design, first one factor at a time and then factorially. It should measure and compare complete empirical covariance functions, match active fraction as a separate control, vary boundary treatment, and replicate the full run using independent seed streams and a second field-generation method. A proposed convergence criterion is a change in spanning probability below 0.02 over two successive lattice-size doublings with a Monte Carlo Wilson half-width no greater than 0.03. This criterion is a plan and has not been met by the present 64 × 64 run.

For a material-level test, use independent growth batches as the statistical units. Measure precursor cation composition, strain and structural phase alongside a locally resolved reduction proxy whose calibration, resolution and uncertainty are reported. Pair transport with independent diamagnetic screening and control for current, filtering, contacts and substrate. Compare held-out-batch prediction using established covariates against the same model plus preregistered spatial statistics. A predictive improvement would establish a useful predictor, not by itself a causal or percolative mechanism. Such synthesis, cryogenic and spatial-characterization work requires qualified facilities.

## 6. Conclusion

In a newly specified correlated-field threshold model, four-neighbour spanning becomes partial at the highest tested matched nearest-neighbour correlation for both of two covariance spectra, while eight-neighbour spanning is already near-certain at zero correlation. The result demonstrates sensitivity to spatial and graph assumptions. It does not independently reproduce the supplied manuscript's original numbers and does not show percolative superconductivity in nickelates. Further computational robustness and material-calibrated, independent-batch prediction tests are required.

## Data and code availability

The code, fixed parameters, derived CSV, test suite, plot script, and test output are included in the accompanying local reproducibility package. The original percolation implementation and raw replicate fields were not provided with the supplied inputs. The tight-binding notebook and its fresh reproduction are preserved separately in the same package. No external repository or DOI has been created.

## References

[1]: https://doi.org/10.1038/s41467-018-06707-y "Emergence of superconductivity in the cuprates via a universal percolation process"
[2]: https://doi.org/10.1038/s41598-020-60037-y "Film-thickness-driven superconductor to insulator transition in cuprate superconductors"
[3]: https://doi.org/10.1103/PhysRevB.78.014506 "Evidence of spatially inhomogeneous pairing on the insulating side of a disorder-tuned superconductor-insulator transition"
[4]: https://doi.org/10.3389/fphy.2021.719534 "Impact of cation stoichiometry on the crystalline structure and superconductivity in nickelates"
[5]: https://doi.org/10.1002/advs.202309092 "Toward reliable synthesis of superconducting infinite layer nickelate thin films by topochemical reduction"
