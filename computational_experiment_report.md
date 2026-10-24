# Computational experiment report: correlated-field connectivity and the nickelate slab notebook

## What was actually run

Two supplied research directions were kept separate. The uploaded layered tight-binding notebook was executed from its computational cell in a clean output directory and its recorded representative outputs matched. An independent six-test audit checked Hamiltonian Hermiticity, dispersion bounds, occupations and fixed filling, the fixed-chemical-potential control, Fermi-function stability, and one k-grid sequence. A separate eight-condition spatial-field sensitivity test was newly implemented and run here because the manuscript's original generator and raw replicate maps were not provided. It is explicitly a new model experiment, not a reproduction of the manuscript's 500-map results or its separately reported padded-field test.

The supplied manuscript reports a 64 × 64 binary threshold-percolation model with mean proxy 0.50, assumed optimum 0.50, standard deviation 0.13 and half-window 0.08. It reports spanning probability 0.150 at nominal correlation length 1 and 0.988 at length 3 from 500 realizations per condition. Those values and Wilson-interval arithmetic are present in the supplied text, but their executable implementation, full parameter data and random maps were not supplied, so the results cannot be independently regenerated. The linked Google Doc protocol was not accessible in this session.

## New, explicitly defined spatial sensitivity experiment

### Model and parameters

For each 64 × 64 map, a mean-zero Gaussian field was generated spectrally on a 192 × 192 periodic square and the central 64 × 64 crop was retained. The cropped map was recentered and rescaled to sample mean 0 and standard deviation 1. Two nonnegative isotropic spectra were compared:

- Gaussian: `S(q) = exp(−ell² |q|² / 2)`.
- Lorentzian-squared: `S(q) = [1 + ell² |q|²]⁻²`.

For each spectrum, the scale parameter `ell` was numerically calibrated so the theoretical full-grid nearest-neighbour correlation equalled 0, 0.35, 0.60, or 0.80. Empirical correlations on the cropped and normalized maps were also recorded. The local proxy was `x_i = 0.50 + 0.13 z_i`; a site was active when `|x_i − 0.50| ≤ 0.08`. Each condition used 300 independent fields and a recorded child seed derived from base seed 53821. Four-neighbour and eight-neighbour left-to-right spanning were measured on every identical active map. Wilson 95% intervals describe replicate-sampling uncertainty only. The full parameters and all aggregate outputs are in `spatial_covariance_sensitivity.csv`; the fixed configuration is at the top of `spatial_covariance_sensitivity.py`.

The one-site Gaussian marginal predicts active fraction `erf(0.08/(sqrt(2)×0.13)) = 0.461699` before per-map finite-field standardization. Across the eight simulated conditions, empirical fractions ranged from 0.46035 to 0.46181. The actual means and standard deviations were imposed exactly per cropped map. Nearest-neighbour correlations after crop and standardization closely matched their targets.

### Results

For **four-neighbour** connectivity, both covariance families produced no spanning fields at target correlations 0, 0.35, and 0.60, except one of 300 Lorentzian-squared maps at 0.60. At target correlation 0.80, the Gaussian-spectrum case spanned in 61/300 maps, `p = 0.203`, Wilson 95% interval 0.162–0.252; the Lorentzian-squared case spanned in 45/300 maps, `p = 0.150`, interval 0.114–0.195. Those final intervals overlap. Thus the tested covariance families agree qualitatively that the high-correlation condition moves into a partial-spanning regime, but this small design does not resolve a between-family difference in that regime.

For **eight-neighbour** connectivity, the baseline uncorrelated fields already spanned at high rates: 300/300 and 297/300 for the Gaussian and Lorentzian-squared spectra, respectively. Across all tested correlations, eight-neighbour spanning probabilities ranged from 0.960 to 1.000. The same fields therefore give sharply different connectivity conclusions under two adjacency rules. The local fraction does not change when adjacency changes. This is direct evidence that the connectivity convention is a consequential assumption, not a cosmetic coding choice.

The plot below uses Wilson 95% intervals and shows both conventions. At probabilities exactly zero or one, displayed error lengths are clipped at zero solely to avoid roundoff-sized negative Matplotlib error bars; stored intervals remain unchanged.

![Spatial covariance and adjacency sensitivity](/home/ubuntu/nickelate_research_audit/spatial_covariance_sensitivity.png)

### Tests and what they cover

Eight independent tests pass for the new spatial code: empty/full limits; a spanning horizontal path and a broken path; diagonal-only connectivity distinction; breadth-first-search agreement with a separately written union-find algorithm on 20 seeded masks under both adjacency conventions; covariance calibration; cropped-field normalization; Wilson interval bounds; and the white-noise active fraction relative to the analytic marginal. The tight-binding notebook audit separately passes six checks recorded in `audit_results.json`.

These tests do **not** establish lattice-size convergence, a thermodynamic-limit transition, accuracy of the missing original generator, or validity of the model as a description of nickelates. The covariance families were matched only by one scalar nearest-neighbour correlation, not by the entire covariance function. Only 64 × 64 maps were simulated in this new sensitivity run. No correlated-field Monte Carlo standard errors beyond the displayed Wilson intervals, model-averaged estimates, or physical proxy calibration are claimed.

## Tight-binding notebook audit and reproduced values

The slab uses `epsilon(k) = −2t(cos kx + cos ky)`, a finite open chain along layers, and a uniform in-plane momentum mesh. Each `K × K` point contributes to an `L × L` Hermitian matrix, which is diagonalized using `numpy.linalg.eigh`. Layer occupation includes spin degeneracy two. The notebook solves the chemical potential separately for every scenario to maintain target mean filling.

At eight layers, filling 1, `T/t = 0.04`, `t_perp/t = 0.20`, and interface onsite shift `+0.30t`, the reproduced 45 × 45 run gives interface `Δn_0 = −0.126066`, profile RMS `0.047690`, and `Σ_l Δn_l = −8.88 × 10⁻¹⁶`. This sum is a consequence of the imposed fixed-total-density condition. If the reference chemical potential is instead held fixed, mean filling becomes 0.982516 and `Σ_l Δn_l = −0.139869`.

Across the representative k-grid sweep, the interface response changes from −0.130943 at 25 × 25 to −0.125817 at 75 × 75. Between 45 × 45 and 75 × 75 it changes by about 0.20%; the profile RMS changes by about 0.19%. This is a local convergence check for that filling, thickness and perturbation, not a general convergence claim. The profile is one depleted interface layer compensated by near-uniform gains in the remaining seven layers. Its sign-blind log-linear fit, approximately 68 layers with `R² ≈ 0.068`, is not an identifiable screening length. The imposed potential does not feed back from density through electrostatics; therefore this notebook does not calculate screening.

## Stage plan and preregisterable decisive computational experiment

### Stage A — baseline reproduction

The current notebook baseline has been executed and compared with its recorded output. For the percolation manuscript, Stage A is **blocked** until the actual implementation, exact configuration, versioned seed generation, raw per-map results, and figure-generation code are obtained. A table in manuscript text is not a computational reproduction. Once those files are supplied, compare every table cell and figure point to a fresh run under the documented environment; retain differences as discrepancies rather than tuning to the target values.

### Stage B — controlled parameter sweep

First vary one factor at a time around the specified central point: spatial standard deviation, window width, target nearest-neighbour correlation, lattice size, and covariance family. Then use a factorial design over `σ ∈ {0.10, 0.13, 0.16}`, `w ∈ {0.06, 0.08, 0.10}`, correlation target `ρ1 ∈ {0, 0.2, 0.4, 0.6, 0.8}`, covariance family (at least Gaussian, exponential-like, and power-law-tailed), and `N ∈ {32, 64, 128, 256}`. Include four-neighbour as the primary rule and eight-neighbour only as an explicit alternative. Use a target active fraction control generated by adjusting `w` or threshold per map according to a preregistered method; do not quietly mix fixed-window and exactly fixed-fraction analyses. Separate exploratory from confirmatory seeds.

### Stage C — numerical and statistical verification

Use independent SeedSequence streams per condition; record them, failed-map counts and warnings. Verify field mean/variance, empirical full covariance, marginal active fraction, threshold calibration, and agreement of independent connectivity solvers. Declare numerical convergence when increasing map size changes the primary spanning estimate by no more than 0.02 absolute **and** does not change the direction of the preregistered effect across two successive size doublings, with Monte Carlo half-width no larger than 0.03. These are proposed tolerances, not evidence that the current run has met them. For confirmatory estimates, choose replicate count prospectively so Wilson half-width is at most 0.03 near the expected transition; use independent seed batches to check between-seed variability. Report both Monte Carlo error and variation over model choices; do not combine these into a single interval unless a defined hierarchical procedure is adopted.

### Stage D — alternatives and model comparison

Compare covariance families at matched marginal active fractions and matched measured nearest-neighbour correlations, then report unmatched full covariance functions. Compare fixed half-window against exact-per-map active fraction, four- versus eight-neighbour adjacency, and central crop versus periodic boundary conditions. For the tight-binding project, compare fixed density with fixed chemical potential, signed profile shapes, and (only after defining physical support) exponential, plateau, and damped-oscillation fits with residual checks and information criteria. Do not fit a length when no model is adequate.

### Stage E — independent replication

Run confirmatory realizations using a separately generated seed stream and an independent implementation. The union-find test currently cross-checks connectivity on binary masks; it does not independently reproduce the random-field spectral generator. A stronger replication should independently implement the covariance or use a direct covariance sampler on smaller grids and compare empirical covariance and spanning statistics. Archive code, config, seeds, raw per-map outcomes, and software environment.

### Prespecified result labels for the primary mathematical hypothesis

For a selected range of `ρ1` that brackets a partial-spanning regime, call the mathematical hypothesis supported only if the matched-active-fraction correlation contrast is at least 0.20 absolute, its confirmatory interval excludes zero, the sign replicates across independent seed batches, and it persists under at least two covariance families and the declared finite-size criterion. Call it weak/ambiguous if it appears only in one family, one grid size, or one seed batch. Call it unsupported if the interval includes effects both below and above the predeclared meaningful contrast after convergence. Call it contradicted in that parameter region if a well-powered, converged run shows a stable opposite-sign contrast. Call it a numerical artifact if it disappears under solver agreement, normalization checks, covariance validation, or grid convergence. These criteria are proposal for the next confirmatory run; they were not retroactively applied to claim the current exploratory run is conclusive.

## Material-level follow-up: a predictive question, not a percolation assertion

The strongest scientifically meaningful next question is whether spatially resolved reduction-state statistics improve prediction of superconducting response across **independent growth batches** beyond the batch mean and established covariates. Prior studies establish superconductivity in nickelate thin films and discuss sample variation, reduction, precursor stoichiometry and structural factors, but the examined sources do not establish that this threshold-percolation field is a nickelate mechanism [1] [13] [14]. Related cuprate and disordered-superconductor literature already treats inhomogeneous local superconductivity and percolative emergence, so generic novelty is not defensible [15] [16] [17].

A feasible future study, contingent on qualified thin-film synthesis and characterization facilities, would use a blocked, preregistered design with several independent growth batches and matched precursor films assigned to under-reduced, target, and extended-reduction endpoints. Measure cation composition; structural phase fraction, strain and out-of-plane lattice spacing; and one independently calibrated local reduction proxy with stated spatial resolution and uncertainty. Use local STEM-EELS or other spatial spectroscopy only where the required sensitivity and resolution are demonstrated; bulk-averaged XPS or diffraction alone cannot be relabelled a local oxygen map. Pair four-probe transport with an independent diamagnetic-screening measurement; include instrument-current, filtering, contact, and substrate controls. Independent growth batch, not pixels, devices or repeated anneals of one specimen, is the inferential unit.

Before outcomes are inspected, define the continuous transport and magnetic endpoints, detection thresholds, exclusions, primary spatial statistic, covariates, and minimum useful predictive gain. Compare a base model (mean proxy plus precursor composition, phase and strain) with a prespecified added spatial-statistics model using leave-one-growth-batch-out validation. Validate the spatial proxy against a reference method and test repeatability with blinded duplicate maps. A gain in prediction is not proof that classical percolation causes superconductivity; that stronger claim requires physically grounded connectivity/coupling and spatially resolved superconducting evidence.

## Current conclusion

The supplied notebook supports reproducible algebraic density responses in a small noninteracting slab, not an electrostatic screening length. The manuscript's reported spatial correlation contrast is plausible as a conditional threshold-model result but cannot be independently reproduced from supplied inputs. A separate, executed sensitivity experiment confirms that the result can depend strongly on adjacency convention and that two covariance spectra matched by nearest-neighbour correlation give similar but not identical four-neighbour estimates in the tested partial-spanning regime. The current evidence therefore justifies a revised research direction and a reproducible computational follow-up, not a material-specific superconductivity conclusion.

## References

[1]: https://doi.org/10.1038/s41586-019-1496-5 "Superconductivity in an infinite-layer nickelate"
[13]: https://doi.org/10.3389/fphy.2021.719534 "Impact of cation stoichiometry on the crystalline structure and superconductivity in nickelates"
[14]: https://doi.org/10.1002/advs.202309092 "Toward reliable synthesis of superconducting infinite layer nickelate thin films by topochemical reduction"
[15]: https://doi.org/10.1038/s41467-018-06707-y "Emergence of superconductivity in the cuprates via a universal percolation process"
[16]: https://doi.org/10.1038/s41598-020-60037-y "Film-thickness-driven superconductor to insulator transition in cuprate superconductors"
[17]: https://doi.org/10.1103/PhysRevB.78.014506 "Evidence of spatially inhomogeneous pairing on the insulating side of a disorder-tuned superconductor-insulator transition"
