# Spatial Correlation and Connectivity Rules Govern Spanning in a Thresholded Random Field

Tejasv Sahu and Aaradhaya Gupta  
Independent student researchers, India  
ORCID: 0009-0008-4301-3633 (T. Sahu) and 0009-0001-1160-6687 (A. Gupta)

## Abstract

Can two fields with the same mean and similar fractions of favorable sites differ in their probability of containing a connected path? We studied this question on 64 × 64 lattices using two random-field spectra, four target nearest-neighbor correlations, and 300 independent realizations per condition. We classified sites through a symmetric window around a fixed proxy value and tested each resulting map for a left-to-right crossing. Mean active fractions ranged from 0.46035 to 0.46181 across conditions. With edge-sharing connectivity, neither spectral family produced a crossing at zero correlation. At a target correlation of 0.80, 61 of 300 Gaussian-spectrum maps crossed (20.3%, 95% Wilson interval 16.2–25.2%), as did 45 of 300 Lorentzian-squared-spectrum maps (15.0%, 11.4–19.5%). Allowing diagonal contact changed the zero-correlation counts to 300 and 297, respectively. These results support a spatial effect on finite-grid crossing at comparable mean occupancy. Its magnitude depends on how connectivity is defined. Application to nickelate films requires a calibrated local proxy and independent transport and magnetic measurements.

**Index terms** — Site percolation, spatial correlation, random fields, nickelate thin films.

## I. Introduction and Hypothesis

A film average discards the positions of the sites that contribute to it. Favorable sites could form compact, isolated regions or a narrow path joining opposite edges. That distinction motivates a question about reduction heterogeneity in infinite-layer nickelate films. If a local material state is favorable only within a finite window, how much can the arrangement of sites inside that window influence connectivity?

Percolation has been used to interpret the emergence of superconducting response in cuprates [1], thickness-dependent transitions in cuprate films [2], and spatially inhomogeneous pairing in disordered indium-oxide films [3]. Nickelate studies provide experimental reasons to account for precursor cation stoichiometry and the reduction process [4], [5]. Whether the spatial statistics of reduction add predictive information beyond those measurements remains open.

The computational hypothesis is directional. At fixed proxy mean, dispersion, and window width, increasing the target nearest-neighbor correlation from 0 to 0.80 should increase the probability of an edge-sharing crossing in both spectral families while their mean active fractions remain comparable. The corresponding null hypothesis is equal endpoint crossing probabilities within each family. A secondary question concerns the same maps: how many additional crossings occur when corner contact is allowed? Because every edge-sharing path survives this change, the magnitude of the increase is the informative quantity.

## II. R&D Methodology

**Field generation.** We considered two non-negative power spectra, Gaussian and Lorentzian squared. Writing q for the spatial wavevector and ℓ for the spectral scale, the spectra are exp(−ℓ²|q|²/2) and (1 + ℓ²|q|²)⁻², respectively. Within each family, bisection selected ℓ to give a target nearest-neighbor correlation of 0, 0.35, 0.60, or 0.80 on the full periodic grid. Matching this statistic leaves the covariance at larger separations free to differ.

The spectral generator applied the square root of the selected power spectrum to Gaussian noise and inverse transformed the result on a 192 × 192 grid. We retained the central 64 × 64 crop, then recentered and rescaled it to zero mean and unit variance. Cropping separates the measured edges from the periodic boundary of the generating grid. The design comprised 2,400 fields in eight conditions, with seed streams derived from NumPy SeedSequence(53821).

**Local activity.** At site i, the proxy was xᵢ = 0.50 + 0.13zᵢ, where zᵢ is the standardized field. A site was active when |xᵢ − 0.50| ≤ 0.08. For a unit-normal variable, the corresponding reference probability is erf[0.08/(√2 × 0.13)] = 0.461699. Each crop has the specified mean and dispersion, but the number of its values falling within the window can fluctuate. The design therefore controls the proxy moments exactly and compares active fractions empirically.

**Outcome and uncertainty.** A realization counted as spanning if an active component touched both the left and right crop boundaries. Four-neighbor connectivity allowed shared edges. Eight-neighbor connectivity also allowed shared corners. We evaluated both rules on each mask and estimated crossing probability as the number of spanning masks divided by 300. Wilson 95% intervals describe Monte Carlo uncertainty for each condition. The unit of replication was the field. The source reports eight passed implementation checks, including comparison of breadth-first search with union-find on 20 seeded masks under both adjacency rules. These checks address algorithm behavior within the model.

## III. Results and Analysis

**Table I. Left-to-right crossings out of 300 maps per condition.** G denotes the Gaussian spectrum and L the Lorentzian-squared spectrum. Columns 4 and 8 give the number of neighbors allowed.

| Target correlation | G (4) | G (8) | L (4) | L (8) |
|---:|---:|---:|---:|---:|
| 0.00 | 0 | 300 | 0 | 297 |
| 0.35 | 0 | 295 | 0 | 296 |
| 0.60 | 0 | 288 | 1 | 289 |
| 0.80 | 61 | 298 | 45 | 292 |

The endpoint increase was 20.3 percentage points for the Gaussian spectrum and 15.0 points for the Lorentzian-squared spectrum. For comparison, the upper Wilson 95% bound for 0/300 is approximately 1.3%. Both families show the predicted direction, although even at 0.80 most fields lacked an edge-sharing crossing. Four correlation settings provide little resolution of the transition between rare and more frequent crossings.

Condition means for the active fraction spanned only 0.00146, or 0.146 percentage points. Thus the large endpoint change in crossings occurred with little change in average occupancy. The counts alone do not isolate an effect at identical occupancy in every realization. The 5.3-percentage-point difference between spectral families at 0.80 is descriptive. Overlap of their individual intervals does not test that difference.

Diagonal contact changed the zero-correlation result on the very same masks. It added 300 crossings in the Gaussian condition and 297 in the Lorentzian-squared condition, where neither had an edge-sharing crossing. Across all eight conditions, eight-neighbor counts ranged from 288 to 300. At 0.80, diagonal contact added 237 crossings in the Gaussian family and 247 in the Lorentzian-squared family. A rerun with the archived seeds reproduced the active-site count and both crossing outcomes for all 2,400 fields.

![Spanning probability under four- and eight-neighbor connectivity](spatial_covariance_sensitivity.png)

**Fig. 1.** Left-to-right crossing probability on 64 × 64 maps. Bars give Wilson 95% intervals for 300 independent fields per point.

**Earlier aggregate comparison.** A separate 64 × 64 study reported mean active fractions of 0.4613, 0.4622, and 0.4584 at correlation lengths of 1, 2, and 3 lattice units. Its reported crossing probabilities were 0.150, 0.980, and 0.988. These aggregate observations show the same broad direction, but a length in lattice units and a nearest-neighbor correlation specify different properties. Pooling the studies would require an explicit mapping between generators. The earlier 48 × 48 scan also omitted the mean-proxy column from its table header. **[DATA NEEDED]** The original generator, realization-level records, and verified row-to-mean mapping are required to reproduce and interpret that scan.

## IV. Discussion and Limitations

A narrow connection can decide whether a component reaches the far edge while changing the active fraction by only a few sites. The adjacency comparison demonstrates this directly: at zero target correlation, allowing movement across corners turned almost every mask into a crossing. Under edge-sharing connectivity, spatial correlation mattered most at the highest target tested. A claim about connectivity in an imaged film would therefore need a physical justification for which contacts can transmit current at the spatial resolution of the image.

Several questions remain open even at the geometric level. Crop standardization alters the finite-sample marginal distribution, and mean occupancy can conceal variation among maps. An informative control would shuffle each binary mask while preserving its exact count of active sites, then measure the paired change in crossing probability. To establish size dependence, repeat the design over successive lattice-size doublings with independent seeds. The proposed convergence criterion is a change below 0.02 in crossing probability over each of two successive doublings, with a Wilson half-width no larger than 0.03. The present sensitivity study used one size. **[DATA NEEDED]** Exact-occupancy comparisons, empirical covariance curves, and replicated size-series results are needed to evaluate these extensions.

**Material-level hypothesis.** Spatial variance and correlation of a calibrated local reduction-state proxy may improve prediction of joint transport and diamagnetic response across independently grown nickelate batches. A prospective study would compare a baseline model using mean reduction state, precursor composition, phase fraction, and strain with an otherwise identical model that adds prespecified spatial statistics. Entire growth batches should be held out together. Measurements from one film must stay in the same partition. Before collecting evaluation data, specify the proxy calibration, spatial resolution, transport and magnetic detection criteria, prediction metric, and minimum useful improvement. A gain that persists across held-out batches would support predictive value. Failure to reach the specified improvement with adequate precision would count against the hypothesis.

**[DATA NEEDED]** No calibrated reduction maps, paired transport and magnetic observations, independent batch counts, or held-out prediction errors are supplied. The binary model includes neither local superconducting response nor coupling, current, or magnetic screening. Its crossing result cannot establish a superconducting pathway in a material.

## V. Conclusion

At similar mean active fractions, increasing the target nearest-neighbor correlation from 0 to 0.80 raised edge-sharing crossing frequency from zero observed crossings to 61/300 and 45/300 in the two spectral families. Allowing corner contact made crossings common even at zero correlation. The calculations support the geometric hypothesis within the tested finite-grid model. A direct nickelate test would ask whether measured spatial structure improves prediction of both transport and magnetic response in independent growth batches.

## Data and Code Availability

The archived implementation, configuration, replicate records and aggregate results used for Table I are available at commit 267535fd908d3e235b5c46d0fefe744221cf16e7. Repository: https://github.com/tejasvswipe/superconductvity-nickelate-research

## References

[1] D. Pelc et al., “Emergence of superconductivity in the cuprates via a universal percolation process,” Nature Communications, vol. 9, art. 4327, 2018. DOI: 10.1038/s41467-018-06707-y.

[2] H.-B. Jang et al., “Film-thickness-driven superconductor to insulator transition in cuprate superconductors,” Scientific Reports, vol. 10, art. 3236, 2020. DOI: 10.1038/s41598-020-60037-y.

[3] K. H. S. B. Tan, K. A. Parendo, and A. M. Goldman, “Evidence of spatially inhomogeneous pairing on the insulating side of a disorder-tuned superconductor-insulator transition,” Physical Review B, vol. 78, art. 014506, 2008. DOI: 10.1103/PhysRevB.78.014506.

[4] “Impact of cation stoichiometry on the crystalline structure and superconductivity in nickelates,” Frontiers in Physics, vol. 9, art. 719534, 2021. DOI: 10.3389/fphy.2021.719534.

[5] A. Gutiérrez-Llorente et al., “Toward reliable synthesis of superconducting infinite layer nickelate thin films by topochemical reduction,” Advanced Science, vol. 11, art. 2309092, 2024. DOI: 10.1002/advs.202309092.

