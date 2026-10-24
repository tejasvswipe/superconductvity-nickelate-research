# Final research diagnosis and recommendation

## Verdict: E — a different research question offers a stronger contribution

The starting spatial-percolation hypothesis has a meaningful but conditional mathematical result: in a correlated random-field threshold model, changing spatial arrangement can change the probability of an edge-to-edge active component even when the mean proxy and active fraction are similar. The supplied manuscript reports a large contrast, but its generator, raw fields and per-realization outputs were not supplied, so that particular numerical contrast could not be reproduced independently. The newly implemented sensitivity model supports the plausibility of correlation-dependent geometry in a limited regime, while showing that adjacency convention can dominate the outcome. This is not evidence that nickelate films contain an optimal local reduction window or that an active geometric cluster carries supercurrent.

The generic physical idea of spatially inhomogeneous pairing followed by percolative emergence is established in cuprate and disordered-superconductor literature [15] [16] [17]. Therefore “correlation can change connectivity” alone is not a compelling novelty claim. The more valuable question is material-specific and predictive: **do spatial statistics of a calibrated local reduction-state proxy improve prediction of transport-plus-magnetic superconducting response across independent nickelate growth batches, beyond mean reduction proxy, precursor cation stoichiometry, phase fraction and strain?** This is a testable extension of the supplied protocol, not a result already established by it.

## What the supplied and reproduced evidence supports

The uploaded tight-binding notebook executes and reproduces its recorded outputs. In its one-orbital, noninteracting slab, fixed-density calculations enforce zero layer-summed density change. A fixed-chemical-potential control gives a nonzero net change. The Hamiltonian has no electrostatic feedback, so the notebook calculates response to an imposed one-body potential and cannot calculate a screened potential or a material screening length. Its representative log-linear fit is poor and sign-blind; the resulting nominal decay length is not identifiable.

The submitted percolation manuscript states that 64 × 64 correlated Gaussian maps at fixed mean and approximately equal active fraction yield different spanning probabilities for correlation lengths 1 and 3. Its pasted intervals agree with Wilson interval arithmetic. But absent source code and raw maps, arithmetic agreement is not independent simulation verification. The new, separate code compares two declared spectral covariances matched by nearest-neighbour correlation on a 64 × 64 crop. At target correlation 0.80, four-neighbour spanning was 61/300 (`p=0.203`, Wilson 95% CI 0.162–0.252) for the Gaussian spectrum and 45/300 (`p=0.150`, CI 0.114–0.195) for the Lorentzian-squared spectrum. The intervals overlap. Under eight-neighbour adjacency, spanning was already 297–300 out of 300 at zero target correlation and remained high throughout. This establishes sensitivity to modeling choices, not a universal threshold or nickelate physics.

The cited nickelate studies answer narrower questions. Li et al. (2021) relates nominal cation beam-flux ratio to structural indicators and superconducting transport in its specific Nd-based process; it does not establish an independent-batch reduction predictor. Gutiérrez-Llorente et al. (2024) follows successive reduction steps on one Pr-based piece, not independent growth-batch repeatability [13] [14]. The linked Google Doc remained inaccessible because no document text/export was available through the connected session. Any protocol details beyond those quoted in the pasted manuscript have not been assumed.

## Main weaknesses to repair

The original percolation code, full raw dataset, and exact field covariance are missing. The reported correlation-length label is not enough to reconstruct a covariance. The new run reveals that four- versus eight-neighbour adjacency changes the answer dramatically. Finite-size convergence and full-covariance matching remain untested. The percolation paper's local proxy, optimum and window are dimensionless assumptions without calibration to oxygen concentration, structure, local pairing, or transport.

The tight-binding slab uses fixed particle number as a computational constraint; zero-sum layer charge therefore cannot be interpreted as Coulomb screening. The temperature and hoppings are dimensionless model settings, not fitted nickelate parameters. No substrate, multiorbital hybridization, disorder, interactions, electrostatics, or magnetic response is present. A material-level interface-screening conclusion exceeds this model.

## Recommended next steps

First, obtain the inaccessible Google Doc as pasted text or export and request the missing original percolation code/configuration/raw outcomes. Second, complete the declared covariance × window × dispersion × lattice-size sweep with independent seed batches, empirical covariance checks, matched active fractions, an independent field-generation implementation, and the preregistered convergence and decision rules in the computational report. Do not spend on broader model complexity before this validation.

For a material study, use independent growth batches as the replication units; map a validated local proxy with calibration and declared spatial resolution; quantify composition, phase fraction and strain; and pair four-probe transport with an independent diamagnetic response. Compare held-out-batch predictive performance for established covariates alone versus the same set plus prespecified spatial statistics. Repeated devices, pixels, or treatments from one film are not independent growth replicates. Specialized synthesis, cryogenic transport, magnetic screening and nanoscale chemical mapping require suitable facilities and qualified supervision.

## What is delivered

The package includes the executed notebook reproduction and its audit, the new spatial sensitivity code/CSV/figure and tests, a derivation, literature review and source comparison, a ranked hypothesis matrix, the computational experiment report, and a revised manuscript that labels the original paper's reported results as unreproduced. It does not claim independent reproduction of the unavailable percolation output or material validation.

## References

[13]: https://doi.org/10.3389/fphy.2021.719534 "Impact of cation stoichiometry on the crystalline structure and superconductivity in nickelates"
[14]: https://doi.org/10.1002/advs.202309092 "Toward reliable synthesis of superconducting infinite layer nickelate thin films by topochemical reduction"
[15]: https://doi.org/10.1038/s41467-018-06707-y "Emergence of superconductivity in the cuprates via a universal percolation process"
[16]: https://doi.org/10.1038/s41598-020-60037-y "Film-thickness-driven superconductor to insulator transition in cuprate superconductors"
[17]: https://doi.org/10.1103/PhysRevB.78.014506 "Evidence of spatially inhomogeneous pairing on the insulating side of a disorder-tuned superconductor-insulator transition"
