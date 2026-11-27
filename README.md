# Spatial correlation and connectivity in a thresholded random field

Tejasv Sahu and Aaradhya Gupta study how the arrangement of active sites changes the probability of a left-to-right crossing on a finite grid. A 2,400-field sensitivity study tests two spectra and two connectivity rules. A further 7,200 original–shuffle pairs test spatial arrangement at identical occupancy within each pair. Edge-sharing crossings became less frequent as the grid grew, and neither spectral family met the proposed convergence criterion.

[Read the paper](revised_manuscript.md) · [Download the PDF](revised_manuscript.pdf) · [LaTeX source](revised_manuscript.tex)

## Data behind the paper

| Analysis | Implementation | Results |
|---|---|---|
| Fixed-window sensitivity | [Generator](spatial_covariance_sensitivity.py) and [configuration](spatial_config.json) | [Condition summaries](spatial_covariance_sensitivity.csv) and [individual maps](spatial_covariance_sensitivity_replicates.csv) |
| Finite-size and shuffle controls | [Experiment](finite_size_control.py) | [Batch and pooled summaries](finite_size_results/summary.csv) and [paired outcomes](finite_size_results/maps.csv) |
| Earlier correlation-length comparison | [Original archived manuscript and aggregate tables](spatial_percolation_manuscript_and_data.txt) | The original generator and individual maps were unavailable. |

## Reproduce the computations

Run the following from the repository root in a Python environment. The finite-size experiment used NumPy 2.3.5. The archived sensitivity run records its environment in [run metadata](run_metadata.json). The simulation scripts write their result CSVs into the checkout.

```sh
python -m pip install numpy==2.3.5 matplotlib
python -m unittest -v test_spatial_percolation.py
python spatial_covariance_sensitivity.py
python plot_spatial_sensitivity.py
python finite_size_control.py
```

Compile the PDF with `tectonic revised_manuscript.tex` from the repository root so the two figure paths resolve.

## Supporting material

The [model derivation](model_derivation.md), [literature review](literature_review.md), and [separate tight-binding notebook](nickelate_interface_model_colab.ipynb) provide project background. The paper evaluates geometric crossing in simulated fields. Its proposed experiment on nickelate films requires linked measurements from independently grown batches.
