# structure-function

Baseline white-matter microstructure against a family of EEG measures, with FDR control.

## What the method does

A per-subject structural value (in the paper, cingulum bundle fractional anisotropy from diffusion imaging) is rank-correlated with each member of a pre-specified EEG family (in the paper, baseline whole-head coherence in six bands). P-values are corrected within the family by Benjamini-Hochberg. A leave-one-out range of rho is reported for every test so that a correlation carried by a single subject is visible.

In the paper, higher cingulum FA was associated with lower baseline coherence in the beta and gamma bands (rho about -0.8, n = 11, FDR-corrected), while delta, theta and alpha were not significant. The association held for coherence but not for band power, and was global rather than specific to a tract or a region.

## Files

| File | Purpose |
|---|---|
| `structure_function_correlation.m` | One structural vector and one EEG matrix in; rho, p, FDR q and leave-one-out range per measure out. |
| `example_usage.m` | Synthetic six-band family with two true associations. |

## Requirements

MATLAB R2020b or later with the Statistics and Machine Learning Toolbox (`corr`).

## Usage

```matlab
out = structure_function_correlation(fa, eeg_by_band, band_names, 'Plot', true);
```

