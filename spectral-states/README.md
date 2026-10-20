# scc-dbs-eeg-states

Spectral brain-state clustering of resting scalp EEG during subcallosal cingulate deep brain stimulation (SCC-DBS) for treatment-resistant depression.

This repository holds a stripped-down, self-contained version of the clustering method used in:

> Nauvel T, et al. Scalp EEG brain states track recovery during subcallosal cingulate deep brain stimulation for treatment-resistant depression. In preparation, 2026.

It is written to be read and reused. The full analysis pipeline, which handles file discovery, subject bookkeeping, hemisphere splits and clinical correlations, lives in a private repository.

## What the method does

Resting EEG is divided into 6 s non-overlapping windows. For each window and channel a multitaper power spectrum (1 to 50 Hz) is computed on two 3 s segments and averaged, then log-normalised across frequency and expressed in dB. Each channel is clustered separately: the windows from all visits are pooled, denoised by projecting onto the principal components that explain 95% of the variance, and partitioned with k-means. Clusters are then ordered along the treatment axis. State 1 is the cluster most occupied at baseline, state 3 the cluster most occupied at the final visit, and state 2 is the remainder. Occupancy per visit is the fraction of windows in each state.

Three states were used in the paper. Two states separate baseline from stimulation but cannot resolve a trajectory; three is the minimum partition that gives an origin, an intermediate and an endpoint. See the paper's Methods for the silhouette and BIC analysis behind that choice.

## Files

| File | Purpose |
|---|---|
| `spectral_state_clustering.m` | The method. One channel in, ordered labels and occupancy out. |
| `example_usage.m` | Runs the function on synthetic spectra so the repository works without patient data. |
| `LICENSE` | MIT. |

## Requirements

MATLAB R2020b or later with the Statistics and Machine Learning Toolbox (`kmeans`, `silhouette`). No other dependencies. Spectral estimation is not included; any multitaper or Welch estimate in dB will do as input.

## Usage

```matlab
% S     [nFreq x nWin]  power spectra in dB, one channel, all visits concatenated
% visit [nWin x 1]      integer visit index of each window (1 = baseline)
out = spectral_state_clustering(S, visit, 3);

out.labels      % [nWin x 1]  state of each window, 1 = baseline-dominant
out.occupancy   % [nVisit x 3] fraction of windows per visit in each state
out.silhouette  % mean silhouette of the partition
out.nComp       % principal components retained
out.centroids   % [3 x nFreq] mean spectrum of each state

out = spectral_state_clustering(S, visit, 3, 'VarKeep', 1, 'Plot', true);  % no denoising, draw occupancy
```

Run `example_usage.m` for a complete worked example.

## Data

No patient data are included. The paper's cohort is 11 participants with treatment-resistant depression implanted with SCC-DBS, recorded at eight visits (pre-implantation, stimulation turn-on, and weeks 4 to 24 of chronic stimulation), 37-channel Laplacian-referenced EEG, 1000 Hz, eyes closed, stimulation off during recording. Frequency bands used in the paper: delta 1 to 4, theta 4 to 8, alpha 8 to 12, low beta 12 to 20, high beta 20 to 30, gamma 30 to 50 Hz.

## Related repositories

- `scc-dbs-eeg-networks`: coherence network metrics with proportional and minimum-spanning-tree thresholding (coming).
- `scc-dbs-eeg-precedence`: within-subject occupancy to symptom correlation and cross-lagged precedence test (coming).

## Citation

If you use this code, cite the paper above. A Zenodo DOI will be added on publication.

## Contact

Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at Mount Sinai.

