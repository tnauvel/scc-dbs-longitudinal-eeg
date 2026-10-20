# scc-dbs-longitudinal-eeg

Longitudinal scalp EEG biomarkers of recovery during subcallosal cingulate deep brain stimulation (SCC-DBS) for treatment-resistant depression.

Stripped-down, self-contained MATLAB implementations of the core methods in:

> Nauvel T, et al. Scalp EEG brain states track recovery during subcallosal cingulate deep brain stimulation for treatment-resistant depression. In preparation, 2026.

Each folder holds one method as a flat function with plain array inputs, a README that states the method and what it showed in the paper, and a synthetic example that runs without patient data. The full pipeline (file handling, subject bookkeeping, hemisphere splits, figure generation, Python phase-amplitude coupling toolkit) is in a private repository.

## Study in one paragraph

Eleven participants with treatment-resistant depression were followed longitudinally with resting scalp EEG at eight visits: before implantation, at the post-implantation visit before chronic stimulation began, and at weeks 4, 8, 12, 16, 20 and 24
of stimulation. All recordings were made with stimulation off. Per-channel spectral clustering identified three brain states that changed over time; state occupancy tracked HDRS-17 within
 subjects and preceded symptom change in a cross-lagged model. Coherence networks showed an acute, reversible loss of global efficiency at the post-implantation visit. Baseline FA predicted 
baseline beta and gamma coherence but not clinical outcome, while theta-gamma phase-amplitude coupling tracked outcome and covaried with LFP gamma power recorded at the stimulation target, SCC.

## Methods

| Folder | Method | Language | Paper result |
|---|---|---|---|
| `spectral-states/` | PCA-denoised k-means on windowed spectra, labels ordered along the treatment axis | MATLAB | Three ordered states; occupancy trajectories (Fig. 1) |
| `temporal-precedence/` | Within-subject occupancy to HDRS correlation, Fisher-z pooling, cross-lagged panel model | MATLAB | EEG leads symptom change, beta 0.34 vs 0.12 (Fig. 1) |
| `coherence-networks/` | BCT graph metrics with proportional and MST thresholds; dwPLI control | MATLAB | Acute global-efficiency drop, recovered by week 24 (Fig. 2) |
| `structure-function/` | Spearman FA vs EEG family with FDR and leave-one-out range | MATLAB | FA vs EEG metrics, rho about -0.8 (Fig. 3) |
| `phase-amplitude-coupling/` | Scalp theta-gamma coupling (Tort modulation index, surrogate z) over time; within-subject repeated-measures correlation with HDRS and with SCC LFP gamma power | Python | PAC tracks outcome and covaries with SCC gamma (Fig. 4) |

## Requirements

MATLAB R2020b or later, Statistics and Machine Learning Toolbox. `coherence-networks` also needs the Brain Connectivity Toolbox. `phase-amplitude-coupling` is Python 3 (numpy, scipy, pingouin); see its own README. Spectral estimation (Chronux multitaper), preprocessing (EEGLAB) and the dwPLI computation (FieldTrip) are not reproduced here; the READMEs state the parameters used.

## Data

No patient data are included. Frequency bands used throughout: delta 1 to 4, theta 4 to 8, alpha 8 to 12, low beta 12 to 20, high beta 20 to 30, gamma 30 to 50 Hz.

## Authorship and AI assistance

The methods, the analysis decisions and the original implementation are the author's (2022 to 2026). The standalone functions in this repository were refactored from that pipeline with the assistance of Claude (Anthropic) in October 2026, then reviewed and validated by the author against the original pipeline output: on one subject the clustering function returns the published state labels exactly when started from the published centroids (37 of 37 channels), and the cross-lagged model reproduces the published coefficient (0.354 vs 0.34). Each function carries the same statement in its header.

## Citation

Cite the paper above. A Zenodo DOI will be added on publication.

## Licence

MIT. See `LICENSE`.

Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at Mount Sinai.

