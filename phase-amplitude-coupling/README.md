# phase-amplitude-coupling

Scalp theta-gamma phase-amplitude coupling across treatment, and its within-subject relationship to symptoms and to activity at the stimulation target.

## What the method does

For each channel and session, epochs are bandpass filtered independently (zero-phase FIR), Hilbert transformed and trimmed of their filter transient, then pooled. Coupling is the Tort et al. (2010) modulation index between theta phase (4 to 8 Hz) and gamma amplitude (30 to 50 Hz) over 18 phase bins, expressed as a z-score against 200 surrogates in which each epoch's gamma envelope is re-paired with the theta phase of a different epoch. The median z across channels summarises a session. The amplitude filter is required to be at least twice as wide as the highest phase frequency so that coupling sidebands survive filtering.

Session z-scores are then related to a second repeatedly measured variable (HDRS-17, or SCC LFP gamma power recorded by the implanted device in the matching week) with the repeated-measures correlation of Bakdash and Marusich (2017), on within-subject standardised series.

In the paper, scalp theta-gamma coupling tracked 24-week outcome between subjects and covaried within subject with SCC gamma power (mean rho +0.41, 7 of 9 subjects positive, p = 0.0065), but not with scalp gamma power or SCC theta, suggesting the scalp measure indexes activity at the target rather than cortical amplitude.

## Files

| File | Purpose |
|---|---|
| `pac_tort_mi.py` | Filtering, Tort MI, epoch-shuffle surrogates, bandwidth guard; `pac_channel` and `pac_session` drivers. |
| `rmcorr.py` | Repeated-measures correlation, dependency-free, checked against `pingouin.rm_corr`. |
| `example_usage.py` | Synthetic coupled vs uncoupled channel, and a synthetic within-subject covariation. |

## Requirements

Python 3.9 or later, numpy, scipy. No other dependencies.

## Usage

```python
from pac_tort_mi import pac_channel, pac_session
from rmcorr import rm_corr

r = pac_channel(epochs, fs, phase_band=(4, 8), amp_band=(30, 50), n_surr=200)   # epochs: (n_epochs, n_samples)
z, median_z = pac_session(data, fs)                                            # data: (n_channels, n_epochs, n_samples)
res = rm_corr(pac_z, lfp_gamma, subject, standardize=True)                      # paired observations across visits
```

Run `python example_usage.py` for a complete worked example.

