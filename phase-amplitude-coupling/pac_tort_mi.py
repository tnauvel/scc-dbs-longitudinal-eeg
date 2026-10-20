"""
pac_tort_mi.py - phase-amplitude coupling on epoched scalp EEG.

METHOD
------
For one channel, each epoch is bandpass filtered independently (zero-phase
FIR, 3 cycles of the lower cutoff), Hilbert transformed, and the filter
transient is trimmed from both ends before the phase and amplitude samples
are pooled across epochs. Coupling is the Tort et al. (2010) Kullback-Leibler
modulation index over 18 phase bins, normalised to [0, 1]. Significance is
assessed against a surrogate distribution in which the amplitude envelope of
every epoch is re-paired with the phase of a different epoch (a derangement),
which breaks the phase-amplitude correspondence while leaving each signal's
own spectrum and waveform intact. The channel z-score against that null is
the quantity carried forward; the median z across channels summarises a
session.

Two guards are built in. The amplitude filter must be at least twice as wide
as the highest phase frequency, otherwise the modulation sidebands are
removed before they can be measured (Aru et al. 2015). Epochs are never
concatenated before filtering, so no phase continuity is manufactured across
epoch boundaries.

REFERENCES
----------
Tort ABL et al. (2010) J Neurophysiol 104:1195-1210
Canolty RT et al. (2006) Science 313:1626-1628
Aru J et al. (2015) Curr Opin Neurobiol 31:51-61
Donoghue T, Dominguez J, Voytek B (2022) eLife 11:e77548

AUTHORSHIP
----------
Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at
Mount Sinai. Method and original implementation (pac_toolkit.py): T. Nauvel,
2026. This release version was reduced to the computational core with the
assistance of Claude (Anthropic), October 2026, and checked by the author to
return identical modulation indices and surrogate z-scores to the original
on the same inputs. Licence: MIT.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import filtfilt, firwin, hilbert


# ---------------------------------------------------------------- guards
def check_bandwidth(phase_band, amp_band):
    """(ok, max_valid_phase_hz, message). Amplitude filter must be >= 2 x
    the highest phase frequency and must not touch the phase band."""
    p_lo, p_hi = float(phase_band[0]), float(phase_band[1])
    a_lo, a_hi = float(amp_band[0]), float(amp_band[1])
    max_phase = (a_hi - a_lo) / 2.0
    if a_lo <= p_hi:
        return False, max_phase, 'phase and amplitude bands overlap or touch'
    if p_hi > max_phase:
        return False, max_phase, (f'amplitude band is {a_hi - a_lo:g} Hz wide, '
                                  f'phase estimable only up to {max_phase:g} Hz')
    return True, max_phase, 'ok'


# ------------------------------------------------------------- filtering
def _fir_taps(fs, lo, n_pts, cycles=3.0):
    n = int(round(cycles * fs / float(lo)))
    n = min(n, (n_pts // 3) - 1)          # filtfilt padding limit
    n = max(n, 15)
    return n + 1 if n % 2 == 0 else n


def bandpass_epoch(x, fs, band, cycles=3.0):
    """Zero-phase linear-phase FIR bandpass of one epoch."""
    lo, hi = float(band[0]), float(band[1])
    nyq = fs / 2.0
    if hi >= nyq:
        raise ValueError(f'band upper edge {hi} Hz >= Nyquist {nyq} Hz')
    taps = firwin(_fir_taps(fs, lo, len(x), cycles), [lo / nyq, hi / nyq],
                  pass_zero=False)
    return filtfilt(taps, [1.0], x)


def auto_edge_s(phase_band, epoch_pts, fs, cycles=3.0, max_frac=0.30):
    """Seconds trimmed from each end: `cycles` periods of the lowest phase
    frequency, capped at max_frac of the epoch."""
    return min(cycles / float(phase_band[0]), max_frac * epoch_pts / fs)


# ------------------------------------------------- phase and amplitude
def epoch_phase_amp(epochs, fs, phase_band, amp_band, edge_s=None):
    """epochs: (n_epochs, epoch_pts) for one channel.
    Returns per-epoch lists of phase (rad) and amplitude, transient trimmed."""
    epochs = np.atleast_2d(np.asarray(epochs, dtype=float))
    n_ep, n_pts = epochs.shape
    if edge_s is None:
        edge_s = auto_edge_s(phase_band, n_pts, fs)
    k = int(round(edge_s * fs))
    if 2 * k >= n_pts:
        raise ValueError('edge trim exceeds epoch length')
    phases, amps = [], []
    for ep in epochs:
        phases.append(np.angle(hilbert(bandpass_epoch(ep, fs, phase_band)))[k:n_pts - k])
        amps.append(np.abs(hilbert(bandpass_epoch(ep, fs, amp_band)))[k:n_pts - k])
    return phases, amps


# ------------------------------------------------------ modulation index
def tort_mi(phase, amp, n_bins=18):
    """Tort (2010) KL modulation index in [0, 1]; NaN if a phase bin is empty."""
    phase = np.asarray(phase)
    amp = np.asarray(amp)
    edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    idx = np.clip(np.digitize(phase, edges) - 1, 0, n_bins - 1)
    mean_amp = np.zeros(n_bins)
    for b in range(n_bins):
        sel = idx == b
        if not sel.any():
            return np.nan
        mean_amp[b] = amp[sel].mean()
    total = mean_amp.sum()
    if total <= 0:
        return np.nan
    p = np.clip(mean_amp / total, 1e-300, None)
    h = -(p * np.log(p)).sum()
    return float((np.log(n_bins) - h) / np.log(n_bins))


# ------------------------------------------------------------ surrogates
def _derangement(n, rng, tries=100):
    if n < 2:
        return np.arange(n)
    for _ in range(tries):
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p
    return np.roll(rng.permutation(n), 1)


def epoch_shuffle_null(phases, amps, n_surr=200, n_bins=18, seed=42):
    """Null MI distribution from re-pairing each epoch's phase with another
    epoch's amplitude envelope."""
    rng = np.random.default_rng(seed)
    ph_cat = np.concatenate(phases)
    null = np.empty(n_surr)
    for s in range(n_surr):
        perm = _derangement(len(amps), rng)
        null[s] = tort_mi(ph_cat, np.concatenate([amps[i] for i in perm]), n_bins)
    return null


def surrogate_stats(mi, null):
    """(z, p): z against the null; one-sided exceedance p with +1 correction."""
    null = np.asarray(null)
    null = null[np.isfinite(null)]
    if not np.isfinite(mi) or null.size == 0:
        return np.nan, np.nan
    sd = null.std(ddof=1)
    z = (mi - null.mean()) / sd if sd > 0 else np.nan
    p = (np.sum(null >= mi) + 1.0) / (null.size + 1.0)
    return float(z), float(p)


# ------------------------------------------------------- one-call driver
def pac_channel(epochs, fs, phase_band=(4, 8), amp_band=(30, 50),
                n_surr=200, n_bins=18, seed=42):
    """Full pipeline for one channel. Returns dict(mi, z, p, null_mean, null_sd)."""
    ok, _, msg = check_bandwidth(phase_band, amp_band)
    if not ok:
        raise ValueError(msg)
    ph, am = epoch_phase_amp(epochs, fs, phase_band, amp_band)
    mi = tort_mi(np.concatenate(ph), np.concatenate(am), n_bins)
    null = epoch_shuffle_null(ph, am, n_surr, n_bins, seed)
    z, p = surrogate_stats(mi, null)
    return dict(mi=mi, z=z, p=p, null_mean=float(np.nanmean(null)),
                null_sd=float(np.nanstd(null, ddof=1)))


def pac_session(data, fs, phase_band=(4, 8), amp_band=(30, 50), **kw):
    """data: (n_channels, n_epochs, epoch_pts). Returns (z per channel,
    median z) - the session summary used in the paper."""
    z = np.array([pac_channel(data[c], fs, phase_band, amp_band, **kw)['z']
                  for c in range(data.shape[0])])
    return z, float(np.nanmedian(z))

