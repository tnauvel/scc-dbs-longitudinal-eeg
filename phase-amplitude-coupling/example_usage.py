"""
example_usage.py - synthetic demonstration of pac_tort_mi and rmcorr.

Part 1 builds 60 three-second epochs of a 6 Hz rhythm whose 30-50 Hz
activity is amplitude-modulated by the 6 Hz phase (coupled channel) and a
second channel with the same spectra but no coupling, then runs pac_channel.
Part 2 builds eight visits for ten subjects in which a session PAC z-score
covaries within subject with a second measure, and runs rm_corr.
"""
import numpy as np
from pac_tort_mi import pac_channel, pac_session
from rmcorr import rm_corr

rng = np.random.default_rng(0)
fs, epoch_s, n_ep = 1000, 3.0, 60
t = np.arange(int(fs * epoch_s)) / fs


def epoch(coupled):
    theta = np.sin(2 * np.pi * 6 * t + rng.uniform(0, 2 * np.pi))
    gamma = np.sum([np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
                    for f in rng.uniform(30, 50, 8)], axis=0)
    env = 1 + 0.6 * theta if coupled else 1.0
    return 3 * theta + 0.3 * env * gamma + 0.8 * rng.standard_normal(t.size)


coupled = np.array([epoch(True) for _ in range(n_ep)])
flat = np.array([epoch(False) for _ in range(n_ep)])
for name, ep in [('coupled channel', coupled), ('uncoupled channel', flat)]:
    r = pac_channel(ep, fs, (4, 8), (30, 50), n_surr=200)
    print(f'{name:18s} MI {r["mi"]:.4f}  z {r["z"]:6.2f}  p {r["p"]:.4f}')

# Part 2: within-subject covariation
subj, pac_z, other = [], [], []
for s in range(10):
    base = rng.uniform(1, 4)
    x = base + rng.standard_normal(8)
    y = 0.4 * x + rng.uniform(-2, 2) + 0.5 * rng.standard_normal(8)
    subj += [f'S{s:02d}'] * 8
    pac_z += list(x)
    other += list(y)
res = rm_corr(pac_z, other, subj, standardize=True)
print(f'rm_corr r = {res["r"]:.3f}, p = {res["p"]:.4f}, dof = {res["dof"]}, '
      f'{res["n_subjects"]} subjects, {res["n_obs"]} observations')

