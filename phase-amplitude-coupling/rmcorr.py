"""
rmcorr.py - repeated-measures correlation (Bakdash & Marusich 2017).

METHOD
------
The common within-subject association between two repeatedly measured
variables, estimated by ANCOVA with subject as a factor: y is regressed on
x and subject indicator variables, and r_rm is the signed square root of the
proportion of the subject-adjusted variance explained by x. Degrees of
freedom are N - k - 1 for N observations and k subjects. Each subject's
series can be standardised first ('standardize=True'), which the paper used
to satisfy the common-slope assumption across subjects.

REFERENCE
---------
Bakdash JZ, Marusich LR (2017) Front Psychol 8:456

AUTHORSHIP
----------
Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at
Mount Sinai. Original analysis used pingouin.rm_corr; this dependency-free
implementation was written with the assistance of Claude (Anthropic),
October 2026, and checked against pingouin on the same inputs. Licence: MIT.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def rm_corr(x, y, subject, standardize=False):
    """
    x, y     : 1-D arrays of paired observations
    subject  : 1-D array of subject labels, same length
    Returns dict(r, p, dof, n_subjects, n_obs).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    subject = np.asarray(subject)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, subject = x[ok], y[ok], subject[ok]
    subs = np.unique(subject)
    if standardize:
        for s in subs:
            m = subject == s
            if m.sum() > 1:
                x[m] = (x[m] - x[m].mean()) / x[m].std(ddof=1)
                y[m] = (y[m] - y[m].mean()) / y[m].std(ddof=1)
    n, k = x.size, subs.size
    dof = n - k - 1
    if dof < 1:
        return dict(r=np.nan, p=np.nan, dof=dof, n_subjects=k, n_obs=n)
    # subject dummies (full set, no intercept) + x
    D = (subject[:, None] == subs[None, :]).astype(float)
    X_full = np.column_stack([D, x])
    X_sub = D
    b_full, *_ = np.linalg.lstsq(X_full, y, rcond=None)
    b_sub, *_ = np.linalg.lstsq(X_sub, y, rcond=None)
    ss_err = np.sum((y - X_full @ b_full) ** 2)
    ss_sub = np.sum((y - X_sub @ b_sub) ** 2)
    ss_x = ss_sub - ss_err                       # variance explained by x after subject
    r = np.sign(b_full[-1]) * np.sqrt(max(ss_x, 0) / (ss_x + ss_err))
    t = r * np.sqrt(dof / max(1 - r ** 2, 1e-300))
    p = 2 * stats.t.sf(abs(t), dof)
    return dict(r=float(r), p=float(p), dof=int(dof), n_subjects=int(k), n_obs=int(n))

