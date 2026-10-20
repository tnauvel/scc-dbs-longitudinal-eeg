# temporal-precedence

Does an EEG measure track symptom change within subjects, and does it lead it?

## What the method does

Two questions, answered on the same per-subject by per-visit arrays. First, the within-subject Spearman correlation between the EEG measure (in the paper, state-1 occupancy) and the concurrent HDRS-17 score, pooled across subjects by Fisher's z with a circular-shift permutation test. Second, a cross-lagged panel model on within-subject standardised series: score at visit T+1 regressed on score at T and EEG at T, and the reverse. If the EEG-to-score coefficient is larger and significant while the reverse is not, the EEG measure leads symptom change.

In the paper the pooled correlation between state-1 occupancy and HDRS-17 was positive and significant, and the cross-lagged model gave an EEG-to-score path of 0.34 (p = 0.003) against a reverse path of 0.12 (n.s.).

## Files

| File | Purpose |
|---|---|
| `occupancy_behavior_crosslag.m` | Two matrices in (EEG measure, symptom score), correlations and cross-lagged coefficients out. |
| `example_usage.m` | Synthetic cohort in which the EEG measure leads the score by one visit. |

## Requirements

MATLAB R2020b or later with the Statistics and Machine Learning Toolbox (`corr`, `fitlm`).

## Usage

```matlab
% occ, score: [nSubj x nVisit], NaN where a visit is missing
out = occupancy_behavior_crosslag(occ, score, 'NPerm', 10000, 'Plot', true);
```

