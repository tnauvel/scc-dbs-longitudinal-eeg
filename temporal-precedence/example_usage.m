%% example_usage.m
% Synthetic demonstration of occupancy_behavior_crosslag: 11 subjects, 8
% visits, an EEG measure that falls one visit ahead of the symptom score.

rng(4);
nSubj = 11; nVisit = 8;
t   = (0:nVisit-1)';
occ = zeros(nSubj, nVisit); score = zeros(nSubj, nVisit);
for s = 1:nSubj
    rate = 0.4 + 0.4 * rand;                                  % subject-specific recovery speed
    e = 0.8 * exp(-rate * t) + 0.05 * randn(nVisit, 1);       % EEG measure (e.g. state-1 occupancy)
    h = 24 * exp(-rate * max(t - 1, 0)) + 2 * randn(nVisit, 1); % HDRS lags EEG by one visit
    occ(s, :) = e; score(s, :) = h;
end

out = occupancy_behavior_crosslag(occ, score, 'NPerm', 2000, 'Plot', true);
fprintf('pooled r = %.2f [%.2f, %.2f], permutation p = %.3f\n', out.r_pooled, out.ci_pooled, out.p_perm);
fprintf('EEG -> score beta = %.2f (p = %.3f); score -> EEG beta = %.2f (p = %.3f)\n', ...
    out.crosslag.beta_eeg2score, out.crosslag.p_eeg2score, out.crosslag.beta_score2eeg, out.crosslag.p_score2eeg);

