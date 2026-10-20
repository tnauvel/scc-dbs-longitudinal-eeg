%% example_usage.m
% Synthetic demonstration of structure_function_correlation: one structural
% measure, six band-wise EEG measures, a true negative association in two
% of them, FDR within the six-band family.

rng(2);
nSubj = 11;
fa    = 0.45 + 0.08 * rand(nSubj, 1);
bands = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta', 'gamma'};
eeg   = 0.25 + 0.05 * randn(nSubj, 6);
eeg(:, 5) = 0.90 - 1.40 * fa + 0.015 * randn(nSubj, 1);   % high beta tracks FA
eeg(:, 6) = 0.80 - 1.20 * fa + 0.015 * randn(nSubj, 1);   % gamma tracks FA

out = structure_function_correlation(fa, eeg, bands, 'Plot', true);
disp(table(out.names, out.rho, out.p, out.q, out.rho_loo_min, out.rho_loo_max, ...
     'VariableNames', {'measure', 'rho', 'p', 'q_FDR', 'loo_min', 'loo_max'}));

