%% example_usage.m
% Worked example for spectral_state_clustering on synthetic data.
% Builds three spectral states whose occupancy shifts across eight visits,
% runs the clustering, and draws the recovered occupancy.

rng(1);
f       = 1:0.25:50;                 % 197 frequency bins
nF      = numel(f);
nVisit  = 8;
nWinPer = 40;                        % windows per visit

% Three template spectra in dB: broadband 1/f plus a state-specific bump
base    = -10 * log10(f);
tmpl    = [base + 6 * exp(-((f - 6)  / 1.5).^2);     % state 1: theta-heavy ("sick")
           base + 4 * exp(-((f - 10) / 2.0).^2);     % state 2: alpha      ("intermediate")
           base + 5 * exp(-((f - 22) / 4.0).^2)];    % state 3: beta       ("well")

% Occupancy schedule: state 1 dominates early, state 3 late
p = [0.80 0.15 0.05;
     0.60 0.30 0.10;
     0.45 0.40 0.15;
     0.30 0.45 0.25;
     0.20 0.40 0.40;
     0.15 0.30 0.55;
     0.10 0.25 0.65;
     0.05 0.20 0.75];

S = []; visit = []; truth = [];
for v = 1:nVisit
    st = randsample(3, nWinPer, true, p(v, :));
    for w = 1:nWinPer
        S(:, end + 1) = tmpl(st(w), :)' + 1.5 * randn(nF, 1);  %#ok<SAGROW>
    end
    visit = [visit; v * ones(nWinPer, 1)];    %#ok<AGROW>
    truth = [truth; st];                      %#ok<AGROW>
end
S = S - mean(S, 1);                           % log-normalise across frequency, as in the paper

%% Cluster and compare with the generating states
out = spectral_state_clustering(S, visit, 3, 'Plot', true);

fprintf('Retained %d principal components, silhouette %.2f\n', out.nComp, out.silhouette);
fprintf('Label agreement with generating states: %.1f%%\n', 100 * mean(out.labels == truth));
disp('Recovered occupancy (rows = visits, cols = states):');
disp(round(100 * out.occupancy, 1));

