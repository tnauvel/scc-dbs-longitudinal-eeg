function out = spectral_state_clustering(S, visit, k, varargin)
% SPECTRAL_STATE_CLUSTERING  Discretise a channel's windowed power spectra into
% k ordered brain states (sick -> intermediate -> well).
%
% METHOD
%   1. Log-power spectra are mean-centred and projected onto the leading
%      principal components that together explain >= VarKeep of the variance
%      (typically 3), then projected back into frequency space. This removes
%      high-dimensional noise so that distances are meaningful.
%   2. k-means with multiple restarts partitions the denoised spectra.
%   3. Labels are ordered along the treatment axis: state 1 is the cluster
%      most occupied at the baseline visit, state k the cluster most occupied
%      at the endpoint visit, remaining clusters fill the middle by baseline
%      occupancy (descending).
%   4. Occupancy per visit and mean silhouette are returned.
%
% REFERENCES
%   Alagapan S et al. (2023) Nature 622:130-138 (sick vs stable-response states)
%   Kaufman L, Rousseeuw PJ (1990) Finding Groups in Data (silhouette)
%   Nauvel T et al., in preparation (SCC-DBS scalp EEG states)
%
% INPUTS
%   S        [nFreq x nWin] spectra in dB (10*log10 power), one channel
%   visit    [nWin x 1]     integer visit index for each window
%   k        number of states (paper uses 3)
%   Name-value:
%     'VarKeep'   variance fraction retained by PCA (default 0.95; 1 = none)
%     'Baseline'  visit index treated as baseline (default min(visit))
%     'Endpoint'  visit index treated as endpoint (default max(visit))
%     'Replicates' k-means restarts (default 20)
%     'Plot'      true to draw the example occupancy figure (default false)
%     'Start'     [k x nFreq] initial centroids in the input space. Used for
%                 exact reproduction of a published partition: k-means
%                 started at that partition's centroids returns it unchanged.
%
% OUTPUT struct
%   .labels     [nWin x 1] ordered state label (1..k)
%   .occupancy  [nVisit x k] fraction of windows per visit in each state
%   .silhouette mean silhouette of the ordered partition
%   .nComp      number of principal components retained
%   .centroids  [k x nFreq] state spectra in dB
%
% Requires Statistics and Machine Learning Toolbox (kmeans, silhouette).
%
% AUTHORSHIP
%   Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at
%   Mount Sinai. Method and original implementation: T. Nauvel, 2022-2026
%   (full analysis pipeline in a private repository, see the README).
%   This release version was refactored into a standalone function with the
%   assistance of Claude (Anthropic), October 2026, then reviewed and
%   validated by the author against the original pipeline:
%   DBS912, 37 channels: started at the published centroids the function returns the published labels exactly (37/37); random restarts agree 96%.
%   Licence: MIT.

p = inputParser;
p.addParameter('VarKeep', 0.95);
p.addParameter('Baseline', min(visit));
p.addParameter('Endpoint', max(visit));
p.addParameter('Replicates', 20);
p.addParameter('Plot', false);
p.addParameter('Start', []);
p.parse(varargin{:});
o = p.Results;

visit = visit(:);
X  = S';                             % [nWin x nFreq]
mu = mean(X, 1);
Xc = X - mu;

% --- PCA denoising (rank chosen by cumulative variance) -------------------
[U, Sg, V] = svd(Xc, 'econ');
ev    = diag(Sg).^2;
nComp = find(cumsum(ev) / sum(ev) >= o.VarKeep - 1e-12, 1, 'first');
if isempty(nComp), nComp = numel(ev); end           % VarKeep >= 1: keep everything
Xden  = U(:, 1:nComp) * Sg(1:nComp, 1:nComp) * V(:, 1:nComp)' + mu;

% --- k-means ---------------------------------------------------------------
rng(0);
if isempty(o.Start)
    raw = kmeans(Xden, k, 'Replicates', o.Replicates, 'MaxIter', 1000, ...
                 'Distance', 'sqeuclidean');
else
    raw = kmeans(Xden, k, 'Start', o.Start, 'MaxIter', 1000, ...
                 'Distance', 'sqeuclidean');
end

% --- order labels along the baseline -> endpoint axis ---------------------
visits = unique(visit);
occRaw = zeros(numel(visits), k);
for v = 1:numel(visits)
    w = visit == visits(v);
    occRaw(v, :) = accumarray(raw(w), 1, [k 1])' / sum(w);
end
iB = find(visits == o.Baseline);
iE = find(visits == o.Endpoint);

order = zeros(1, k);
[~, order(1)] = max(occRaw(iB, :));                 % state 1 = baseline-dominant
rest = setdiff(1:k, order(1));
[~, ie] = max(occRaw(iE, rest));
order(k) = rest(ie);                                 % state k = endpoint-dominant
rest = setdiff(rest, order(k));
[~, ib] = sort(occRaw(iB, rest), 'descend');
order(2:k-1) = rest(ib);                             % middle by baseline occupancy

labels = zeros(size(raw));
for s = 1:k
    labels(raw == order(s)) = s;
end

% --- outputs ---------------------------------------------------------------
out.labels     = labels;
out.occupancy  = occRaw(:, order);
out.silhouette = mean(silhouette(Xden, labels));
out.nComp      = nComp;
out.centroids  = zeros(k, size(S, 1));
for s = 1:k
    out.centroids(s, :) = mean(X(labels == s, :), 1);
end

% --- example figure --------------------------------------------------------
if o.Plot
    figure('Color', 'w');
    area(visits, 100 * out.occupancy);
    xlabel('Visit'); ylabel('Occupancy (%)');
    legend(arrayfun(@(s) sprintf('State %d', s), 1:k, 'uni', 0), ...
           'Location', 'eastoutside');
    title(sprintf('Ordered state occupancy (silhouette = %.2f, %d PCs)', ...
                  out.silhouette, nComp));
end
end

