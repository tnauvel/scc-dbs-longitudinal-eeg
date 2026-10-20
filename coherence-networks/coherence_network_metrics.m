function out = coherence_network_metrics(W, labels, varargin)
% COHERENCE_NETWORK_METRICS  Graph-theoretic summary of an EEG coherence matrix.
%
% METHOD
%   The symmetric coherence matrix is thresholded either proportionally
%   (keep the strongest fraction of edges) or reduced to its maximum spanning
%   tree (threshold-free). Global efficiency, characteristic path length,
%   clustering coefficient and Louvain modularity are computed on the
%   binarised graph; mean strength is the mean coherence of retained edges.
%   Regional coherence means (left, right, cross-hemispheric, midline) are
%   computed from the unthresholded matrix using the 10-20 label convention
%   (odd number = left, even = right, 'z' = midline).
%
% REFERENCES
%   Rubinov M, Sporns O (2010) NeuroImage 52:1059-1069 (Brain Connectivity Toolbox)
%   Tewarie P et al. (2015) NeuroImage 104:177-188 (minimum spanning tree)
%   van den Heuvel MP et al. (2017) NeuroImage 152:437-449 (proportional thresholds)
%   Humphries MD, Gurney K (2008) PLoS ONE 3:e0002051 (small-world sigma)
%
% INPUTS
%   W       [n x n] symmetric coherence matrix, zeros on the diagonal.
%           For a measure that can be negative (e.g. debiased wPLI), pass
%           'Negative','zero' to zero negatives (literature default).
%   labels  cell array of n channel labels, e.g. {'FP1','FZ','C4',...}
%   Name-value:
%     'Threshold'  fraction of edges kept (default 0.25) or 'mst'
%     'Negative'   'abs' | 'zero' | 'none' (default 'none')
%     'SmallWorld' true to compute sigma against degree-preserving nulls
%     'NRand'      number of null networks (default 20)
%     'Plot'       true to draw the example thresholded-adjacency figure
%
% OUTPUT struct
%   .global_efficiency, .path_length, .clustering, .modularity,
%   .mean_strength, .sigma (NaN unless SmallWorld)
%   .region  struct with .left .right .cross .midline mean coherence
%   .A       binarised adjacency used for the metrics
%
% Requires Brain Connectivity Toolbox on the path (threshold_proportional,
% efficiency_bin, distance_bin, charpath, clustering_coef_bu,
% community_louvain, randmio_und).
%
% AUTHORSHIP
%   Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at
%   Mount Sinai. Method and original implementation: T. Nauvel, 2022-2026
%   (full analysis pipeline in a private repository, see the README).
%   This release version was refactored into a standalone function with the
%   assistance of Claude (Anthropic), October 2026, then reviewed and
%   validated by the author against the original pipeline:
%   metric definitions follow the Brain Connectivity Toolbox calls used in the original pipeline; thresholding and regional means checked on a synthetic 19-channel network.
%   Licence: MIT.

p = inputParser;
p.addParameter('Threshold', 0.25);
p.addParameter('Negative', 'none');
p.addParameter('SmallWorld', false);
p.addParameter('NRand', 20);
p.addParameter('Plot', false);
p.parse(varargin{:});
o = p.Results;

n = size(W, 1);
W = (W + W') / 2;
W(1:n+1:end) = 0;
switch lower(o.Negative)
    case 'abs',  W = abs(W);
    case 'zero', W(W < 0) = 0;
end

% --- thresholding ----------------------------------------------------------
if ischar(o.Threshold) && strcmpi(o.Threshold, 'mst')
    G   = graph(1 - W, 'omitselfloops');            % distance = 1 - coherence
    T   = minspantree(G);
    Wt  = zeros(n);
    e   = T.Edges.EndNodes;
    for i = 1:size(e, 1)
        Wt(e(i,1), e(i,2)) = W(e(i,1), e(i,2));
        Wt(e(i,2), e(i,1)) = W(e(i,1), e(i,2));
    end
else
    Wt = threshold_proportional(W, o.Threshold);
end
A = double(Wt > 0);

% --- global metrics --------------------------------------------------------
out.global_efficiency = efficiency_bin(A);
D = distance_bin(A);
out.path_length       = charpath(D, 0, 0);
out.clustering        = mean(clustering_coef_bu(A));
[~, out.modularity]   = community_louvain(Wt);
out.mean_strength     = mean(Wt(A > 0));
out.sigma             = NaN;

if o.SmallWorld
    Cr = zeros(o.NRand, 1); Lr = zeros(o.NRand, 1);
    for r = 1:o.NRand
        Ar = randmio_und(A, 10);
        Cr(r) = mean(clustering_coef_bu(Ar));
        Lr(r) = charpath(distance_bin(Ar), 0, 0);
    end
    out.sigma = (out.clustering / mean(Cr)) / (out.path_length / mean(Lr));
end

% --- regional means from the unthresholded matrix --------------------------
side = zeros(n, 1);                                    % -1 left, 0 midline, +1 right
for i = 1:n
    lab = upper(char(labels{i}));
    d   = regexp(lab, '\d+$', 'match', 'once');
    if isempty(d) || lab(end) == 'Z'
        side(i) = 0;
    elseif mod(str2double(d), 2) == 1
        side(i) = -1;
    else
        side(i) = 1;
    end
end
[I, J] = find(triu(true(n), 1));
pairVal = W(sub2ind([n n], I, J));
out.region.left    = mean(pairVal(side(I) == -1 & side(J) == -1));
out.region.right   = mean(pairVal(side(I) ==  1 & side(J) ==  1));
out.region.cross   = mean(pairVal(side(I) .* side(J) == -1));
out.region.midline = mean(pairVal(side(I) ==  0 | side(J) ==  0));
out.A = A;

% --- example figure --------------------------------------------------------
if o.Plot
    figure('Color', 'w');
    imagesc(Wt); axis square; colorbar;
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'YTick', 1:n, ...
             'YTickLabel', labels, 'XTickLabelRotation', 90, 'FontSize', 6);
    title(sprintf('Thresholded coherence  (E_{glob} = %.2f, Q = %.2f)', ...
                  out.global_efficiency, out.modularity));
end
end

