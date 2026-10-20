%% example_usage.m
% Synthetic demonstration of coherence_network_metrics on a 19-channel 10-20
% montage. Requires the Brain Connectivity Toolbox on the MATLAB path.

rng(3);
labels = {'FP1','FP2','F7','F3','FZ','F4','F8','T7','C3','CZ','C4','T8', ...
          'P7','P3','PZ','P4','P8','O1','O2'};
n = numel(labels);
W = 0.15 + 0.10 * rand(n);  W = (W + W') / 2;  W(1:n+1:end) = 0;
% add a stronger within-hemisphere module on each side
for i = 1:n
    for j = i+1:n
        li = labels{i}(end); lj = labels{j}(end);
        if ~isletter(li) && ~isletter(lj) && mod(str2double(li), 2) == mod(str2double(lj), 2)
            W(i, j) = W(i, j) + 0.25;  W(j, i) = W(i, j);
        end
    end
end

out = coherence_network_metrics(W, labels, 'Threshold', 0.25, 'SmallWorld', true, 'Plot', true);
fprintf('global efficiency %.3f, path length %.2f, clustering %.3f, modularity %.3f, sigma %.2f\n', ...
    out.global_efficiency, out.path_length, out.clustering, out.modularity, out.sigma);
fprintf('regional coherence: left %.3f, right %.3f, cross %.3f, midline %.3f\n', ...
    out.region.left, out.region.right, out.region.cross, out.region.midline);

% Threshold-free alternative
out_mst = coherence_network_metrics(W, labels, 'Threshold', 'mst');
fprintf('MST: global efficiency %.3f, modularity %.3f\n', out_mst.global_efficiency, out_mst.modularity);

