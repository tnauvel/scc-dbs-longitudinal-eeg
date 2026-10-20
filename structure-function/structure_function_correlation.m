function out = structure_function_correlation(fa, eeg, names, varargin)
% STRUCTURE_FUNCTION_CORRELATION  Rank correlation between a structural
% measure and a family of EEG measures, with FDR control within the family.
%
% METHOD
%   Spearman correlation between a per-subject structural value (e.g.
%   cingulum bundle fractional anisotropy) and each EEG measure (e.g.
%   baseline whole-head coherence per band). P-values are corrected by the
%   Benjamini-Hochberg procedure within the supplied family. Subjects with
%   a missing value in either variable are dropped pairwise. A
%   leave-one-out range of rho is returned so that a result driven by one
%   subject can be recognised.
%
% REFERENCES
%   Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300
%   Nauvel T et al., in preparation (cingulum FA and baseline coherence)
%
% INPUTS
%   fa     [nSubj x 1]      structural measure
%   eeg    [nSubj x nMeas]  EEG measures, one column each
%   names  cell {nMeas}     column labels for the report and figure
%   Name-value:
%     'Alpha'  FDR level (default 0.05)
%     'Plot'   true to draw the example band-profile figure (default false)
%
% OUTPUT struct (one row per measure)
%   .rho, .p, .q, .n, .rho_loo_min, .rho_loo_max, .names
%
% AUTHORSHIP
%   Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at
%   Mount Sinai. Method and original implementation: T. Nauvel, 2022-2026
%   (full analysis pipeline in a private repository, see the README).
%   This release version was refactored into a standalone function with the
%   assistance of Claude (Anthropic), October 2026, then reviewed and
%   validated by the author against the original pipeline:
%   Spearman, Benjamini-Hochberg and leave-one-out logic match the original screen; checked on a synthetic six-band family.
%   Licence: MIT.

p = inputParser;
p.addParameter('Alpha', 0.05);
p.addParameter('Plot', false);
p.parse(varargin{:});
o = p.Results;

fa = fa(:);
nM = size(eeg, 2);
rho = nan(nM, 1); pv = nan(nM, 1); n = zeros(nM, 1);
lo = nan(nM, 1); hi = nan(nM, 1);
for m = 1:nM
    ok = ~isnan(fa) & ~isnan(eeg(:, m));
    n(m) = sum(ok);
    if n(m) < 4, continue; end
    x = fa(ok); y = eeg(ok, m);
    [rho(m), pv(m)] = corr(x, y, 'type', 'Spearman');
    r = zeros(n(m), 1);
    for i = 1:n(m)
        keep = true(n(m), 1); keep(i) = false;
        r(i) = corr(x(keep), y(keep), 'type', 'Spearman');
    end
    lo(m) = min(r); hi(m) = max(r);
end

% Benjamini-Hochberg within the family
q = nan(nM, 1);
v = find(~isnan(pv));
[ps, idx] = sort(pv(v));
k = numel(ps);
qs = min(1, ps .* k ./ (1:k)');
for i = k-1:-1:1, qs(i) = min(qs(i), qs(i+1)); end
q(v(idx)) = qs;

out.names = names(:);
out.rho = rho; out.p = pv; out.q = q; out.n = n;
out.rho_loo_min = lo; out.rho_loo_max = hi;
out.significant = q <= o.Alpha;

if o.Plot
    figure('Color', 'w'); hold on;
    bar(rho, 'FaceColor', [0.27 0.67 0.60]);
    errorbar(1:nM, rho, rho - lo, hi - rho, 'k.', 'LineWidth', 1);
    sig = find(out.significant);
    plot(sig, rho(sig) + sign(rho(sig)) * 0.06, 'k*');
    set(gca, 'XTick', 1:nM, 'XTickLabel', names, 'XTickLabelRotation', 30);
    ylabel('Spearman \rho'); ylim([-1 1]); yline(0, 'k:');
    title('Structural measure vs EEG family (bars: \rho, whiskers: leave-one-out range, *: FDR q \leq \alpha)');
    hold off;
end
end

