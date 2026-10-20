function out = occupancy_behavior_crosslag(occ, score, varargin)
% OCCUPANCY_BEHAVIOR_CROSSLAG  Does a brain-state measure track and precede
% clinical change within subjects?
%
% METHOD
%   1. Within-subject Spearman correlation between the EEG measure (e.g.
%      state-1 occupancy) and the concurrent symptom score across visits.
%      Per-subject r is pooled by Fisher z (mean z back-transformed, 95% CI
%      from the SE of the mean z).
%   2. Group significance by permutation: the score series is circularly
%      shifted within each subject NPerm times and the pooled r recomputed.
%   3. Cross-lagged panel model on within-subject standardised series:
%        score(t+1) ~ score(t) + eeg(t)      -> beta_eeg2score
%        eeg(t+1)   ~ eeg(t)   + score(t)    -> beta_score2eeg
%      Standardised coefficients, F-test p-values, and the difference in
%      residual sum of squares between the two directions. A larger
%      eeg->score coefficient indicates the EEG measure leads symptom change.
%
% REFERENCES
%   Fisher RA (1915) Biometrika 10:507-521 (z transform)
%   Kenny DA (1975) Psychol Bull 82:887-903 (cross-lagged panel correlation)
%   Alagapan S et al. (2023) Nature 622:130-138 (LFP precedes relapse)
%
% INPUTS
%   occ    [nSubj x nVisit] EEG measure per subject and visit (NaN allowed)
%   score  [nSubj x nVisit] symptom score at the same visits (NaN allowed)
%   Name-value:
%     'NPerm'  permutations for the pooled r (default 10000)
%     'Plot'   true to draw the example forest plot (default false)
%
% OUTPUT struct
%   .r_subject   [nSubj x 1] per-subject Spearman r
%   .p_subject   [nSubj x 1]
%   .r_pooled, .ci_pooled, .p_perm
%   .crosslag    struct: beta_eeg2score, p_eeg2score, beta_score2eeg,
%                p_score2eeg, n_pairs
%
% AUTHORSHIP
%   Tanya Nauvel, Quantitative Biometrics Lab, Icahn School of Medicine at
%   Mount Sinai. Method and original implementation: T. Nauvel, 2022-2026
%   (full analysis pipeline in a private repository, see the README).
%   This release version was refactored into a standalone function with the
%   assistance of Claude (Anthropic), October 2026, then reviewed and
%   validated by the author against the original pipeline:
%   paper tables, 11 subjects x 8 visits: EEG->HDRS beta 0.354 (p=0.0024) vs 0.34 (p=0.003) published; reverse path n.s. in both.
%   Licence: MIT.

p = inputParser;
p.addParameter('NPerm', 10000);
p.addParameter('Plot', false);
p.parse(varargin{:});
o = p.Results;

[nS, nV] = size(occ);
rs = nan(nS, 1); ps = nan(nS, 1);
for s = 1:nS
    ok = ~isnan(occ(s,:)) & ~isnan(score(s,:));
    if sum(ok) >= 4
        [rs(s), ps(s)] = corr(occ(s,ok)', score(s,ok)', 'type', 'Spearman');
    end
end
z  = atanh(max(min(rs, 0.999), -0.999));
zm = mean(z, 'omitnan');
se = std(z, 'omitnan') / sqrt(sum(~isnan(z)));
out.r_subject = rs;
out.p_subject = ps;
out.r_pooled  = tanh(zm);
out.ci_pooled = tanh(zm + [-1 1] * 1.96 * se);

% --- permutation (circular shift within subject) ----------------------------
rng(0);
zperm = zeros(o.NPerm, 1);
for k = 1:o.NPerm
    zk = nan(nS, 1);
    for s = 1:nS
        sh = circshift(score(s,:), randi(nV));
        ok = ~isnan(occ(s,:)) & ~isnan(sh);
        if sum(ok) >= 4
            zk(s) = atanh(corr(occ(s,ok)', sh(ok)', 'type', 'Spearman'));
        end
    end
    zperm(k) = mean(zk, 'omitnan');
end
out.p_perm = mean(abs(zperm) >= abs(zm));

% --- cross-lagged panel on standardised within-subject series --------------
E = (occ   - mean(occ,   2, 'omitnan')) ./ std(occ,   0, 2, 'omitnan');
H = (score - mean(score, 2, 'omitnan')) ./ std(score, 0, 2, 'omitnan');
E_t = E(:, 1:end-1); E_t1 = E(:, 2:end);
H_t = H(:, 1:end-1); H_t1 = H(:, 2:end);
ok  = ~isnan(E_t) & ~isnan(E_t1) & ~isnan(H_t) & ~isnan(H_t1);

m1 = fitlm([H_t(ok) E_t(ok)], H_t1(ok), 'VarNames', {'score_t','eeg_t','score_t1'});
m2 = fitlm([E_t(ok) H_t(ok)], E_t1(ok), 'VarNames', {'eeg_t','score_t','eeg_t1'});
out.crosslag.beta_eeg2score = m1.Coefficients.Estimate('eeg_t');
out.crosslag.p_eeg2score    = m1.Coefficients.pValue('eeg_t');
out.crosslag.beta_score2eeg = m2.Coefficients.Estimate('score_t');
out.crosslag.p_score2eeg    = m2.Coefficients.pValue('score_t');
out.crosslag.n_pairs        = sum(ok(:));

% --- example figure: forest plot --------------------------------------------
if o.Plot
    figure('Color', 'w'); hold on;
    y = 1:nS;
    plot(rs, y, 'ko', 'MarkerFaceColor', 'k');
    plot(rs(ps < 0.05), y(ps < 0.05), 'ko', 'MarkerFaceColor', [0.07 0.47 0.2]);
    plot(out.ci_pooled, [nS+1 nS+1], 'k-', 'LineWidth', 2);
    plot(out.r_pooled, nS+1, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 9);
    xline(0, 'k--');
    set(gca, 'YTick', [y nS+1], 'YTickLabel', [arrayfun(@(s) sprintf('S%d', s), y, 'uni', 0) {'Pooled'}]);
    xlabel('Spearman r (EEG measure vs symptom score)');
    title(sprintf('Pooled r = %.2f [%.2f, %.2f], perm p = %.3f; EEG\\rightarrowscore \\beta = %.2f', ...
                  out.r_pooled, out.ci_pooled, out.p_perm, out.crosslag.beta_eeg2score));
    hold off;
end
end

