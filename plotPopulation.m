function plotPopulation()
% PLOTPOPULATION  Figures for Intermediate 1 from pop_pyr.mat and pop_fs.mat.
%   Figure4_parameter_spread.png: how much each parameter varies (report Figure 2)
%   Figure6_parameter_correlations.png: parameter correlations (report Figure 3)
%   Figure5_population_fits.png: accepted models vs the recording (report Figure 4)
% Prints a summary to popSummary_output.txt and saves PNGs.

cells = {'pyr', 'Pyramidal',    '329549955_ephys.nwb', 170;
         'fs',  'Fast spiking', '487537397_ephys.nwb', 310};
diary('popSummary_output.txt'); diary on
rng(2);

f4 = figure('Color', 'w', 'Position', [60 60 1150 430]);
for c = 1:2
    S = load(sprintf('pop_%s.mat', cells{c,1}));  pop = S.pop;
    P = pop.params;  vars = pop.varied;
    base = P(1, :);
    fprintf('\n===== %s: %d accepted models =====\n', cells{c,2}, height(P));
    subplot(1, 2, c); hold on
    for k = 1:numel(vars)
        v = P.(vars{k});
        if any(strcmp(vars{k}, {'EL','VT','wV'}))
            y = 2 .^ ((v - base.(vars{k})) / 3);
            fprintf('  %-6s %7.1f to %7.1f mV (range %.1f mV)\n', vars{k}, min(v), max(v), max(v)-min(v));
        else
            y = v / base.(vars{k});
            fprintf('  %-6s %.1f-fold range\n', vars{k}, max(v)/min(v));
        end
        x = k + 0.25 * (rand(size(y)) - 0.5);
        scatter(x, y, 10, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.35);
        plot(k + [-0.3 0.3], [1 1] * median(y), 'k-', 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log', 'XTick', 1:numel(vars), 'XTickLabel', vars, 'XTickLabelRotation', 45);
    yline(1, ':');  ylim([0.3 3.5]);  xlim([0.4 numel(vars) + 0.6]);
    set(gca, 'YTick', [0.25 0.5 1 2 4], 'YTickLabel', {'0.25x','0.5x','1x','2x','4x'});
    ylabel('Value relative to baseline (log scale)');
    title(sprintf('%s (n = %d)', cells{c,2}, height(P))); box on
    if strcmp(cells{c,1}, 'fs')
        fprintf('  gAHP   zero in %.0f%% of models; nonzero range %.4f to %.4f\n', ...
            100*mean(P.gAHP == 0), min(P.gAHP(P.gAHP > 0)), max(P.gAHP));
    end
end
exportgraphics(f4, 'Figure4_parameter_spread.png', 'Resolution', 300);

f5 = figure('Color', 'w', 'Position', [60 60 1150 650]);
for c = 1:2
    S = load(sprintf('pop_%s.mat', cells{c,1}));  pop = S.pop;  P = pop.params;
    data = readFiles(cells{c,3});
    amps = unique([data.amp]);  amps = amps(amps > 0);
    pick = randperm(height(P), min(20, height(P)));
    names = P.Properties.VariableNames;  pc = struct();
    for k = 1:numel(names)
        pc.(names{k}) = repelem(P.(names{k})(pick)', numel(amps));
    end
    sw = hhModel(pc, repmat(amps, 1, numel(pick)));
    rates = reshape(cellfun(@(v, t0, t1, t) sum(diff([false; v > -20 & t >= t0 & t <= t1]) == 1), ...
        {sw.v}, {sw.tOn}, {sw.tOff}, {sw.t}) / 1, numel(amps), []);
    dr = arrayfun(@(a) mean(arrayfun(@(s) sum(diff([false; s.v > -20 & s.t >= s.tOn & s.t <= s.tOff]) == 1), ...
        data([data.amp] == a))), amps);

    subplot(2, 2, c); hold on
    hm = plot(amps, rates, '-', 'Color', [0.85 0.3 0.3 0.35]);
    hd = plot(amps, dr, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    xlabel('Current (pA)'); ylabel('Firing rate (Hz)');
    title(sprintf('%s: f-I curves', cells{c,2}));
    legend([hd hm(1)], {'Data', sprintf('%d accepted models', numel(pick))}, 'Location', 'northwest');
    box on

    subplot(2, 2, c + 2); hold on
    a = cells{c,4};
    jd = find([data.amp] == a, 1);
    for i = 1:min(5, numel(pick))
        s = sw((i-1)*numel(amps) + find(amps == a));
        plot(s.t - s.tOn, s.v - 130*(i), 'Color', [0.85 0.3 0.3]);
    end
    plot(data(jd).t - data(jd).tOn, data(jd).v, 'k');
    text(1.0, max(data(jd).v) + 15, 'Data', 'HorizontalAlignment', 'right');
    text(1.0, -130 + 45, 'Models', 'HorizontalAlignment', 'right', 'Color', [0.85 0.3 0.3]);
    ylim([-130*min(5, numel(pick)) - 100, max(data(jd).v) + 30]);
    xlim([-0.05 1.05]); set(gca, 'YTick', []);
    xlabel('Time from step onset (s)');
    title(sprintf('%s, %d pA step', cells{c,2}, a)); box on
end
exportgraphics(f5, 'Figure5_population_fits.png', 'Resolution', 300);

f6 = figure('Color', 'w', 'Position', [60 60 1150 480]);
for c = 1:2
    S = load(sprintf('pop_%s.mat', cells{c,1}));  pop = S.pop;  P = pop.params;
    vars = pop.varied;  X = zeros(height(P), numel(vars));
    for k = 1:numel(vars)
        v = P.(vars{k});
        if any(strcmp(vars{k}, {'EL','VT','wV'})), X(:,k) = v; else, X(:,k) = log(v); end
    end
    R = spearman(X);
    subplot(1, 2, c);
    imagesc(R, [-1 1]); axis square
    colormap(gca, redblue());  colorbar;
    set(gca, 'XTick', 1:numel(vars), 'XTickLabel', vars, 'XTickLabelRotation', 45, ...
             'YTick', 1:numel(vars), 'YTickLabel', vars);
    title(cells{c,2});
    U = triu(abs(R), 1);  [~, order] = sort(U(:), 'descend');
    fprintf('\n%s strongest parameter correlations:\n', cells{c,2});
    for q = 1:4
        [i, j] = ind2sub(size(R), order(q));
        fprintf('  %-6s vs %-6s  rho = %+.2f\n', vars{i}, vars{j}, R(i, j));
    end
end
exportgraphics(f6, 'Figure6_parameter_correlations.png', 'Resolution', 300);
diary off
end

function R = spearman(X)
Rk = zeros(size(X));
for k = 1:size(X, 2)
    [~, ~, Rk(:, k)] = unique(X(:, k));
end
R = corrcoef(Rk);
end

function m = redblue()
x = linspace(0, 1, 128)';
m = [[x; ones(128,1)], [x; flipud(x)], [ones(128,1); flipud(x)]];
end
