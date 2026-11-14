% runNE.m  Preliminary norepinephrine test on the baseline pyramidal model.
% NE is modeled as a fractional block r of the AHP conductance:
%   gAHP_NE = (1 - r) * gAHP,   r = 0, 0.25, 0.5, 0.75
% Adaptation is measured at a fixed step (170 pA, the data's hero sweep)
% so every condition is compared at the same input.

p0    = paramsPyr();
rList = [0 0.25 0.5 0.75];
amps  = [-110 -50 90 110 130 150 170 190 210 230 250];
aFix  = 170;
cols  = [0 0 0; 0.85 0.55 0.1; 0.85 0.3 0.1; 0.7 0 0];

diary('ne_output.txt'); diary on
fprintf('%-8s %10s %12s %12s %16s\n', 'Block', 'Rheobase', 'f-I slope', 'Adaptation', 'Spikes@170pA');
res = struct();
for k = 1:numel(rList)
    p = p0;  p.gAHP = (1 - rList(k)) * p0.gAHP;
    sw = hhModel(p, amps);
    F  = extractFeatures(sw);
    s  = sw([sw.amp] == aFix);
    up = find(diff([0; (s.v > -20) & s.t >= s.tOn & s.t <= s.tOff]) == 1);
    isi = diff(s.t(up));
    adapt = mean(diff(isi) ./ (isi(2:end) + isi(1:end-1)));
    res(k).r = rList(k);  res(k).F = F;  res(k).adapt = adapt;
    res(k).nsp = numel(up);  res(k).sw170 = s;
    fprintf('%-8s %8d pA %9.3f %12.3f %16d\n', sprintf('%d%%', 100*rList(k)), ...
            F.rheobase, F.fi_slope, adapt, numel(up));
end
diary off
save('ne_results.mat', 'res');

fig = figure('Name', 'NE preliminary', 'Color', 'w', 'Position', [100 100 1150 380]);

subplot(1, 3, 1); hold on
for k = 1:numel(rList)
    plot(res(k).F.amps(res(k).F.amps > 0), res(k).F.rate(res(k).F.amps > 0), 'o-', ...
        'Color', cols(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%d%% block', 100*rList(k)));
end
xlabel('Current (pA)'); ylabel('Firing rate (Hz)'); title('A  f-I curve'); legend('Location', 'northwest'); box on

subplot(1, 3, 2); hold on
plot(100*rList, [res.adapt], 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
xlabel('g_{AHP} block (%)'); ylabel('Adaptation index'); title(sprintf('B  Adaptation at %d pA', aFix));
xlim([-5 80]); box on

subplot(1, 3, 3); hold on
for k = [1 numel(rList)]
    s = res(k).sw170;
    plot(s.t - s.tOn, s.v, 'Color', cols(k,:), 'DisplayName', sprintf('%d%% block', 100*rList(k)));
end
xlim([-0.1 1.1]); ylim([-100 50]);
xlabel('Time from step onset (s)'); ylabel('V_m (mV)'); title(sprintf('C  Response to %d pA', aFix));
legend('Location', 'southeast'); box on

exportgraphics(fig, 'Figure3_NE_pyramidal.png', 'Resolution', 300);
