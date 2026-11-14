function runModel(cellName)
% RUNMODEL  Simulate the baseline model for one cell, measure it with the
% same extractFeatures.m used on the recording, and compare the two.
%
%   runModel('pyr')     pyramidal cell (Allen 329549957)
%   runModel('fs')      fast spiking cell (Allen 487537399)

switch lower(cellName)
    case 'pyr'
        p = paramsPyr();  file = '329549955_ephys.nwb';  label = 'Pyramidal';
        showAmps = [-110 130 250];
    case 'fs'
        p = paramsFS();   file = '487537397_ephys.nwb';  label = 'Fast spiking';
        showAmps = [-90 270 330];
    otherwise
        error('Use runModel(''pyr'') or runModel(''fs'')');
end

data  = readFiles(file);
Fdat  = extractFeatures(data);

amps  = unique([data.amp]);
tic;  model = hhModel(p, amps);  fprintf('Simulated %d steps in %.1f s\n', numel(amps), toc);
Fmod  = extractFeatures(model);

names = {'rheobase','fi_slope','adaptation','rin','tau','vrest', ...
         'threshold','amplitude','width','ahp'};
units = {'pA','Hz/pA','','MOhm','ms','mV','mV','mV','ms','mV'};
fprintf('\n===== %s: data vs model =====\n', label);
fprintf('%-12s %10s %10s %8s\n', 'Feature', 'Data', 'Model', 'Units');
for k = 1:numel(names)
    fprintf('%-12s %10.3g %10.3g %8s\n', names{k}, Fdat.(names{k}), Fmod.(names{k}), units{k});
end

figure('Name', [label ' fit'], 'Color', 'w', 'Position', [100 100 1100 650]);
subplot(2, 3, [1 4]); hold on
plot(Fdat.fi_amps, Fdat.fi_rate, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Data');
plot(Fmod.fi_amps, Fmod.fi_rate, 'rs--', 'LineWidth', 1.5, 'DisplayName', 'Model');
xlabel('Current (pA)'); ylabel('Firing rate (Hz)'); title('f-I curve'); legend('Location', 'northwest'); box on

for k = 1:numel(showAmps)
    subplot(2, 3, k + 1 + (k > 2));
    jd = find([data.amp]  == showAmps(k), 1);
    jm = find([model.amp] == showAmps(k), 1);
    plot(data(jd).t - data(jd).tOn,  data(jd).v,  'k', 'DisplayName', 'Data'); hold on
    plot(model(jm).t - model(jm).tOn, model(jm).v, 'r', 'DisplayName', 'Model');
    xlim([-0.1 1.1]); ylim([-110 50]);
    xlabel('Time from step onset (s)'); ylabel('V_m (mV)');
    title(sprintf('%d pA', showAmps(k))); box on
    if k == 1, legend('Location', 'southeast'); end
end
end
