function pop = popFit(cellName, nRounds, nPerRound)
% POPFIT  Find a population of parameter sets that all fit one Allen cell.
%
%   pop = popFit('pyr');            % pyramidal cell
%   pop = popFit('fs', 8, 400);     % fast spiking, 8 rounds of 400 models
%
% Method (random walk through the space of acceptable models):
%   1. Start the pool with the baseline model (paramsPyr / paramsFS).
%   2. Each round, pick models from the pool at random and nudge every
%      varied parameter a little (conductances by ~15%, voltages by ~1 mV).
%   3. Simulate the new models and keep every one whose features all fall
%      within tolerance of the recorded cell.
%   4. Add the kept models to the pool and repeat.
% Because new models only come from models that already fit, the pool
% spreads out through the region of parameter space that fits the cell.
%
% Results are saved to pop_<cell>.mat.

if nargin < 2, nRounds   = 8;   end
if nargin < 3, nPerRound = 400; end
batch = 150;
rng(1);

switch lower(cellName)
    case 'pyr'
        base = paramsPyr();  file = '329549955_ephys.nwb';
        logVars = {'gNa','gK','gCa','gAHP','gh','gL','tauCa','area'};
        linVars = {'EL','VT'};
        amps    = [-110 -50 90 110 130 170 210 250];
        heroAmp = 170;
        tol = struct('rheo',{'abs',20}, 'r170',{'abs',2}, 'r250',{'abs',3}, ...
                     'adapt',{'abs',0.05}, 'rin',{'rel',0.15}, 'tau',{'rel',0.20}, ...
                     'vrest',{'abs',3}, 'thr',{'abs',3}, 'amp',{'rel',0.10}, ...
                     'width',{'rel',0.15}, 'ahp',{'abs',4});
    case 'fs'
        base = paramsFS();   file = '487537397_ephys.nwb';
        base.gCa = 0.1;
        logVars = {'gNa','gK','gKv3','gL','area','q','tw1','gCa'};
        linVars = {'EL','VT','wV'};
        amps    = [-90 -50 190 230 270 310 330];
        heroAmp = 270;
        tol = struct('rheo',{'abs',40}, 'r310',{'rel',0.20}, 'r330',{'rel',0.20}, ...
                     'adapt',{'abs',0.05}, 'rin',{'rel',0.15}, ...
                     'vrest',{'abs',3}, 'thr',{'abs',3}, 'amp',{'rel',0.15}, ...
                     'width',{'rel',0.25}, 'ahp',{'abs',4});
    otherwise
        error('Use popFit(''pyr'') or popFit(''fs'')');
end
tolNames = fieldnames(tol);

data = readFiles(file);
D = popFeatures(pickSweeps(data, amps), amps, heroAmp);
fprintf('\nTargets for %s (from the recording):\n', cellName);
for k = 1:numel(tolNames)
    fprintf('  %-6s %8.3g  (+/- %g%s)\n', tolNames{k}, D.(tolNames{k}), ...
        tol(2).(tolNames{k}), ifelse(strcmp(tol(1).(tolNames{k}),'rel'),' x value',''));
end

names = fieldnames(base);
pool  = struct2table(base);
if strcmp(cellName, 'fs'), pool.gAHP = 0; end
poolF = {};

for r = 1:nRounds
    tic;
    idx = randi(height(pool), nPerRound, 1);
    P = pool(idx, :);
    for k = 1:numel(logVars)
        P.(logVars{k}) = P.(logVars{k}) .* exp(0.15 * randn(nPerRound, 1));
    end
    for k = 1:numel(linVars)
        P.(linVars{k}) = P.(linVars{k}) + 1.0 * randn(nPerRound, 1);
    end
    if strcmp(cellName, 'fs')
        g = P.gAHP;  z = (g == 0);
        on = z & (rand(nPerRound, 1) < 0.5);
        g(on) = exp(log(0.002) + rand(sum(on), 1) * (log(0.05) - log(0.002)));
        nz = ~z;
        g(nz) = g(nz) .* exp(0.3 * randn(sum(nz), 1));
        g(nz & rand(nPerRound, 1) < 0.1) = 0;
        P.gAHP = g;
    end

    keep = false(nPerRound, 1);  Fr = cell(nPerRound, 1);
    for b0 = 1:batch:nPerRound
        b1 = min(b0 + batch - 1, nPerRound);  nb = b1 - b0 + 1;
        pc = struct();
        for k = 1:numel(names)
            pc.(names{k}) = repelem(P.(names{k})(b0:b1)', numel(amps));
        end
        sw = hhModel(pc, repmat(amps, 1, nb));
        for i = 1:nb
            Fi = popFeatures(sw((i-1)*numel(amps) + (1:numel(amps))), amps, heroAmp);
            Fr{b0+i-1} = Fi;
            keep(b0+i-1) = passes(Fi, D, tol, tolNames);
        end
    end
    pool  = [pool; P(keep, :)];
    poolF = [poolF; Fr(keep)];
    fprintf('Round %d: kept %d of %d, pool = %d models (%.0f s)\n', ...
        r, sum(keep), nPerRound, height(pool), toc);
end

fprintf('\nRange of each varied parameter across the %d accepted models:\n', height(pool));
for k = 1:numel(logVars)
    v = pool.(logVars{k});  v = v(v > 0);
    fprintf('  %-6s %10.4g to %-10.4g (%.1f-fold)\n', logVars{k}, min(v), max(v), max(v)/min(v));
end
for k = 1:numel(linVars)
    v = pool.(linVars{k});
    fprintf('  %-6s %8.1f to %-8.1f mV\n', linVars{k}, min(v), max(v));
end
if strcmp(cellName, 'fs')
    fprintf('  gAHP   zero in %.0f%% of models, max %.4f\n', 100*mean(pool.gAHP == 0), max(pool.gAHP));
end

pop.params  = pool;
pop.feats   = poolF;
pop.targets = D;
pop.tol     = tol;
pop.varied  = [logVars, linVars];
pop.amps    = amps;
save(sprintf('pop_%s.mat', lower(cellName)), 'pop');
end

function S = pickSweeps(data, amps)
for k = 1:numel(amps)
    S(k) = data(find([data.amp] == amps(k), 1));
end
end

function ok = passes(F, D, tol, names)
ok = true;
for k = 1:numel(names)
    x = F.(names{k});  d = D.(names{k});  lim = tol(2).(names{k});
    if strcmp(tol(1).(names{k}), 'rel'), lim = lim * abs(d); end
    if isnan(x) || abs(x - d) > lim, ok = false; return, end
end
end

function F = popFeatures(sw, amps, heroAmp)
F = struct('rheo',NaN,'adapt',NaN,'rin',NaN,'tau',NaN,'vrest',NaN, ...
           'thr',NaN,'amp',NaN,'width',NaN,'ahp',NaN);
n = numel(sw);  rate = zeros(1, n);  tsp = cell(1, n);
for j = 1:n
    s  = sw(j);
    in = s.t >= s.tOn & s.t <= s.tOff;
    up = find(diff([false; (s.v > -20) & in]) == 1);
    tsp{j}  = s.t(up);
    rate(j) = numel(up) / (s.tOff - s.tOn);
    if amps(j) > 0, F.(sprintf('r%d', amps(j))) = rate(j); end
end

base = zeros(1, n);
for j = 1:n
    w = sw(j).t > sw(j).tOn - 0.11 & sw(j).t < sw(j).tOn - 0.01;
    base(j) = mean(sw(j).v(w));
end
F.vrest = median(base);
neg = find(amps < 0);
dV = zeros(size(neg));
for m = 1:numel(neg)
    s = sw(neg(m));
    dV(m) = mean(s.v(s.t > s.tOff - 0.1 & s.t <= s.tOff)) - base(neg(m));
end
p = polyfit(amps(neg), dV, 1);  F.rin = p(1) * 1e3;
s = sw(amps == -50);
k0 = find(s.t >= s.tOn, 1);
seg = movmean(s.v(k0:k0 + round(0.3*s.fs)), max(1, round(s.fs / 2000))) - base(amps == -50);
F.tau = (find(seg <= 0.632 * min(seg), 1) - 1) / s.fs * 1e3;

pos = find(amps > 0 & rate > 0);
if isempty(pos), return, end
F.rheo = amps(pos(1));
isi = diff(tsp{amps == heroAmp});
if numel(isi) > 1
    F.adapt = mean(diff(isi) ./ (isi(2:end) + isi(1:end-1)));
end

s = sw(pos(1));  v = s.v;  fs = s.fs;
dvdt = gradient(movmean(v, max(1, round(fs / 22000)))) * fs / 1000;
in = s.t >= s.tOn & s.t <= s.tOff;
i0 = find(v > -20 & in, 1);
[vp, ipk] = max(v(i0:min(i0 + round(0.003*fs), numel(v))));  ipk = ipk + i0 - 1;
lo = ipk - round(0.002*fs);
[~, imx] = max(dvdt(lo:ipk));  it = imx + lo - 1;
while it > lo && dvdt(it) >= 20, it = it - 1; end
vt = v(it);  half = (vt + vp) / 2;
ir = find(v(it:ipk) >= half, 1) + it - 1;
iff = find(v(ipk:end) < half, 1) + ipk - 1;
F.thr = vt;  F.amp = vp - vt;
F.width = (iff - ir) / fs * 1e3;
F.ahp = min(v(ipk:min(ipk + round(0.05*fs), numel(v)))) - vt;
end

function out = ifelse(c, a, b)
if c, out = a; else, out = b; end
end
