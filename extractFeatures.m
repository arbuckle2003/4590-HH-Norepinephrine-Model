function F = extractFeatures(sw)
% EXTRACTFEATURES  Electrophysiological features from long square sweeps.
%   F = extractFeatures(sw)   where sw comes from readFiles.m
% The same function is used on recorded data and on model output, so
% both are measured exactly the same way.
%
% Features (definitions as in the Methods):
%   rheobase    smallest positive step (pA) that evokes >= 1 spike
%   fi_slope    slope of firing rate vs current at and above rheobase (Hz/pA)
%   adaptation  mean of (ISI(k+1)-ISI(k))/(ISI(k+1)+ISI(k)) on the "hero"
%               sweep: lowest step >= rheobase + 40 pA with >= 4 spikes
%   threshold   V where dV/dt first exceeds 20 mV/ms (first spike at rheobase)
%   amplitude   peak minus threshold (mV)
%   width       full width at half amplitude (ms)
%   ahp         trough after the spike minus threshold (mV)
%   vrest       median V in the 100 ms before the step (mV)
%   rin         input resistance from small hyperpolarizing steps (MOhm)
%   tau         membrane time constant from the same steps (ms)

amps = [sw.amp];
n    = numel(sw);

sp  = cell(1, n);
for k = 1:n, sp{k} = detectSpikes(sw(k)); end
nsp  = cellfun(@(s) numel(s.tPk), sp);
dur  = [sw.tOff] - [sw.tOn];
rate = nsp ./ dur;  rate(isnan(rate)) = 0;
F.nspikes = nsp;  F.amps = amps;  F.rate = rate;

sup = find(nsp > 0 & amps > 0);
F.rheobase = min(amps(sup));

ua = unique(amps(amps >= F.rheobase));
mr = arrayfun(@(a) mean(rate(amps == a)), ua);
F.fi_amps = ua;  F.fi_rate = mr;
p = polyfit(ua, mr, 1);
F.fi_slope = p(1);

cand = find(amps >= F.rheobase + 40 & nsp >= 4);
if isempty(cand), cand = find(nsp >= 4); end
[~, j] = min(amps(cand));  h = cand(j);
isi = diff(sp{h}.tPk);
F.adaptation = mean(diff(isi) ./ (isi(2:end) + isi(1:end-1)));
F.hero_amp   = amps(h);
F.hero_rate  = rate(h);

r  = sup(find(amps(sup) == F.rheobase, 1));
s1 = sp{r};
F.threshold = s1.vThr(1);
F.peak      = s1.vPk(1);
F.amplitude = s1.vPk(1) - s1.vThr(1);
F.width     = s1.width(1);
F.ahp       = s1.vTr(1) - s1.vThr(1);

base = zeros(1, n);
for k = 1:n
    w = sw(k).t > sw(k).tOn - 0.11 & sw(k).t < sw(k).tOn - 0.01;
    base(k) = median(sw(k).v(w));
end
F.vrest = median(base);

neg = find(amps < 0 & amps >= -110);
if numel(neg) < 2, neg = find(amps < 0); end
dV  = zeros(size(neg));  taus = zeros(size(neg));
for m = 1:numel(neg)
    s  = sw(neg(m));
    ss = s.t > s.tOff - 0.1 & s.t <= s.tOff;
    dV(m)   = mean(s.v(ss)) - base(neg(m));
    taus(m) = fitTau(s);
end
p = polyfit(amps(neg), dV, 1);
F.rin = p(1) * 1e3;
F.tau = median(taus);
end

function S = detectSpikes(s)
S = struct('tPk', [], 'vPk', [], 'vThr', [], 'width', [], 'vTr', []);
if isnan(s.tOn), return, end
fs   = s.fs;  v = s.v(:);
dvdt = gradient(movmean(v, 9)) * fs / 1000;
in   = s.t(:) >= s.tOn & s.t(:) <= s.tOff + 0.002;
above = (v > -20) & in;
up = find(diff([false; above]) == 1);
dn = find(diff([above; false]) == -1);
prevPk = 1;
for m = 1:numel(up)
    [vpk, ipk] = max(v(up(m):dn(m)));  ipk = ipk + up(m) - 1;
    lo = max(ipk - round(0.002 * fs), prevPk + 1);
    [~, imx] = max(dvdt(lo:ipk));  imx = imx + lo - 1;
    ithr = imx;
    while ithr > lo && dvdt(ithr) >= 20, ithr = ithr - 1; end
    vthr = v(ithr);
    half  = (vpk + vthr) / 2;
    irise = find(v(ithr:ipk) >= half, 1) + ithr - 1;
    ifall = find(v(ipk:end) <  half, 1) + ipk - 1;
    if m < numel(up), nxt = up(m+1); else, nxt = min(ipk + round(0.05*fs), numel(v)); end
    S.tPk(end+1)   = s.t(ipk);
    S.vPk(end+1)   = vpk;
    S.vThr(end+1)  = vthr;
    S.width(end+1) = (ifall - irise) / fs * 1e3;
    S.vTr(end+1)   = min(v(ipk:nxt));
    prevPk = ipk;
end
end

function tau = fitTau(s)
k0  = find(s.t >= s.tOn, 1);
k1  = find(s.t >= s.tOn + 0.3, 1);
vs  = movmean(s.v, 21);
[~, im] = min(vs(k0:k1));  im = im + k0 - 1;
step = max(1, round(s.fs / 1e4));
idx = k0:step:im;
t   = s.t(idx) - s.t(k0);  y = vs(idx);
t = t(:); y = y(:);
if numel(t) < 5, tau = NaN; return, end
f   = @(q) sum((q(1) + q(2) * exp(-t / q(3)) - y).^2);
q0  = [y(end), y(1) - y(end), 0.01];
q   = fminsearch(f, q0, optimset('Display', 'off'));
tau = abs(q(3)) * 1e3;
end
