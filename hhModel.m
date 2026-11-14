function sweeps = hhModel(p, amps, dt)
% HHMODEL  Single-compartment Hodgkin-Huxley type neuron driven by 1 s
% current steps. All step amplitudes are simulated at once (one column
% of V per amplitude), and the output has the same format as readFiles.m
% so extractFeatures.m can be run on model and data in the same way.
%
%   sw = hhModel(paramsPyr(), [-110 -50 110 130 170 210 250]);
%
% p fields can be scalars (one model) or 1 x nS rows (one model per
% column), so many parameter sets can be simulated at once.
%
% Units: V in mV, t in ms, conductances in mS/cm^2, currents in uA/cm^2,
% C in uF/cm^2, [Ca] in uM, membrane area in cm^2, step amplitudes in pA.
%
% C dV/dt = I_inj - I_Na - I_Kdr - I_Kv3 - I_Ca - I_AHP - I_h - I_L

if nargin < 3, dt = 0.025; end
tOn  = 200;  tOff = 1200;  T = 1300;
rec  = 2;

amps = amps(:)';
nS   = numel(amps);
Iinj = amps * 1e-6 ./ p.area;

V  = p.EL .* ones(1, nS);
[am, bm, ah, bh, an, bn] = naK_rates(V, p);
m  = am ./ (am + bm);   h = ah ./ (ah + bh);   n = an ./ (an + bn);
w  = kv3_inf(V, p);
r  = h_inf(V, p);
Ca = p.Ca0 .* ones(1, nS);

N    = round(T / dt);
Vout = zeros(floor((N - 1) / rec) + 1, nS, 'single');
for k = 0:N-1
    t = k * dt;
    if t >= tOn && t < tOff, I = Iinj; else, I = 0 * Iinj; end

    [am, bm, ah, bh, an, bn] = naK_rates(V, p);
    winf = kv3_inf(V, p);
    tauw = p.tw0 + p.tw1 ./ (1 + exp((V + 15) / 12));
    sinf = 1 ./ (1 + exp(-(V - p.sV) ./ p.sk));
    rinf = h_inf(V, p);

    INa  = p.gNa  .* m.^3 .* h .* (V - p.ENa);
    IK   = p.gK   .* n.^4 .* (V - p.EK);
    IKv3 = p.gKv3 .* w .* (V - p.EK);
    ICa  = p.gCa  .* sinf.^2 .* (V - p.ECa);
    IAHP = p.gAHP .* Ca ./ (Ca + p.Kd) .* (V - p.EK);
    Ih   = p.gh   .* r .* (V - p.Eh);
    IL   = p.gL   .* (V - p.EL);

    V = V + dt * (I - INa - IK - IKv3 - ICa - IAHP - Ih - IL) ./ p.C;

    tm = 1 ./ (am + bm);  m = m + (am .* tm - m) .* (1 - exp(-dt ./ tm));
    th = 1 ./ (ah + bh);  h = h + (ah .* th - h) .* (1 - exp(-dt ./ th));
    tn = 1 ./ (an + bn);  n = n + (an .* tn - n) .* (1 - exp(-dt ./ tn));
    w  = w + (winf - w) .* (1 - exp(-dt ./ tauw));
    r  = r + (rinf - r) .* (1 - exp(-dt ./ p.taur));

    Ca = Ca + dt * (-p.phi .* ICa - (Ca - p.Ca0) ./ p.tauCa);

    if mod(k, rec) == 0
        Vout(k / rec + 1, :) = V;
    end
end

fs = 1000 / (dt * rec);
tv = (0:size(Vout, 1) - 1)' / fs;
sweeps = struct('name', {}, 't', {}, 'v', {}, 'i', {}, 'fs', {}, ...
                'amp', {}, 'tOn', {}, 'tOff', {});
for j = 1:nS
    iv = zeros(size(tv));  iv(tv >= tOn/1000 & tv < tOff/1000) = amps(j);
    sweeps(j) = struct('name', sprintf('Model_%d', j), 't', tv, ...
        'v', double(Vout(:, j)), 'i', iv, 'fs', fs, 'amp', amps(j), ...
        'tOn', tOn/1000, 'tOff', tOff/1000);
end
end

function [am, bm, ah, bh, an, bn] = naK_rates(V, p)
u  = V - p.VT;
am = 0.32  * vtrap(-(u - 13), 4);
bm = 0.28  * vtrap(u - 40, 5);
ah = 0.128 * exp(-(u - 17) / 18);
bh = 4 ./ (1 + exp(-(u - 40) / 5));
an = 0.032 * vtrap(-(u - 15), 5);
bn = 0.5   * exp(-(u - 10) / 40);
am = p.q .* am;  bm = p.q .* bm;  ah = p.q .* ah;
bh = p.q .* bh;  an = p.q .* an;  bn = p.q .* bn;
end

function y = vtrap(x, s)
y = x ./ (exp(x / s) - 1);
small = abs(x / s) < 1e-6;
y(small) = s * (1 - x(small) / s / 2);
end

function winf = kv3_inf(V, p)
winf = 1 ./ (1 + exp(-(V - p.wV) ./ p.wk));
end

function rinf = h_inf(V, p)
rinf = 1 ./ (1 + exp((V - p.rV) / 7));
end
