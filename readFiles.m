function sweeps = readFiles(nwbFile)
% READFILES  Read all "Long Square" sweeps from an Allen Cell Types NWB file
% (NWB 1.x / HDF5) using only built-in MATLAB HDF5 functions.
%
%   sw = readFiles('329549955_ephys.nwb');
%
% Each element of sweeps has:
%   .name      sweep group name, e.g. 'Sweep_42'
%   .t         time (s)
%   .v         membrane potential (mV)
%   .i         injected current (pA)
%   .fs        sampling rate (Hz)
%   .amp       step amplitude (pA)
%   .tOn/.tOff step onset and offset (s)
%
% Units are decided from the size of the data (some older Allen files have
% wrong 'conversion' attributes), and blank (NaN) samples are trimmed off.

info  = h5info(nwbFile, '/acquisition/timeseries');
names = {info.Groups.Name};

sweeps = struct('name', {}, 't', {}, 'v', {}, 'i', {}, 'fs', {}, ...
                'amp', {}, 'tOn', {}, 'tOff', {});

for k = 1:numel(names)
    sw = names{k}(find(names{k} == '/', 1, 'last') + 1:end);

    stimName = readString(nwbFile, ...
        {['/stimulus/presentation/' sw '/aibs_stimulus_name'], ...
         ['/acquisition/timeseries/' sw '/aibs_stimulus_name']});
    if ~contains(lower(stimName), 'long square')
        continue
    end

    v = readSI(nwbFile, ['/acquisition/timeseries/' sw '/data'], 'volt');
    i = readSI(nwbFile, ['/stimulus/presentation/' sw '/data'],  'amp');
    v = v * 1e3;
    i = i * 1e12;

    n    = min(numel(v), numel(i));
    v    = v(1:n);  i = i(1:n);
    good = ~isnan(v) & ~isnan(i);
    last = find(good, 1, 'last');
    v    = v(1:last);  i = i(1:last);
    v(isnan(v)) = 0;  i(isnan(i)) = 0;

    fs = double(h5readatt(nwbFile, ...
        ['/acquisition/timeseries/' sw '/starting_time'], 'rate'));
    t  = (0:numel(v)-1)' / fs;

    base = median(i(1:round(0.005 * fs)));
    on   = abs(i - base) > 1;
    [s, e] = longestRun(on);
    if isempty(s)
        tOn = NaN; tOff = NaN; ampDetected = 0;
    else
        tOn = t(s); tOff = t(e); ampDetected = median(i(s:e)) - base;
    end

    amp = readNumber(nwbFile, ['/stimulus/presentation/' sw '/aibs_stimulus_amplitude_pa']);
    if isnan(amp)
        amp = ampDetected;
    end

    sweeps(end+1) = struct('name', sw, 't', t, 'v', v, 'i', i, 'fs', fs, ...
                           'amp', round(amp), 'tOn', tOn, 'tOff', tOff);
end

[~, idx] = sort([sweeps.amp]);
sweeps   = sweeps(idx);
fprintf('Loaded %d long square sweeps (%d to %d pA)\n', numel(sweeps), ...
        sweeps(1).amp, sweeps(end).amp);
end

function d = readSI(f, p, expected)
d = double(h5read(f, p));
d = d(:);
if strcmp(expected, 'volt')
    if median(abs(d), 'omitnan') > 1
        d = d / 1e3;
    end
else
    if max(abs(d), [], 'omitnan') > 1e-3
        d = d / 1e12;
    end
end
end

function s = readString(f, paths)
s = '';
for p = paths
    try
        x = h5read(f, p{1});
        if iscell(x), x = x{1}; end
        s = strtrim(char(x(:)'));
        return
    catch
    end
end
end

function x = readNumber(f, p)
try
    x = double(h5read(f, p));
catch
    x = NaN;
end
end

function [s, e] = longestRun(mask)
mask = mask(:)';
d = diff([0 mask 0]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
if isempty(starts), s = []; e = []; return, end
[~, k] = max(ends - starts);
s = starts(k);  e = ends(k);
end
