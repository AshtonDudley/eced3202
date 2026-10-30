%% THD of amplifier input/output from a Siglent SDS1104X-E CSV export
% CH1 = amplifier input, CH2 = amplifier output
% Method: block-average to reduce data size, estimate the fundamental
% frequency with a sine fit, then least-squares fit DC + H harmonics.
% The fit does not require an integer number of cycles in the record.

clear; clc;

fname = 'thd_data.csv';   % path to the exported CSV
nHarm = 10;               % harmonics 1..nHarm used in the fit and THD
decim = 1000;             % block-average factor (500 MSa/s -> 500 kSa/s)
nHead = 12;               % 11 metadata lines + 1 column header line

%% Read header and data
fid = fopen(fname, 'r');
assert(fid > 0, 'Could not open %s', fname);

dt = [];
for k = 1:nHead
    ln = fgetl(fid);
    if startsWith(ln, 'Sample Interval')
        tok = regexp(ln, 'Analog:([-+0-9.Ee]+)', 'tokens', 'once');
        dt  = str2double(tok{1});
    end
end

% textscan continues from the current file position (after the header)
C = textscan(fid, '%f %f %f', 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);
D = C{1};
clear C

if isempty(dt)
    dt = (D(end,1) - D(1,1)) / (size(D,1) - 1);
end
fprintf('Loaded %d samples, fs = %.3g Sa/s\n', size(D,1), 1/dt);

%% Block-average (acts as a low-pass filter, reduces quantization noise)
n    = floor(size(D,1) / decim) * decim;
blk  = @(v) mean(reshape(v(1:n), decim, []), 1).';
vin  = blk(D(:,2));
vout = blk(D(:,3));
fs   = 1 / (dt * decim);
clear D

%% Fundamental frequency (estimated from the output, used for both)
f0 = estimateF0(vout, fs);
fprintf('Fundamental: %.4f Hz (%.2f cycles in record)\n', ...
        f0, numel(vout) * f0 / fs);

%% Harmonic fit and THD
[Ain,  phIn,  ~] = harmonicFit(vin,  fs, f0, nHarm);
[Aout, phOut, ~] = harmonicFit(vout, fs, f0, nHarm);

% Fundamental-only fit gives the residual used for THD+N
[~, ~, resIn]  = harmonicFit(vin,  fs, f0, 1);
[~, ~, resOut] = harmonicFit(vout, fs, f0, 1);

thdIn   = sqrt(sum(Ain(2:end).^2))  / Ain(1);
thdOut  = sqrt(sum(Aout(2:end).^2)) / Aout(1);
thdnIn  = sqrt(mean(resIn.^2))  / (Ain(1)  / sqrt(2));
thdnOut = sqrt(mean(resOut.^2)) / (Aout(1) / sqrt(2));

fprintf('\n%-8s %12s %12s %12s %12s\n', '', 'Fund (Vpk)', 'THD (%)', 'THD (dB)', 'THD+N (%)');
fprintf('%-8s %12.4f %12.4f %12.2f %12.4f\n', 'Input',  Ain(1),  100*thdIn,  20*log10(thdIn),  100*thdnIn);
fprintf('%-8s %12.4f %12.4f %12.2f %12.4f\n', 'Output', Aout(1), 100*thdOut, 20*log10(thdOut), 100*thdnOut);

gain = Aout(1) / Ain(1);
dph  = rad2deg(angle(exp(1i * (phOut(1) - phIn(1)))));
fprintf('\nGain: %.3f V/V (%.2f dB), phase: %.1f deg\n', gain, 20*log10(gain), dph);

fprintf('\nHarmonic levels relative to fundamental (dBc):\n');
fprintf('%4s %10s %10s\n', 'h', 'Input', 'Output');
for h = 2:nHarm
    fprintf('%4d %10.1f %10.1f\n', h, ...
            20*log10(Ain(h)/Ain(1)), 20*log10(Aout(h)/Aout(1)));
end

%% Spectrum plot (trimmed to whole cycles, rectangular window)
figure('Name', 'Amplifier spectrum');
hold on; grid on;
plotSpectrum(vin,  fs, f0, 'Input (CH1)');
plotSpectrum(vout, fs, f0, 'Output (CH2)');
xlim([0, (nHarm + 1) * f0]);
xlabel('Frequency (Hz)'); ylabel('Level (dBV rms)');
legend show;
title(sprintf('f_0 = %.1f Hz, THD_{out} = %.3f %%', f0, 100*thdOut));

%% Local functions
function f0 = estimateF0(x, fs)
    x  = x - mean(x);
    N  = numel(x);
    w  = 0.5 - 0.5*cos(2*pi*(0:N-1).'/(N-1));
    Nf = 2^nextpow2(16*N);
    X  = abs(fft(x .* w, Nf));
    [~, k] = max(X(2:floor(Nf/2)));
    fc = k * fs / Nf;                       % coarse estimate
    df = fs / N;
    t  = (0:N-1).' / fs;
    f0 = fminbnd(@(f) sineResid(x, t, f), max(fc - df, df/10), fc + df);
end

function r = sineResid(x, t, f)
    M = [cos(2*pi*f*t), sin(2*pi*f*t), ones(size(t))];
    r = norm(x - M * (M \ x));
end

function [A, ph, res] = harmonicFit(x, fs, f0, H)
    t = (0:numel(x)-1).' / fs;
    M = ones(numel(x), 2*H + 1);
    for h = 1:H
        M(:, 2*h)     = cos(2*pi*h*f0*t);
        M(:, 2*h + 1) = sin(2*pi*h*f0*t);
    end
    c   = M \ x;
    a   = c(2:2:end);
    b   = c(3:2:end);
    A   = hypot(a, b);       % peak amplitude per harmonic
    ph  = atan2(-b, a);      % phase per harmonic (rad)
    res = x - M * c;
end

function plotSpectrum(x, fs, f0, name)
    nCyc = floor(numel(x) * f0 / fs);
    M    = round(nCyc * fs / f0);
    seg  = x(1:M) - mean(x(1:M));
    S    = 2 * abs(fft(seg)) / M;             % peak amplitude
    f    = (0:M-1).' * fs / M;
    idx  = 1:floor(M/2);
    plot(f(idx), 20*log10(S(idx)/sqrt(2) + eps), '.-', 'DisplayName', name);
end
