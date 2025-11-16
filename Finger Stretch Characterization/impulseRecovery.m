clear; close all; clc;

% Load resistance curve and DAQ data
load interpolatedImpedance_data.mat      % freq_interp, ohm_interp
load('28.4g_30-350.mat')                 % outputVec, inputVec, tvec

Fs = 10000;
Ts = 1/Fs;

% ================== 1) Chirp, current, force ==================

[voltage, t, current, instFreq] = chirpWithEnvelope( ...
    30, 350, 0.8, 5, Fs, freq_interp, ohm_interp, false);

BL = 2.7;
F_full = BL * current(:);                % full 5 s force waveform

figure;
plot(t, F_full);
title('Force Input F(t) [N]');
xlabel('Time (s)');

% ================== 2) Velocity -> displacement (meters) =======

raw = outputVec{1}(:);                   % your motion channel
raw = raw - mean(raw);
raw = raw * 125;                         % mm/s

raw = lowpass(raw, 1000, Fs);
b = fir1(500, 25/(Fs/2), 'high');
raw = filter(b, 1, raw);

figure;
plot(raw);
title('Velocity (mm/s)');
xlabel('Sample');

pos_mm = cumsum(raw) * Ts;              % mm
pos_mm = pos_mm - mean(pos_mm);

figure;
plot(pos_mm);
title('Displacement (mm)');
xlabel('Sample');

pos_full = pos_mm / 1000;               % meters

% ================== 3) Align displacement to force =============

N_F = length(F_full);
N_pos = length(pos_full);

% Cross-correlate full signals
[cc, lags] = xcorr(pos_full, F_full, 'none');
[~, L] = max(abs(cc));
lag = lags(L);

% Use the FULL F_full, and cut a matching window from pos_full
if lag >= 0
    % force starts earlier in time than displacement
    start_idx = lag + 1;
    stop_idx  = start_idx + N_F - 1;
    if stop_idx > N_pos
        stop_idx = N_pos;
        N_F = stop_idx - start_idx + 1;
        F_full = F_full(1:N_F);
    end
    pos = pos_full(start_idx:stop_idx);
    F_time = F_full(1:length(pos));
else
    % displacement starts earlier than force
    shift = -lag;
    start_idx_F = shift + 1;
    stop_idx_F  = min(start_idx_F + N_pos - 1, N_F);
    F_time = F_full(start_idx_F:stop_idx_F);
    pos    = pos_full(1:length(F_time));
end

N = min(length(pos), length(F_time));
pos = pos(1:N);
F_time = F_time(1:N);

figure;
plot((0:N-1)/Fs, F_time, (0:N-1)/Fs, pos*1000);
legend('Force (N)', 'Displacement (m)');
xlabel('Time (s)');
title('Aligned Force Input and Output Displacement');

% ================== 4) FRF in frequency domain =================

Fspec = fft(F_time);
Xspec = fft(pos);
f = (0:N-1)' * (Fs/N);

figure;
plot(f, abs(Fspec)); xlim([0 500]);
xlabel('Frequency (Hz)');
ylabel('|F(f)|');
title('Force Spectrum');

figure;
plot(f, abs(Xspec)); xlim([0 500]);
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
title('Displacement Spectrum');

H = Xspec ./ Fspec;                      % FRF in m/N

figure;
plot(f, abs(H)); xlim([0 500]);
xlabel('Frequency (Hz)');
ylabel('|H(f)|');
title('FRF |H(f)| = |X(f)/F(f)| (m/N)');

% ================== 5) Dynamic extraction: m, c, k =============

idx = f > 20 & f < 350;
f_fit = f(idx);
w = 2*pi*f_fit;
H_fit = H(idx);

Y  = 1 ./ H_fit;                         % inverse FRF
Yr = real(Y);
Yi = imag(Y);

figure;
subplot(2,1,1);
plot(f_fit, Yr);
xlabel('Frequency (Hz)');
ylabel('Real(1/H)');
title('Real(1/H) = k - m \omega^2');

subplot(2,1,2);
plot(f_fit, Yi);
xlabel('Frequency (Hz)');
ylabel('Imag(1/H)');
title('Imag(1/H) = c \omega');

% Solve for k, m from Yr = k - m w^2
A  = [ones(length(w),1), -w.^2];
sol = A \ Yr;
k = sol(1);
m = sol(2);

% Solve for c from Yi = c w
c = Yi \ w;

m
c
k

% ================== 6) Compare measured vs model FRF ===========

Hmodel = 1 ./ (k - m*w.^2 + 1j*c*w);

figure;
subplot(2,1,1);
plot(w/(2*pi), 20*log10(abs(H_fit)), 'b', ...
     w/(2*pi), 20*log10(abs(Hmodel)), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Measured','Model');
title('FRF Magnitude Fit');

subplot(2,1,2);
plot(w/(2*pi), unwrap(angle(H_fit)), 'b', ...
     w/(2*pi), unwrap(angle(Hmodel)), 'r');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('Measured','Model');
title('FRF Phase Fit');

% ================== 7) Wiertlewski & Hayward Dynamic Extraction ==================
% Z(jw) = 1/(j*w*H(jw))
Z = 1 ./ (1j * w .* H_fit);   % impedance over fitting band

Zr = real(Z);
Zi = imag(Z);

% ---- 1) Viscosity b = average real(Z) ----
b_est = mean(Zr);

% ---- 2) Fit imaginary part: Im(Z) = m*w - k*(1/w) ----
X1 = w;           % coefficient for mass (m*w)
X2 = -1 ./ w;     % coefficient for stiffness (-k/w)

A = [X1, X2];
theta = A \ Zi;   % least-squares fit for [m; k]

m_est = theta(1);
k_est = theta(2);

% ---- 3) Zero-crossing refinement: Im(Z(w0)) = 0 ----
s = sign(Zi);
zc_idx = find(diff(s) ~= 0, 1);

if ~isempty(zc_idx)
    w0 = interp1(Zi(zc_idx:zc_idx+1), w(zc_idx:zc_idx+1), 0);
    m_refined = k_est / (w0^2);
else
    m_refined = m_est;   % fallback
end

fprintf('\n===== Wiertlewski & Hayward Parameter Estimates =====\n');
fprintf('Viscosity b       = %.6f N·s/m\n', b_est);
fprintf('Stiffness k       = %.6f N/m\n',   k_est);
fprintf('Mass (LSQ) m      = %.6f kg\n',    m_est);
fprintf('Mass (zero-x) m0  = %.6f kg\n',    m_refined);

% ---- Plots: real and imaginary impedance ----
figure;
subplot(2,1,1);
plot(f_fit, Zr);
xlabel('Frequency (Hz)');
ylabel('Real(Z)');
title('Real Impedance (Viscosity Term)');

subplot(2,1,2);
plot(f_fit, Zi);
xlabel('Frequency (Hz)');
ylabel('Imag(Z)');
title('Imaginary Impedance (Mass–Stiffness Term)');



