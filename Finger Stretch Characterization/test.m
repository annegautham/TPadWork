clear; close all; clc;

% ----------------- LOAD DATA -----------------
load interpolatedImpedance_data.mat      % freq_interp, ohm_interp
load('28.4g_30-350.mat')                 % outputVec, inputVec, tvec

Fs = 10000;
Ts = 1/Fs;

% ================== 1) Chirp, current, force ==================
[voltage, t, current, instFreq] = chirpWithEnvelope( ...
    30, 350, 0.8, 5, Fs, freq_interp, ohm_interp, false);

BL = 2.7;
F_full = BL * current(:);        % theoretical force [N]

figure; plot(t, F_full);
title('Force Input F(t) [N]'); xlabel('Time (s)');

% ================== 2) Velocity -> displacement ================
raw = outputVec{1}(:);
raw = raw - mean(raw);
raw = raw * 125;                 % mm/s (LDV calib)

raw = lowpass(raw, 1000, Fs);
b_hp = fir1(500, 25/(Fs/2), 'high');
raw = filter(b_hp, 1, raw);      % still mm/s

figure; plot(raw); title('Velocity (mm/s)');

pos_mm = cumsum(raw)*Ts;
pos_mm = pos_mm - mean(pos_mm);

figure; plot(pos_mm); title('Displacement (mm)');

pos_full     = pos_mm / 1000;    % m
vel_mps_full = raw    / 1000;    % m/s

% ================== 3) Alignment (xcorr) =======================
N_F   = length(F_full);
N_pos = length(pos_full);

[cc,lags] = xcorr(pos_full, F_full,'none');
[~,L] = max(abs(cc));
lag = lags(L);

if lag >= 0
    start_idx = lag + 1;
    stop_idx  = start_idx + N_F - 1;
    if stop_idx > N_pos
        stop_idx = N_pos;
        N_F = stop_idx - start_idx + 1;
        F_full = F_full(1:N_F);
    end
    pos   = pos_full(start_idx:stop_idx);
    vel   = vel_mps_full(start_idx:stop_idx);
    F_time = F_full(1:length(pos));
else
    shift      = -lag;
    start_idx_F = shift + 1;
    stop_idx_F  = min(start_idx_F + N_pos - 1, N_F);
    F_time = F_full(start_idx_F:stop_idx_F);
    pos    = pos_full(1:length(F_time));
    vel    = vel_mps_full(1:length(F_time));
end

N = min([length(pos), length(F_time), length(vel)]);
pos   = pos(1:N);
vel   = vel(1:N);
F_time = F_time(1:N);

time = (0:N-1)'/Fs;
figure; plot(time, F_time, time, pos);
legend('Force (N)', 'Displacement (m)');
title('Aligned Force Input and Output Displacement');
xlabel('Time (s)');

% ================== 4) Spectra & FRF ===========================
Fspec = fft(F_time);
Xspec = fft(pos);
Vspec = fft(vel);
f = (0:N-1)' * (Fs/N);

H = Xspec ./ Fspec;      % compliance (m/N), mostly sanity check

figure; plot(f, abs(H)); xlim([0 500]);
xlabel('Hz'); ylabel('|H| (m/N)');
title('FRF |H(f)| = X/F');

% ================== 5) Impedance Z = F / V =====================
Z = Fspec ./ Vspec;      % mechanical impedance [N·s/m]

% keep 30–350 Hz band
idx = (f > 30 & f < 350);
fZ = f(idx);
wZ = 2*pi*fZ;
Z_fit = Z(idx);

Zr = real(Z_fit);
Zi = imag(Z_fit);

% ---- enforce sign convention: low-frequency Real(Z) > 0 ----
lowBand = (fZ > 35 & fZ < 120);
if mean(Zr(lowBand)) < 0
    Z_fit = -Z_fit;
    Zr = -Zr;
    Zi = -Zi;
    fprintf('Flipped Z sign so low-frequency Real(Z) > 0.\n');
end

figure;
subplot(2,1,1);
plot(fZ, Zr); xlabel('Hz'); ylabel('Real(Z)');
title('Real(Z)');

subplot(2,1,2);
plot(fZ, Zi); xlabel('Hz'); ylabel('Imag(Z)');
title('Imag(Z)');

% ================== 6) W&H global parameter extraction =========
% Z(jw) = b + j(m w - k / w)

% ---- 6a) viscosity b from real part, away from resonance ----
band_b = (fZ > 40  & fZ < 140) | (fZ > 230 & fZ < 320);
b_est = mean(Zr(band_b));

% ---- 6b) mass & stiffness from imaginary part around resonance ----
% find resonance roughly from max |Imag(Z)|
[~,i_res] = max(abs(Zi));
f0  = fZ(i_res);
w0  = 2*pi*f0;

band_mk = (fZ > 0.7*f0) & (fZ < 1.3*f0);   % narrow around resonance
w_mk  = 2*pi*fZ(band_mk);
Zi_mk = Zi(band_mk);

A = [w_mk, -1./w_mk];   % Imag(Z) = m*w - k/w
theta = A \ Zi_mk;
m_est = theta(1);
k_est = theta(2);

% ---- 6c) refined mass from zero crossing of Imag(Z) ----
s = sign(Zi);
zc_idx = find(diff(s) ~= 0, 1);
if ~isempty(zc_idx)
    w0_zero = interp1(Zi(zc_idx:zc_idx+1), wZ(zc_idx:zc_idx+1), 0);
    m_zero  = k_est / (w0_zero^2);
else
    m_zero = m_est;
end

fprintf('\n===== Wiertlewski & Hayward Parameter Estimates (approx) =====\n');
fprintf('Viscosity b       = %.6f N·s/m\n',        b_est);
fprintf('Stiffness k       = %.6f N/m   (%.3f N/mm)\n', k_est, k_est/1000);
fprintf('Mass (LSQ) m      = %.6f kg\n',          m_est);
fprintf('Mass (zero-x) m0  = %.6f kg\n\n',        m_zero);

% ==============================================================
% === 7) Frequency-Dependent Dynamic Parameters (approx) =======
% ==============================================================

% Using Imag(Z) = m w - k / w, derive:
% m_dyn(ω) = (Imag(Z) + k_est/ω) / ω  (should hover near m)
% k_dyn(ω) = m_zero*ω^2 - ω*Imag(Z)   (should hover near k)

m_dyn = (Zi + k_est./wZ) ./ wZ;       % kg
k_dyn =  m_zero*(wZ.^2) - wZ.*Zi;     % N/m
b_dyn = Zr;                           % N·s/m

figure;
plot(fZ, b_dyn, 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('b_{dyn} (N·s/m)');
title('Approx. Dynamic Damping b(\omega)');
grid on;

figure;
plot(fZ, m_dyn*1000, 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('m_{dyn} (g)');
title('Approx. Dynamic Mass m(\omega)');
grid on;

figure;
plot(fZ, k_dyn/1000, 'LineWidth', 1.2); hold on;
yline(k_est/1000, '--r', 'k_{global}', 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('k_{dyn} (N/mm)');
title('Approx. Dynamic Stiffness k(\omega)');
grid on;
