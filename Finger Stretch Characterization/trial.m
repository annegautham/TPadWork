clear; close all; clc;

% ================================
%  Load data
% ================================
load interpolatedImpedance_data.mat      % freq_interp, ohm_interp
load('28.4g_30-350.mat')                 % outputVec, inputVec, tvec

Fs = 10000;
Ts = 1/Fs;

% ================================
% 1) Generate chirp, compute current & force
% ================================
[voltage, t, current, instFreq] = chirpWithEnvelope( ...
        30, 350, 0.8, 5, Fs, freq_interp, ohm_interp, false);

BL = 2.7;
F_full = BL * current(:);    % N

figure; plot(t, F_full);
title('Input Force (theoretical)'); xlabel('Time (s)');

% ================================
% 2) LDV velocity → displacement
% ================================
raw = outputVec{1}(:);            % LDV channel
raw = raw - mean(raw);
raw = raw * 125;                  % mm/s

raw = lowpass(raw,1000,Fs);
b_hp = fir1(500,25/(Fs/2),'high');
raw = filter(b_hp,1,raw);

vel_mps_full = raw / 1000;        % m/s
pos_mm = cumsum(raw)*Ts;
pos_mm = pos_mm - mean(pos_mm);
pos_full = pos_mm / 1000;         % m

% ================================
% 3) ALIGN velocity to force (safe)
% ================================
[cc, lags] = xcorr(vel_mps_full, F_full, 'none');
[~, idxMax] = max(abs(cc));
lag = lags(idxMax);

if lag >= 0
    % velocity occurs AFTER force
    start_pos = lag + 1;
    start_force = 1;
else
    % force occurs AFTER velocity
    start_pos = 1;
    start_force = -lag + 1;
end

% compute maximum usable length
Lmax = min(length(vel_mps_full)-start_pos+1, ...
           length(F_full)-start_force+1);

vel_sync = vel_mps_full(start_pos:start_pos+Lmax-1);
pos_sync = pos_full(start_pos:start_pos+Lmax-1);
F_sync   = F_full(start_force:start_force+Lmax-1);

t_sync = (0:Lmax-1)'/Fs;

figure;
plot(t_sync, F_sync, t_sync, pos_sync);
legend('Force (N)','Displacement (m)');
title('Aligned Force & Displacement');

% ================================
% 4) Spectra and FRF
% ================================
Fspec = fft(F_sync);
Xspec = fft(pos_sync);
Vspec = fft(vel_sync);

N = length(Fspec);
f = (0:N-1)' * (Fs/N);

H = Xspec ./ Fspec;   % compliance FRF (m/N)

figure; plot(f, abs(H)); xlim([0 500]);
title('|H(f)| = X/F'); xlabel('Hz');

% ================================
% 5) Impedance Z = F / V
% ================================
Z = Fspec ./ Vspec;   % mechanical impedance

idx = f>30 & f<350;
fZ = f(idx); 
wZ = 2*pi*fZ;
Z_fit = Z(idx);

Zr = real(Z_fit);
Zi = imag(Z_fit);

figure;
subplot(2,1,1); plot(fZ,Zr); title('Real(Z)');
subplot(2,1,2); plot(fZ,Zi); title('Imag(Z)');

% ================================
% 6) Wiertlewski & Hayward (global parameters)
%     Z = b + j( m*w - k/w )
% ================================
b_est = mean(Zr);             % viscosity

A = [wZ, -1./wZ];
theta = A \ Zi;
m_est = theta(1);
k_est = theta(2);

% zero-crossing method
s = sign(Zi);
zc = find(diff(s)~=0,1);
if ~isempty(zc)
    w0 = interp1(Zi(zc:zc+1), wZ(zc:zc+1), 0);
    m_zero = k_est / w0^2;
else
    m_zero = m_est;
end

fprintf('\n===== Wiertlewski & Hayward Parameter Estimates =====\n');
fprintf('Viscosity b       = %.6f N·s/m\n', b_est);
fprintf('Stiffness k       = %.6f N/m   (%.3f N/mm)\n', k_est, k_est/1000);
fprintf('Mass (LSQ) m      = %.6f kg\n', m_est);
fprintf('Mass (zero-x) m0  = %.6f kg\n', m_zero);

% ================================
% 7) Frequency-dependent parameters
% ================================
b_dyn = Zr;                   % Real(Z)
m_dyn = Zi ./ wZ;             % Imag(Z)/w
k_dyn = -wZ .* Zi;            % -w * Imag(Z)

figure; plot(fZ,b_dyn); title('Dynamic Damping b(ω)');
xlabel('Hz'); ylabel('N·s/m');

figure; plot(fZ,m_dyn*1000); title('Dynamic Mass m(ω)');
xlabel('Hz'); ylabel('g');

figure; plot(fZ, k_dyn/1000);
title('Dynamic Stiffness k(ω)'); 
xlabel('Hz'); ylabel('N/mm');
