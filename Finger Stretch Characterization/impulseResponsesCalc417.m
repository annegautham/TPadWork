close all;
% load('0g_30-350.mat');
% load('18g_30-350.mat');
% load('24.4g_30-350.mat');
% load('27g_30-350.mat');
load('28.4g_30-350.mat');
Fs = 10000;
Ts = 1/Fs;
[in, tin] = chirpWithEnvelope(30, 350, 0.8, 5, 10000);



raw = outputVec{1};
raw = raw(:)';
vel = (raw - mean(raw))*125;
vel = lowpass(vel, 1000, Fs);

%high pass
order = 500;
cuttoffFreq = 25/(5000/2);
b = fir1(order, cuttoffFreq, 'high');
vel = filter(b, 1, vel);

%position
figure;
pos = cumsum(vel)*Ts;

pos = pos - mean(pos);
plot(pos)


[cc, lags] = xcorr(pos, in, 'none');
[~, max_idx] = max(abs(cc));
lag = lags(max_idx);  % lag in samples
if lag >= 0
    pos = pos(lag+1:end);
end
N_sync = min(length(in), length(pos));
in  = in(1:N_sync);
pos = pos(1:N_sync);


figure;
hold on;
plot(in*0.5);
plot(pos);

%IR recover
inv = fliplr(in);
ir = conv(inv, pos);
ir_bef = ir;
ir = ir(50000-2500:50000+2500);
figure
hold on;
plot(ir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lin = conv(in, ir);
sf = max(pos)/max(lin);
lin = sf*lin;
lin = lin(2500:end);


figure;
plot(pos);
hold on;
plot(lin);
legend('pos','lin');

figure;
spec(ir, Fs);
figure;
spec(pos,Fs);
figure;
spec(lin,Fs);


% Define time vector for IR
t = (-2500:2500) * Ts;
% FFT of the impulse response
L = length(ir);
H = fft(ir, 2^nextpow2(L));  % Zero-padding for better frequency resolution
f = Fs * (0:(length(ir)/2)) / length(ir);  % Frequency axis (up to Nyquist)


% Define model function
ir_model = @(p, t) p(1) * exp(-p(2)*t) .* sin(p(3)*t + p(4)) .* (t > 0);

% Initial guess: [A, damping*wn, wd, phase]
[~, idx_peak] = max(abs(fft(ir)));
f_guess = f(idx_peak);
p0 = [max(ir), 100, 2*pi*f_guess, 0];


% Cost function
cost = @(p) sum((ir(:) - ir_model(p, t(:))).^2);

% Fit using fminsearch
p_fit = fminsearch(cost, p0);

% Plot fit
figure;
plot(t, ir, 'k', 'DisplayName', 'Measured IR'); hold on;
plot(t, ir_model(p_fit, t), 'r--', 'DisplayName', 'Fitted Model');
legend;
title('Impulse Response Fit');

% Extract fitted parameters
A_fit = p_fit(1);
damping_wn = p_fit(2);
wd = p_fit(3);
phi = p_fit(4);
zeta = damping_wn / sqrt(wd^2 + damping_wn^2);
wn = damping_wn / zeta;

fprintf('Fitted parameters:\n');
fprintf('  A     = %.4f\n', A_fit);
fprintf('  ζ     = %.4f\n', zeta);
fprintf('  ω_n   = %.2f rad/s\n', wn);
fprintf('  f_n   = %.2f Hz\n', wn/(2*pi));
fprintf('  ω_d   = %.2f rad/s\n', wd);



% Calculate mass (m), spring constant (k), and damping coefficient (b)
% Assuming typical system response characteristics
% If you have an estimate for one of the parameters, use it, or assume reasonable values
m = 28.4/1000;  % Example mass (kg), you can change this based on your system
k = m * wn^2;   % Spring constant (N/m)
b = 2 * zeta * sqrt(m * k);  % Damping coefficient (Ns/m)

% Display results
fprintf('Physical Parameters from fit:\n');
fprintf('  Mass (m)       = %.2f kg\n', m);
fprintf('  Spring constant (k) = %.2f N/m\n', k);
fprintf('  Damping coefficient (b) = %.2f Ns/m\n', b);




% Plot magnitude spectrum
figure;
plot(f, abs(H(1:length(f))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of Impulse Response');
grid on;

% Compute fitted IR
ir_fit = ir_model(p_fit, t);

% Only consider t > 0
idx_pos = t > 0;
ir_trimmed = ir(idx_pos);
ir_fit_trimmed = ir_fit(idx_pos);

% Compute correlation coefficient on t > 0
R = corrcoef(ir_trimmed, ir_fit_trimmed);
corr_val = R(1, 2);

% Display result
fprintf('Correlation coefficient between measured and fitted IR (t > 0): %.4f\n', corr_val);

% Optional: Annotate correlation on plot
figure;
plot(t, ir, 'k', 'DisplayName', 'Measured IR'); hold on;
plot(t, ir_fit, 'r--', 'DisplayName', 'Fitted Model');
legend;
title(sprintf('Impulse Response Fit (Corr = %.4f)', corr_val));
xlabel('Time (s)');
ylabel('Amplitude');
