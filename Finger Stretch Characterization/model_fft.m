% --- Pre-existing IR data: 'ir', Fs, t defined as before ---

% STEP 1: Estimate dominant frequency from FFT of measured IR
L = length(ir);
H = fft(ir, 2^nextpow2(L));
f = Fs * (0:(length(H)/2 - 1)) / length(H);
[~, idx_peak] = max(abs(H(1:length(f))));
f_peak = f(idx_peak);
fprintf('Estimated peak frequency from FFT: %.2f Hz\n', f_peak);

% STEP 2: Use that to build better initial guess
wd_guess = 2*pi*f_peak;
p0 = [max(ir), 0.1 * wd_guess, wd_guess, 0];  % [A, damping*wn, wd, phase]

% STEP 3: Fit model
ir_model = @(p, t) p(1) * exp(-p(2)*t) .* sin(p(3)*t + p(4)) .* (t > 0);
cost = @(p) sum((ir(:) - ir_model(p, t(:))).^2);
p_fit = fminsearch(cost, p0);

% STEP 4: Plot and verify
ir_fit = ir_model(p_fit, t);

figure;
plot(t, ir, 'k', 'DisplayName', 'Measured IR'); hold on;
plot(t, ir_fit, 'r--', 'DisplayName', 'Fitted Model');
legend; title('Impulse Response Fit');

% STEP 5: Extract fitted parameters
A_fit = p_fit(1);
damping_wn = p_fit(2);
wd = p_fit(3);
phi = p_fit(4);
zeta = damping_wn / sqrt(wd^2 + damping_wn^2);
wn = damping_wn / zeta;

fprintf('\nFitted Parameters:\n');
fprintf('  A     = %.4f\n', A_fit);
fprintf('  ζ     = %.4f\n', zeta);
fprintf('  ω_n   = %.2f rad/s (%.2f Hz)\n', wn, wn/(2*pi));
fprintf('  ω_d   = %.2f rad/s (%.2f Hz)\n', wd, wd/(2*pi));

% STEP 6: Confirm FFT of model aligns
H_fit = fft(ir_fit, 2^nextpow2(length(ir_fit)));
f_fit = Fs * (0:(length(H_fit)/2 - 1)) / length(H_fit);
[~, idx_peak_model] = max(abs(H_fit(1:length(f_fit))));
fprintf('FFT peak of model IR: %.2f Hz\n', f_fit(idx_peak_model));

% Plot FFTs
figure;
plot(f, abs(H(1:length(f))), 'k', 'DisplayName', 'Measured IR'); hold on;
plot(f_fit, abs(H_fit(1:length(f_fit))), 'r--', 'DisplayName', 'Fitted Model');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Comparison: Measured vs Fitted Model');
legend; grid on;

% Physical parameter estimation
m = 18 / 1000;  % kg
k = m * wn^2;
b = 2 * zeta * sqrt(m * k);

fprintf('\nEstimated Physical Parameters:\n');
fprintf('  Mass (m)              = %.4f kg\n', m);
fprintf('  Spring constant (k)   = %.2f N/m\n', k);
fprintf('  Damping coefficient (b) = %.4f Ns/m\n', b);


