% Signal parameters
fs = 10000; % Sampling frequency (Hz)
t = 0:1/fs:1; % Time vector
f0 = 500; % Fundamental frequency (Hz)
signal = sin(2*pi*f0*t) + 0.3*sin(2*pi*2*f0*t) + 0.2*sin(2*pi*3*f0*t); % Signal with harmonics

% Spectrogram
window = hamming(256);
noverlap = 128;
nfft = 1024;
[s, f, t, p] = spectrogram(signal, window, noverlap, nfft, fs);

% Convert to magnitude spectrum (averaging over time)
avg_power = mean(p, 2); % Average power over time
magnitude = sqrt(avg_power); % Convert to amplitude
frequencies = f; % Frequency vector

% Find fundamental frequency and harmonics
f1_idx = find(frequencies >= f0, 1); % Index of fundamental frequency
A1 = magnitude(f1_idx); % Amplitude of fundamental

% Find harmonic indices
harmonic_indices = f1_idx * (2:5); % Assume up to the 5th harmonic
harmonic_indices = harmonic_indices(harmonic_indices <= length(magnitude)); % Ensure valid indices

% Extract harmonic amplitudes
harmonic_amplitudes = magnitude(harmonic_indices);

% Calculate THD
THD = sqrt(sum(harmonic_amplitudes.^2)) / A1 * 100;

% Display results
disp(['Fundamental Amplitude (A1): ', num2str(A1)]);
disp(['Harmonic Amplitudes: ', num2str(harmonic_amplitudes')]);
disp(['Total Harmonic Distortion (THD): ', num2str(THD), '%']);

figure;
imagesc(t, f, 10*log10(p)); % Plot spectrogram
axis xy;
colormap jet;
colorbar;
hold on;
plot(t, f0*ones(size(t)), 'r--', 'LineWidth', 1.5); % Fundamental
for k = 2:5
    plot(t, k*f0*ones(size(t)), 'g--', 'LineWidth', 1.5); % Harmonics
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram with Harmonics');
