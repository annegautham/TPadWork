% Example signal (50 Hz sine wave with harmonic distortion)
[signal, t] = chirpWithEnvelope(30, 350, 2.4, 5, 10000);

% Add harmonic distortion with decreasing amplitude
base_amplitude = 0.09;                        % Amplitude of the first harmonic
harmonics = [2, 3, 4, 5, 6, 7, 8];            % Harmonic multiples
distorted_signal = signal;

for i = 1:length(harmonics)
    harmonic_amplitude = base_amplitude / i;   % Decreasing amplitude
    distorted_signal = distorted_signal + harmonic_amplitude * sin(2 * pi * harmonics(i) * 50 * t);
end

% Add noise around 171 Hz + other nearby random frequencies
signallr = fliplr(distorted_signal);
noise_level = 0.03;                           % Noise amplitude
frequencies = [171, 165, 180, 190];            % Slightly spread noise frequencies
noise = zeros(size(50000));

% Generate the mixed noise signal
for f = frequencies
    noise = noise + noise_level * sin(2 * pi * f * t + 2 * pi * rand); % Random phase
end

% ---- Shift signallr by 4 seconds to the right ----
fs = 1 / (t(2) - t(1));                      % Sampling frequency
shift_samples = round(4 * fs);                % 4-second shift in samples

% Pad with zeros before and after
signallr_shifted = [zeros(1, shift_samples), signallr, zeros(1, shift_samples)];

% Adjust the time vector to match the new length
t_shifted = (0:(length(signallr_shifted)-1)) / fs;

% Zero pad the noise to match the length of the shifted signal
noise_padded = [zeros(1, shift_samples), noise, zeros(1, length(signallr_shifted) - length(noise) - shift_samples)];

% Combine the shifted signal and padded noise
ns = signallr_shifted + noise_padded;

% Convolve the distorted signal with the noisy signal
c = conv(distorted_signal, ns);

plot(c);