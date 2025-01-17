% Example parameters
fs = 5000; % Sampling frequency in Hz
t = 0:300000; % Time vector for 10 seconds
signal = useTemp;
reference_waveform =signal;
% Normalize the signal and reference waveform
signal = signal / max(abs(signal));
reference_waveform = reference_waveform / max(abs(reference_waveform));

% Compute cross-correlation
[xcorr_result, lags] = xcorr(signal, reference_waveform, 'normalized');

% Find peaks in the cross-correlation
threshold = 0.8; % Adjust based on desired correlation strength
[peaks, locs] = findpeaks(xcorr_result, 'MinPeakHeight', threshold);

% Convert lags to timestamps
detected_lags = lags(locs);
timestamps = detected_lags / fs; % Convert to seconds

% Display results
disp('Timestamps of detected correlations:');
disp(timestamps);

% Optional: Plot results
figure;
subplot(2,1,1);
plot(signal);
title('Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(lags/fs, xcorr_result);
hold on;
plot(timestamps, peaks, 'ro'); % Mark detected peaks
title('Cross-Correlation');
xlabel('Lag (s)');
ylabel('Correlation');
