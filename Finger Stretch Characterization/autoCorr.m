h = ir_bef;
% Compute autocorrelation of the impulse response with itself
[acf, lags] = xcorr(h, 'coeff');  % Autocorrelation of IR

% Find the index of the maximum value of the autocorrelation (which is at lag 0)
[~, center] = max(abs(acf));

% Find peaks in the autocorrelation after the central lag (for harmonics)
[peakVals, locs] = findpeaks(acf(center+1:end), 'MinPeakHeight', 0.05, 'MinPeakDistance', 1000);

% Calculate harmonic lags based on detected peaks (after the central one)
harmonic_lags_samples = lags(center + locs);

% Display harmonic lags
disp('Harmonic lags (in samples):');
disp(harmonic_lags_samples);
% Find the main peak of the impulse response (linear response)
[~, peakIdx] = max(abs(h));  % The index of the maximum value in the IR

% Use the lag between the linear response and the first harmonic
lag_linear_to_harmonic = harmonic_lags_samples(1);  % The first harmonic lag

% Define the range around the peak where the linear response is located
halfWindow = lag_linear_to_harmonic;  % +/- lag from the peak
idxStart = max(1, peakIdx - halfWindow);
idxEnd   = min(length(h), peakIdx + halfWindow);

% Extract the linear response (centered around the peak and +/- the lag)
linear_response = h(idxStart:idxEnd);

% Plot the linear response
figure;
plot(linear_response);
title('Extracted Linear Response from Impulse Response');
xlabel('Samples');
ylabel('Amplitude');
