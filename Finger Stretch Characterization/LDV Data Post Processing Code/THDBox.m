% Parameters
Fs = 10000; % Sampling frequency
Fc = 1000;  % Lowpass cutoff frequency
windowLength = 1024; % Spectrogram window length
overlap = 768; % Overlap for spectrogram
fftPoints = 2048; % FFT points for resolution

% Generate Spectrogram
[s, f, t] = spectrogram(signal, windowLength, overlap, fftPoints, Fs);

% Frequency Range for THD Analysis
freqCenters = 10:20:350; % Frequency centers (10, 30, ..., 350 Hz)
bandWidth = 10; % Half-width of each frequency band

% Initialize THD Results
thdValues = [];

% Compute THD for Each Frequency Band
for fc = freqCenters
    % Define frequency range for the current band
    bandIdx = (f >= (fc - bandWidth)) & (f <= (fc + bandWidth));
    
    % Extract power in the frequency band
    bandPower = abs(s(bandIdx, :)).^2; % Power spectrum
    [maxPower, idxMax] = max(bandPower, [], 1); % Fundamental power
    harmonicPower = sum(bandPower, 1) - maxPower; % Harmonic power
    
    % THD Calculation in dB
    thd_dB = 10 * log10(harmonicPower ./ maxPower); 
    thdValues = [thdValues; thd_dB];
end

% Box Plot of THD in dB Across Frequency Bands
figure;
boxplot(thdValues', freqCenters);
xlabel('Frequency (Hz)');
ylabel('THD (dB)');
title('THD Across Frequency Bands');
grid on;
