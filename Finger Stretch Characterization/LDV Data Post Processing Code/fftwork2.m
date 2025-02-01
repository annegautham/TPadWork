clear;
close all;

% Simulation parameters
Fs = 5000; % Sampling frequency
T = 1 / Fs;
N = 2^14; % Padding for FFT
f = Fs * (0:(floor(N / 2))) / N; % Frequency vector for plotting
resonanceFreq = 185; % Resonance frequency in Hz
chirpFreqStart = 30; % Chirp start frequency
chirpFreqEnd = 350; % Chirp end frequency
chirpAmplitude = 1; % Chirp amplitude

% Generate the chirp signal using signalGen's chirpWithEnvelope function
duration = 5; % Duration of the signal in seconds
% Call the chirpWithEnvelope function and capture the signal and time
[funs] = signalGen(); % Get the signal generation functions
[chirpSignal, t] = funs.chirpWithEnvelope(chirpFreqStart, chirpFreqEnd, chirpAmplitude, duration, Fs);

% Simulate a signal with a resonance pole at 185 Hz
resonanceSignal = simulateResonanceSignal(chirpSignal, resonanceFreq, Fs);

% Plot input and output signal
figure;
subplot(2, 1, 1);
plot(t, chirpSignal);
title('Input Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(t, resonanceSignal);
title('Resonance Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% FFT analysis for input (chirpSignal)
y_chirp = fft(chirpSignal, N);
p2_chirp = abs(y_chirp / N);
p1_chirp = p2_chirp(1:floor(N / 2) + 1);
p1_chirp(2:end-1) = 2 * p1_chirp(2:end-1);

% FFT analysis for output (resonanceSignal)
y_resonance = fft(resonanceSignal, N);
p2_resonance = abs(y_resonance / N);
p1_resonance = p2_resonance(1:floor(N / 2) + 1);
p1_resonance(2:end-1) = 2 * p1_resonance(2:end-1);

% Plot Power Spectrum for Chirp Signal
figure;
plot(f, p1_chirp.^2);
title('Power Spectrum of Chirp Signal');
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0 1000]);

%%%Force through voice coil

% Plot Power Spectrum for Resonance Signal
figure;
plot(f, p1_resonance.^2);
title('Power Spectrum of Resonance Signal');
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0 1000]);

% Calculate the Transfer Function (TF)
H = p1_resonance ./ p1_chirp; % Simple ratio for transfer function magnitude
H_phase = angle(y_resonance) - angle(y_chirp); % Phase difference

% Plot Transfer Function Magnitude
figure;
plot(f, H);
title('Magnitude of Transfer Function');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 1000]);

% Simulate Resonance Signal with Pole at 185 Hz
function outputSignal = simulateResonanceSignal(inputSignal, resonanceFreq, fs)
    % Transfer function for a second-order system with a resonance at 185 Hz
    % H(s) = 1 / (s^2 + 2 * zeta * w_n * s + w_n^2)
    % where w_n = 2*pi*resonanceFreq, zeta = 0.1 (damping factor)
    zeta = 0.1;
    wn = 2 * pi * resonanceFreq;
    s = tf('s');
    H = 1 / (s^2 + 2 * zeta * wn * s + wn^2);
    outputSignal = lsim(H, inputSignal, (0:length(inputSignal)-1) / fs);
end
