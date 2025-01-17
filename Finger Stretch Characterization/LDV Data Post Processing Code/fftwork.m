close all;

dataNF1 = load("LDV Data\chirp30-340_NF1.mat");
dataNF2 = load("LDV Data\chirp30-340_NF2.mat");
dataNF3 = load("LDV Data\chirp30-340_NF3.mat");
dataF1 = load("LDV Data\chirp30-340_F1.mat");
dataF2 = load("LDV Data\chirp30-340_F2.mat");
dataF3 = load("LDV Data\chirp30-340_F3.mat");
dataBox = load("LDV Data\chirp30-340_Box.mat");

rangesNF1 = {22181:38248, 47559:63730, 63750:79648, 259547:275641, 275641:291558};
rangesNF2 = {20700:36288, 62239:78065, 86415:103167,103167:119047, 128481:144726, 144726:160340, 173287:190084, 190084:205850, 230417:246542, 269287:284896};
rangesNF3 = {22761:38662, 46503:63606, 99837:115870, 141154:157099, 164834:181666, 181666:197500, 205996:222812,222812:238710, 238710:254834, 262562:279377, 279377:295127};
rangesBox = {11642:27576, 27576:43400,43400:59482, 59482:75302, 75302:91224, 91224:107939, 107939:122821,122821:141012, 141012:155604, 155604:170680, 170680:186353, 186353:202880};
% rangesF1 = {};
% rangesF2 = {};
% rangesF3 = {};

tempNF1 = dataNF1.temp;
tempNF2 = dataNF2.temp;
tempNF3 = dataNF3.temp;
tempF1 = dataF1.temp;
tempF2 = dataF2.temp;
tempF3 = dataF3.temp;
tempBox = dataBox.temp;

sf = 125;

curr = 2;
Fs = 5000; % sampling freq
T = 1 / Fs;
N = 2^14; %padding
powerAccum = zeros(floor(N / 2) + 1, 1);
fftRealAccum = zeros(N, 1);
fftImagAccum = zeros(N, 1);
f = Fs * (0:(floor(N / 2))) / N;


tempNF1 = lowpass(tempNF1, 1000, Fs);
tempNF2 = lowpass(tempNF2, 1000, Fs);
tempNF3 = lowpass(tempNF2, 1000, Fs);
tempF1 = lowpass(tempF1, 1000, Fs);
tempF2 = lowpass(tempF2, 1000, Fs);
tempF3 = lowpass(tempF2, 1000, Fs);
tempBox = lowpass(tempBox, 1000, Fs);

useRanges = rangesNF2;
numRange = length(useRanges);
useTemp = tempNF2;


for i = 1:numRange
    rangeSignal = useTemp(useRanges{i})*sf;
    rangeSignal = rangeSignal - mean(rangeSignal);
    % truncate or add zero pad to signal
    if length(rangeSignal) > N
        rangeSignal = rangeSignal(1:N); % Truncate
    elseif length(rangeSignal) < N
        rangeSignal = [rangeSignal; zeros(N - length(rangeSignal), 1)]; % Zero-pad
    end
    % compute fft
    y = fft(rangeSignal, N); % fixed length fft
    p2 = abs(y / N);
    p1 = p2(1:floor(N / 2) + 1);
    p1(2:end-1) = 2 * p1(2:end-1);
    % accumulated power spec
    powerAccum = powerAccum + p1.^2;
    fftRealAccum = fftRealAccum + real(y);
    fftImagAccum = fftImagAccum + imag(y);
end


% Averaged Power Spectrum (0–5000 Hz)
avePower = powerAccum / numRange;
% average phase spectrum
aveFFT = complex(fftRealAccum / numRange, fftImagAccum / numRange);
avePhase = angle(aveFFT(1:floor(N / 2) + 1));

figure;
plot(f, avePower);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Averaged Velocity Signal Power Spectrum (0–5000 Hz)');
xlim([0 1000]);

% Plot averaged phase spectrum
figure;
plot(f, unwrap(avePhase));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Averaged Velocity Signal Phase Spectrum (0–5000 Hz)');
xlim([0 1000]);


%%%Reconstructed Average Velocity Signal from Averaged FFTPower + Phase

aveMagnitude = sqrt(avePower);
aveFFT_reconstructed = aveMagnitude .* exp(1i * avePhase);

%mirror positive to get negative
aveFFT_full = [aveFFT_reconstructed; conj(flipud(aveFFT_reconstructed(2:end-1)))];
% zero pad full fft to match N
aveFFT_full = [aveFFT_full; zeros(N - length(aveFFT_full), 1)];

% reconstructed signal
aveVelocity_reconstructed = ifft(aveFFT_full, N);
figure;
plot(real(aveVelocity_reconstructed));
xlabel('Time (samples)');
ylabel('Amplitude');
title('Reconstructed Average Velocity Signal (IFFT)');
xlim([0 N]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Current Velocity / Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Work
velocity = useTemp(useRanges{curr})*sf;
velocity = velocity - mean(velocity);
% % % FIR Highpass Filter Design
order = 700; % wo finger
cutoffFreq = 11 / (Fs/2);
b = fir1(order, cutoffFreq, 'high');
velocityhp = filter(b, 1, velocity);

% velocityhp = velocity;

figure;
plot(velocity);
hold on;
plot(velocityhp);
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity (mm/s)');
legend('Raw','HighPassed');
xticks(0:5000:25000);
xticklabels(0:1:5);


y = fft(velocity, N);
p2 = abs(y / N);
p1 = p2(1:floor(N / 2) + 1);
p1(2:end-1) = 2 * p1(2:end-1);

figure;
plot(f, p1.^2);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Current Velocity Power Spectrum');
xlim([0 1000]);

%%%%%%%%%%%%%Integrating to find position
figure;
dispRaw = cumsum(velocity)*T;
disphp = cumsum(velocityhp)*T;
plot(dispRaw);
hold on;
plot(disphp);
title('Displacement vs. Time');
legend('Raw','HighPassed');
xlabel('Time (s)');
xticks(0:5000:25000);
xticklabels(0:1:5);
ylabel('Displacement (mm)');

% Power Spectrum and Phase Spectrum for Displacement
L = length(dispRaw);
Y_raw = fft(dispRaw);
P2_raw = abs(Y_raw / L);
P1_raw = P2_raw(1:floor(L / 2) + 1);
P1_raw(2:end-1) = 2 * P1_raw(2:end-1);

Y_hp = fft(disphp);
P2_hp = abs(Y_hp / L);
P1_hp = P2_hp(1:floor(L / 2) + 1);
P1_hp(2:end-1) = 2 * P1_hp(2:end-1);
phase_disp_hp = angle(Y_hp);

% displacement frequencies
f_disp = Fs * (0:(floor(L / 2))) / L; % Same as f but for displacement

% Plot Power Spectrum of Displacement (HighPassed)
figure;
plot(f_disp, P1_hp.^2);
title('Current Displacement Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power');
legend('HighPassed Displacement');
ylim([0 1.5e-4]);
xlim([0 500]);

% Plot Phase Spectrum of HighPassed Displacement
figure;
plot(f_disp, unwrap(phase_disp_hp(1:floor(L / 2) + 1)));
hold on;
a = 1;
for k = -a:a
    yline(k * pi, 'g');
    yline(k / 2 * pi, 'r');
end

title('Phase Spectrum of HP Displacement (With Finger)');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
xlim([0 500]);

window = hamming(256);
noverlap = 128;
nfft = 1024;
figure;
spectrogram(useTemp(useRanges{curr}), window, noverlap, nfft, Fs);
[s, f, t, p] = spectrogram(useTemp(useRanges{curr}), window, noverlap, nfft, Fs);

figure; 
spectrogram(real(aveVelocity_reconstructed), window, noverlap, nfft, Fs);
[save, fave, tave, pave] = spectrogram(real(aveVelocity_reconstructed), window, noverlap, nfft, Fs);
