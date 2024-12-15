close all;

% ranges = {5542:21494, 21494:37658, 37658:53513, 53513:69714}; %2.8 NF Ranges
ranges = {10487:26574, 26574:42553, 42553:58608, 58608:74854}; %2.8 F Ranges
% ranges = {13003:29132, 29132:45128, 45128:61195} %3.5 NF Ranges
% ranges = {3160:19232, 19232:35196, 35196:51368, 51368:67149}; %3.5 F Ranges
sf = 125;
numRange = length(ranges);
curr = 1;
Fs = 5000; % sampling freq
T = 1 / Fs;
N = 2^14; %padding
powerAccum = zeros(floor(N / 2) + 1, 1);
fftRealAccum = zeros(N, 1);
fftImagAccum = zeros(N, 1);
f = Fs * (0:(floor(N / 2))) / N;


for i = 1:numRange
    rangeSignal = temp(ranges{i})*sf;
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
xlim([0 5000]);

% Plot averaged phase spectrum
figure;
plot(f, avePhase);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Averaged Velocity Signal Phase Spectrum (0–5000 Hz)');
xlim([0 5000]);


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
velocity = temp(ranges{curr})*sf;
velocity = velocity - mean(velocity);
% % % FIR Highpass Filter Design
% order = 700; % wo finger
% cutoffFreq = 11 / (Fs/2);
% b = fir1(order, cutoffFreq, 'high');
% velocityhp = filter(b, 1, velocity);

velocityhp = velocity;

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
xlim([0 5000]);

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
plot(f_disp, phase_disp_hp(1:floor(L / 2) + 1));
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
