close all

fs = 10000;
Tchirp = 5;
freqStart = 30;
freqEnd = 350;
amplitude = 0.7;

[chirpTest, t_chirp] = chirpWithEnvelope(freqStart, freqEnd, amplitude, Tchirp, fs);
invFilter = flip(chirpTest);
num_trials = numel(tempVec);

for i = 1:num_trials
    y = tempVec{i}(5100:30100);
    t_y = (0:length(y)-1)/fs;
    figure;
    plot(t_y, y, 'b');
    hold on;
    plot(t_chirp, invFilter, 'r');
    title(sprintf('Output y and Inverse Filter for Trial %d', i));
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    h_rec = conv(y, invFilter, 'full');
    h_rec = h_rec(1:length(y)); 
    t_h = (0:length(h_rec)-1)/fs;
    figure;
    plot(t_h, h_rec, 'LineWidth',1.2);
    title(sprintf('Impulse Response for Trial %d', i));
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    H = fft(h_rec);
    nfft = length(H);
    H = H(1:nfft/2+1);
    freq = linspace(0, fs/2, nfft/2+1);
    
    figure;
    plot(freq, abs(H), 'LineWidth',1.2);
    title(sprintf('H(f) for Trial %d', i));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    
    epsilon = 1e-6;
    C = conj(H)./(abs(H).^2 + epsilon);
    c_t = ifft(C, 'symmetric');
    c_t = c_t - mean(c_t);
    
    c_t = c_t(:); % Ensure column vector
    c_t = [c_t(ceil(length(c_t)/2):end); zeros(ceil(length(c_t)/2)-1, 1)];

    t_c = (0:length(c_t)-1)/fs;
    figure;
    plot(t_c, real(c_t), 'LineWidth',1.2);
    title(sprintf('c(t) for Trial %d', i));
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    h_filtered = conv(h_rec, c_t, 'full');
    h_filtered = h_filtered(1:length(y));  
    t_hf = (0:length(h_filtered)-1)/fs;
    figure;
    plot(t_hf, h_filtered, 'LineWidth',1.2);
    title(sprintf('Filtered Impulse for Trial %d', i));
    xlabel('Time (s)');
    ylabel('Amplitude');
end

function [signal, t] = chirpWithEnvelope(freqStart, freqEnd, amplitude, duration, fs)
t = 0:1/fs:duration;
chirpSignal = amplitude * chirp(t, freqStart, duration, freqEnd);
envelope = tukeywin(length(t), min(1, 0.5/duration))';
signal = chirpSignal .* envelope;
end

function [scaledSignal, t] = scaleSignal(signal, amplitude, t)
DAC_max = 4095;
V_ref = 3.3;
signal_max = max(abs(signal));
if signal_max == 0, signal_max = 1; end
scaleFactor = amplitude / V_ref;
scaledSignal = round((signal / signal_max) * (DAC_max/2) * scaleFactor + (DAC_max/2));
scaledSignal = max(0, min(scaledSignal, DAC_max));
end
