
function [signal, t] = sineWithEnvelope(freq, amplitude, duration, fs)
    t = 0:1/fs:duration;
    sineSignal = amplitude * sin(2 * pi * freq * t); 
    envelope = tukeywin(length(t), min(1, 0.5 / duration))';  % Tukey window
    % signal = sineSignal .* envelope;
    signal = sineSignal;
end