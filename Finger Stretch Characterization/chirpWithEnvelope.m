function [signal, t] = chirpWithEnvelope(freqStart, freqEnd, amplitude, duration, fs)
t = 0:1/fs:duration;
chirpSignal = amplitude * chirp(t, freqStart, duration, freqEnd);
envelope = tukeywin(length(t), min(1, 0.5/duration))';
signal = chirpSignal .* envelope;
end

