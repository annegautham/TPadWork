function sendSineWithEnvelope(serialPort, freq, amplitude, duration, fs)
    % Generates and sends a sine wave with a fade-in and fade-out envelope
    % Parameters:
    % - serialPort: Serial port object
    % - freq: Frequency of the sine wave (Hz)
    % - amplitude: Maximum amplitude of the signal (0 to 1)
    % - duration: Duration of the sine wave (seconds)
    % - fs: Sampling frequency (Hz)

    t = 0:1/fs:duration; % Time vector
    sineSignal = amplitude * sin(2 * pi * freq * t); % Generate sine wave

    % Create envelope (triangular shape: fade-in and fade-out)
    envelope = 1 - abs(2 * t / duration - 1); % Linear fade-in and fade-out
    sineSignal = sineSignal .* envelope; % Apply envelope

    sendToTeensy(serialPort, sineSignal, fs);
end
