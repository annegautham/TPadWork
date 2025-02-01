function funs = signalGen()
    % Returns a struct of function handles for generating signals
    funs.sineWithEnvelope = @sineWithEnvelope;
    funs.chirpWithEnvelope = @chirpWithEnvelope;
end

%% Sine wave with Tukey window envelope
function [signal, t] = sineWithEnvelope(freq, amplitude, duration, fs)
    t = 0:1/fs:duration; % Time vector
    sineSignal = amplitude * sin(2 * pi * freq * t); % Generate sine wave
    rampTime = 0.25;
    tukeyAlpha = min(1, (rampTime * 2 / duration)); % Ensure alpha ≤ 1
    envelope = tukeywin(length(t), tukeyAlpha)'; % Tukey window as row vector
    signal = sineSignal .* envelope; % Apply envelope
    titleStr = sprintf('Freq: %.2f Hz, Amplitude: %.1f V', freq, amplitude);
    plotSignal(t, signal, titleStr);
end

%% Chirp signal with Tukey window envelope
function [signal, t] = chirpWithEnvelope(freqStart, freqEnd, amplitude, duration, fs)
    t = 0:1/fs:duration; % Time vector
    chirpSignal = amplitude * chirp(t, freqStart, duration, freqEnd); % Generate chirp
    rampTime = 0.25;
    tukeyAlpha = min(1, (rampTime * 2 / duration)); % Ensure alpha ≤ 1
    envelope = tukeywin(length(t), tukeyAlpha)'; % Tukey window as row vector
    signal = chirpSignal .* envelope; % Apply envelope
    titleStr = sprintf('Freq: %.2f to %.2f Hz, Amplitude: %.1f V', freqStart, freqEnd, amplitude);
    plotSignal(t, signal, titleStr);
end

%% Function to send signal to Teensy
function sendToTeensy(serialPort, signal, fs)
    signal = (signal + 1) / 2; % Scale to 0-1
    signal = uint16(signal * 4095); % Convert to 12-bit DAC values
    for i = 1:length(signal)
        highByte = bitshift(signal(i), -8);
        lowByte = bitand(signal(i), 255);
        write(serialPort, [highByte, lowByte], "uint8");
        pause(1/fs);
    end
end

function plotSignal(t, signal, titleStr)
    figure;
    plot(t, signal, 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(titleStr);
end
