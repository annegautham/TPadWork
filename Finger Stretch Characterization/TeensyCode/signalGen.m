function funs = signalGen(serialPort)
    funs.sineWithEnvelope = @(freq, amplitude, duration, fs) ...
        sendSignal(sineWithEnvelope(freq, amplitude, duration, fs), serialPort, fs);
    
    funs.chirpWithEnvelope = @(freqStart, freqEnd, amplitude, duration, fs) ...
        sendSignal(chirpWithEnvelope(freqStart, freqEnd, amplitude, duration, fs), serialPort, fs);
end

%% **Sine wave with Tukey envelope**
function [signal, t] = sineWithEnvelope(freq, amplitude, duration, fs)
    t = 0:1/fs:duration;
    sineSignal = amplitude * sin(2 * pi * freq * t); 
    envelope = tukeywin(length(t), min(1, 0.5 / duration))';  % Tukey window
    signal = sineSignal .* envelope;  
    signal = sineSignal;
    [signal, t] = scaleSignal(signal, amplitude, t);
end

%% **Chirp signal with Tukey envelope**
function [signal, t] = chirpWithEnvelope(freqStart, freqEnd, amplitude, duration, fs)
    t = 0:1/fs:duration;
    chirpSignal = amplitude * chirp(t, freqStart, duration, freqEnd);
    envelope = tukeywin(length(t), min(1, 0.5 / duration))';
    signal = chirpSignal .* envelope;
    signal = chirpSignal;
    [signal, t] = scaleSignal(signal, amplitude, t);
end

%% **Scale signal to 12-bit DAC (0-4095)**
function [scaledSignal, t] = scaleSignal(signal, amplitude, t)
    DAC_max = 4095;  % 12-bit DAC resolution
    V_ref = 3.3;  % Teensy DAC reference voltage

    % Normalize by actual max amplitude (not forcing full scale)
    signal_max = max(abs(signal));  % Find peak
    if signal_max == 0
        signal_max = 1;  % Avoid divide-by-zero
    end

    % Scale based on requested amplitude
    scaleFactor = (amplitude / V_ref);  % Convert to DAC proportion
    scaledSignal = round((signal / signal_max) * (DAC_max / 2) * scaleFactor + (DAC_max / 2));

    % Ensure values stay in range
    scaledSignal = max(0, min(scaledSignal, DAC_max));
end



%% **Send signal to Teensy**
function sendSignal(signal, serialPort, fs)
    t = (0:length(signal)-1) / fs;  % Time axis
    figure;
    plot(t, signal, 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Amplitude (DAC units)');

    serialObj = serialport(serialPort, 115200);
    pause(2);  
    csvData = sprintf('%d,', signal);
    writeline(serialObj, csvData(1:end-1));
    writeline(serialObj, "END");
    pause(1);
    clear serialObj;
end
