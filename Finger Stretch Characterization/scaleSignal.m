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