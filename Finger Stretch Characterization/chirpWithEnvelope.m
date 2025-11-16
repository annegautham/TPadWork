% function [signal, t] = chirpWithEnvelope(freqStart, freqEnd, amplitude, duration, fs)
%     t = 0:1/fs:duration;
%     chirpSignal = amplitude * chirp(t, freqStart, duration, freqEnd);
%     envelope = tukeywin(length(t), min(1, 0.5 / duration))';
%     signal = chirpSignal .* envelope;
%     % [signal, t] = scaleSignal(signal, amplitude, t);
% end
% 


function [signal, t, current, instFreq] = chirpWithEnvelope( ...
        freqStart, freqEnd, amplitude, duration, fs, ...
        freq_interp, ohm_interp, doPlot)

    if nargin < 8, doPlot = true; end

    % --- time & chirp ---
    t = (0:1/fs:duration).';                          % column
    chirpSignal = amplitude * chirp(t, freqStart, duration, freqEnd, 'linear');
    envelope = tukeywin(length(t), min(1, 0.5 / duration));
    signal = chirpSignal .* envelope;                 % voltage vs time

    % --- instantaneous frequency of linear chirp ---
    instFreq = freqStart + (freqEnd - freqStart) * (t / duration);

    % --- interpolate impedance at f(t) and compute current ---
    % freq_interp, ohm_interp come from your .mat (1x1000 doubles)
    if nargin >= 7 && ~isempty(freq_interp) && ~isempty(ohm_interp)
        fZ = freq_interp(:); Zmag = ohm_interp(:);
        % Keep interpolation safe
        % Option A: extrapolate (use with caution)
        % Z_t = interp1(fZ, Zmag, instFreq, 'pchip', 'extrap');
        % Option B: clamp to valid band (safer around ends)
        instFreq_clamped = min(max(instFreq, min(fZ)), max(fZ));
        Z_t = interp1(fZ, Zmag, instFreq_clamped, 'pchip');

        % Avoid divide-by-zero at sharp resonances
        Z_t(Z_t < 1e-6) = 1e-6;

        current = signal ./ Z_t;                      % i(t) in Amps
    else
        current = [];
    end

    % --- quick plot ---
    if doPlot
        figure('Name','Chirp Voltage & Current');
        tiledlayout(3,1);

        nexttile;
        plot(t, signal, 'LineWidth', 1.2); grid on;
        xlabel('Time (s)'); ylabel('V(t) [V]');
        title('Chirp Voltage');

        nexttile;
        plot(t, instFreq, 'LineWidth', 1.2); grid on;
        xlabel('Time (s)'); ylabel('f(t) [Hz]');
        title('Instantaneous Frequency');

        nexttile;
        if ~isempty(current)
            plot(t, current, 'LineWidth', 1.2); grid on;
            xlabel('Time (s)'); ylabel('i(t) [A]');
            title('Instantaneous Current (V/Z(f(t)))');
        else
            text(0.1,0.5,'(no impedance provided)', 'FontSize', 12);
            axis off
        end
    end
end
