%% Parameters (using fs = 5000 Hz)
fs   = 5000;        % Sampling frequency in Hz
T    = 15;          % Duration of the sweep in seconds
f1   = 22;          % Start frequency (Hz)
f2   = 2500;        % End frequency (Hz), limited by Nyquist (fs/2)

% Create a time vector for the full sweep duration
t = linspace(0, T, fs * T);

%% Generate the Test Sweep Signal (x)
% Here we use a linear chirp as the test signal.
x = chirp(t, f1, T, f2, 'linear');

%% Simulate a System Under Test
% For demonstration purposes, simulate an impulse response (e.g., exponential decay)
h = exp(-0.1 * t);  % Example impulse response

% The measured (recorded) output signal is obtained by convolving the test signal x with h.
useTemp = conv(x, h, 'same');  
% In practice, useTemp is your measured output, not a simulated one.

%% Define Multiple Ranges (in seconds)
% Each row defines a [start_time, end_time] (in seconds) that contains one measurement sweep.
ranges = [0, 5; 5, 10; 10, 15];  

%% Process Each Range: Deconvolve Using the Inverse Filter
% For each range, we assume that the relevant test signal segment (from useTemp)
% should be deconvolved by its time-reversed version to retrieve the impulse response.
for i = 1:size(ranges, 1)
    % Define start and end times (in seconds) for the current range
    start_sec = ranges(i, 1);
    end_sec   = ranges(i, 2);
    
    % Convert these times to sample indices
    start_idx = max(1, round(start_sec * fs));
    end_idx   = min(length(useTemp), round(end_sec * fs));
    
    % Extract the segment from useTemp corresponding to the current sweep
    segment = useTemp(start_idx:end_idx);
    t_segment = t(start_idx:end_idx);
    
    % Create the inverse filter as the time-reversal of the segment.
    % (In the ideal case, the original test signal is known, and f(t) = time-reversed test signal.)
    inverse_filter = flip(segment);
    
    % Deconvolve: convolve the segment with its inverse filter to recover h(t)
    impulse_response = conv(segment, inverse_filter, 'same');
    
    %% Plot the Results for the Current Range
    figure;
    
    % Plot the extracted segment from useTemp
    subplot(3,1,1);
    plot(t_segment, segment);
    title(['Recorded Segment ' num2str(i) ' from useTemp']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot the inverse filter (time-reversed segment)
    subplot(3,1,2);
    plot(t_segment, inverse_filter);
    title(['Inverse Filter (Time-Reversed) for Segment ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot the recovered impulse response
    subplot(3,1,3);
    plot(t_segment, impulse_response);
    title(['Recovered Impulse Response for Segment ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    %% (Optional) Post-Processing: Apply a Time Window
    % For example, if you want to isolate only the linear response portion,
    % you can define a time window (adjust these values as needed).
    window_start = round(0.5 * fs); % Example: start window at 0.5 s into the segment
    window_end   = round(1.5 * fs); % Example: end window at 1.5 s into the segment
    if window_end <= length(impulse_response)
        h_filtered = impulse_response(window_start:window_end);
        t_filtered = t_segment(window_start:window_end);
        
        figure;
        plot(t_filtered, h_filtered, 'm');
        title(['Windowed (Filtered) Impulse Response for Segment ' num2str(i)]);
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
end
