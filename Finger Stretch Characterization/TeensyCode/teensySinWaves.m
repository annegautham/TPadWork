% Parameters
Fs = 50000; % Sampling frequency (Hz)
f = 100;    % Sine wave frequency (Hz)
duration = 10; % Duration in seconds
amplitude = 2047; % Max value for 12-bit DAC (half-scale for sine wave)

% Generate time vector
t = 0:1/Fs:duration - 1/Fs;

% Generate sine wave (scaled for 12-bit DAC)
sineWave = amplitude * (1 + sin(2 * pi * f * t)); % Offset to avoid negative values

% Open serial connection to Teensy
serialPort = 'COM4'; % Replace with your Teensy COM port
baudRate = 115200;
s = serialport(serialPort, baudRate);

% Send data to Teensy
for i = 1:length(sineWave)
    % Scale and send data as 16-bit integers
    value = uint16(sineWave(i));
    write(s, value, 'uint16');
    pause(1/Fs); % Ensure correct timing
end

% Close serial connection after the duration
clear s;
disp('Sine wave transmission complete.');
