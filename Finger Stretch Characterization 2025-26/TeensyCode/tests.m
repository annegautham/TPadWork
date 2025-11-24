clear; close all; clc;

% ===== USER PARAMETERS =====
serialPort = "COM3";          % Teensy port
f1 = 30;                      % Start frequency (Hz)
f2 = 350;                     % End frequency (Hz)
duration = 5;                 % Duration (s)
env = 0.8;                    % Envelope scaling
fs_teensy = 10000;            % Teensy sample rate (command)
daqTime = 10;                 % DAQ record time (s)

% ===== OPEN SERIAL =====
teensy = serialport(serialPort, 115200);
pause(0.5);

% ===== SEND CHIRP COMMAND TO TEENSY =====
cmd = sprintf("CHIRP %.2f %.2f %.2f %.2f %.2f\n", f1, f2, duration, env, fs_teensy);
write(teensy, cmd, "char");
disp("Sent:");
disp(cmd);

% ===== CONFIGURE DAQ =====
dq = daq("ni");
dq.Rate = 10000;
addinput(dq, "Dev1", [13 1], "Voltage");   % Channels 13 (output), 1 (input)
ok aw
% ===== ACQUIRE DATA =====
disp("Recording...");
data = read(dq, seconds(daqTime));
disp("Done.");

t = data.Time;
out = data.("Dev1_ai13");     % Output channel
inp = data.("Dev1_ai1");      % Input channel

% ===== CLOSE SERIAL =====
delete(teensy);
clear teensy;
disp("Serial closed.");

% ===== PLOTS =====
figure;
plot(t, out);
xlabel("Time (s)");
ylabel("Output (V)");
title("Output Channel 13");

figure;
plot(t, inp);
xlabel("Time (s)");
ylabel("Input (V)");
title("Input Channel 1");

% ===== SAVE RAW DATA =====
save("chirp_data.mat", "t", "out", "inp", "-v7.3");

disp("Saved chirp_data.mat");
