
% ===== CONFIGURE DAQ =====
dq = daq("ni");
dq.Rate = 10000;
addinput(dq, "Dev2", [9], "Voltage");   % Channels 13 (output), 1 (input)
daqTime = 10;

% ===== ACQUIRE DATA =====
disp("Recording...");
data = read(dq, seconds(daqTime));
disp("Done.");

t = data.Time;
out = data.("Dev2_ai9");     % Output channel

