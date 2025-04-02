clear
close all

serialPort = "COM4";  % Change to your Teensy port
funs = signalGen(serialPort);

dq = daq("ni");

time = 10;  
num_trials = 1;  % Define number of trials
tvec = cell(1, num_trials);  % Store each time vardector separately
tempVec = cell(1, num_trials);  % Store each signal separately

%% DAQ Configuration
dq.Rate = 5000;
channel_list = [13];
ch = addinput(dq, "Dev1", channel_list, "Voltage");

%% Countdown Timer
disp('Starting in:');
for i = 5:-1:1  
    fprintf('%d...\n', i);
    pause(1);
end

%% Run Chirp Signals and Record Data
for i = 1:num_trials
    fprintf('Running chirp %d/%d...\n', i, num_trials);
    
    funs.lin_chirpWithEnvelope(30, 350, 0.8, 5, 10000);  % Start chirp signal
    disp('Recording...');
    
    Data_dq = read(dq, seconds(time));  % Read DAQ data
    
    disp('Stopped...');
    
    % Store data properly
    tvec{i} = Data_dq.Time * 1e3;  % Convert to milliseconds
    tempVec{i} = Data_dq.(Data_dq.Properties.VariableNames{1});
end

clearSerialPort(serialPort);  % Close serial port properly

%% Plot Results
figure;
hold on;
for i = 1:num_trials
    plot(tvec{i}, tempVec{i});
    hold on;
end
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Chirp Response');
legend(arrayfun(@(x) sprintf('Trial %d', x), 1:num_trials, 'UniformOutput', false));
hold off;

%% Save Data
save('30-350_DAC_Output.mat', 'tempVec', 'tvec', '-v7.3');

%% Function to Clear Serial Port
function clearSerialPort(port)
    serialObj = serialportfind('Port', port); % Find existing serial object
    if ~isempty(serialObj)
        fclose(serialObj);
        delete(serialObj);
        clear serialObj;
    end
    disp("Serial port cleared.");
end
