function sendToTeensy(serialPort, signal)
    % Open serial connection
    if ~isempty(serialportfind('Port', serialPort))
        fclose(serialportfind('Port', serialPort));  % Close any existing connection
    end
    serialObj = serialport(serialPort, 115200);

    % Normalize signal to 12-bit DAC range (0 to 4095)
    signal = round((signal + 1) / 2 * 4095);  % Scale from [-1,1] to [0,4095]
    signal = max(0, min(4095, signal));  % Clip values to 12-bit range
    
    % Send signal data
    for i = 1:length(signal)
        write(serialObj, typecast(uint16(signal(i)), 'uint8'), 'uint8');
        pause(1e-4);  % Small delay to ensure reliable transmission
    end
    
    % Close serial connection
    clear serialObj;
end
