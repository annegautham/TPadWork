signal = temp;
Fs = 5000;
T = 1/Fs;
L = length(signal);
t= (0:L-1)*T;

y = fft(signal);
p2 = abs(y/L);
p1 = p2(1:L/2+1);
p1(2:end-1) = 2*p1(2:end-1);

f = Fs*(0:(L/2))/L;
figure;
plot(f, p1.^2)
title("150 Hz @ 16 V")
xlabel('Frequency (Hz)');
ylabel('Power')

figure;
velocity = temp(1:500)*125;

order = 10;
cutoffFreq = 10 / (Fs/2);
b = fir1(order, cutoffFreq, 'high');
velocityhp = filter(b, 1, velocity);

plot(velocity);
hold on;
plot(velocityhp);
title('Velocity vs. Time');
xlabel('Time (Hz)');
ylabel('Velocity (mm/s)');
legend('Raw','HighPassed');


figure;
dispRaw = cumsum(velocity)*T;
disphp = cumsum(velocityhp)*T;
plot(dispRaw);
hold on;
plot(disphp);
title('Displacement vs. Time');
legend('Raw','HighPassed');
xlabel('Time (Hz)');
ylabel('Displacement (mm)');


