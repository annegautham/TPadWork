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
plot(f, p1.^2)