close all;
Fs = 5000;
T = 1 / Fs;
L = length(disphp);
t = (0:L-1) * T;

f0 = 30;
f1 = 300;
deltaT = 3.1958;
C = (f1 - f0)/deltaT;
i_t = 18 * sin(2 * pi * t .* (f0 + C * t / 2));

I = fft(i_t);
p2 = abs(I/L);
p1 = p2(1:L/2+1);
p1(2:end-1) = 2*p1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f, p1.^2)
xlabel('Frequency (Hz)');
ylabel('Power')
title('30-300, deltaT = 3.1958 Chirp')
xlim([0 300]);

N = length(i_t);
Y_hp = fft(disphp);
H = Y_hp ./ I;

P2_H = abs(H/L);
P1_H = P2_H(1:L/2+1);
P1_H(2:end-1) = 2*P1_H(2:end-1);

figure;
plot(f, P1_H.^2)
xlabel('Frequency (Hz)');
ylabel('Power')
title('H')
xlim([0 500]);
% % Remove DC bias
% 
% % FFT of input and output
% N = length(I_t);
% U = fft(I_t);
% Y = fft(temp);
% 
% % Frequency vector (only positive frequencies)
% f = (0:N-1) * (Fs/N);
% f = f(1:N/2);
% U = U(1:N/2);
% Y = Y(1:N/2);
% 
% % Estimate Transfer Function in Frequency Domain
% H_raw = Y ./ U;
% 
% % Use only magnitude response for fitting
% H_mag = abs(H_raw);
% 
% % Angular frequency vector
% omega = 2 * pi * f(:); % Ensuring it's a column vector
% 
% % Manual Rational Polynomial Fit: 2nd order system approximation
% % Assume H(s) = (b0 + b1*s + b2*s^2) / (1 + a1*s + a2*s^2)
% X_num = [ones(length(omega), 1), 1i*omega, (1i*omega).^2]; % Numerator terms
% X_den = [-1i*omega, -(1i*omega).^2]; % Denominator terms
% 
% % Calculate numerator and denominator coefficients
% b_coeffs = real(X_num \ H_mag); % Fit numerator [b0, b1, b2]
% a_coeffs = real(X_den \ (H_mag - X_num*b_coeffs)); % Fit denominator [a1, a2]
% 
% % Construct numerator and denominator polynomials
% numerator = [b_coeffs(3), b_coeffs(2), b_coeffs(1)]; % Flip order
% denominator = [1, a_coeffs(2), a_coeffs(1)]; % Standard form
% 
% % Find poles and zeros
% poles = roots(denominator);
% zeros = roots(numerator);
% 
% % Display poles and zeros
% disp('Poles of the system:');
% disp(poles);
% disp('Zeros of the system:');
% disp(zeros);
% 
% % Calculate Resonant Frequency (highest imaginary part of poles)
% [~, idx] = max(abs(imag(poles)));
% resonant_freq = abs(imag(poles(idx))) / (2 * pi);
% disp(['Resonant Frequency (Hz): ', num2str(resonant_freq)]);
% 
% % Frequency response of the estimated transfer function
% H_est = polyval(numerator, 1i * omega) ./ polyval(denominator, 1i * omega);
% 
% % Plot Raw H(f) vs Estimated Transfer Function
% figure;
% subplot(2, 1, 1);
% plot(f, abs(H_raw), 'b', 'LineWidth', 1.5);
% hold on;
% plot(f, abs(H_est), 'r--', 'LineWidth', 1.5);
% title('Magnitude of Transfer Function');
% xlabel('Frequency (Hz)');
% ylabel('|H(f)|');
% legend('Raw H(f)', 'Estimated H(f)');
% grid on;
% 
% subplot(2, 1, 2);
% plot(f, angle(H_raw), 'b', 'LineWidth', 1.5);
% hold on;
% plot(f, angle(H_est), 'r--', 'LineWidth', 1.5);
% title('Phase of Transfer Function');
% xlabel('Frequency (Hz)');
% ylabel('Phase (Radians)');
% legend('Raw H(f)', 'Estimated H(f)');
% grid on;