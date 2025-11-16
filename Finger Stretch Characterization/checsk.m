figure;
subplot(2,1,1);
plot(f_fit, real(Z)); grid on;
ylabel('Real(Z)');

subplot(2,1,2);
plot(f_fit, imag(Z)); grid on;
ylabel('Imag(Z)');
xlabel('Frequency (Hz)');


figure;
plot(f_fit, unwrap(angle(H_fit)));
xlabel('Hz'); ylabel('phase (rad)'); grid on;
title('FRF Phase');
