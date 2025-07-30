clear; clc;

% params
M = 10;     
m = 2;
k = 10; 
c = 0.01;  
F0 = 0.2; 

% Mass damping and stiff matrices
MM = [M, 0;
      0, m];
CC = [2*c, -c;
      -c,  c];
KK = [2*k, -k;
      -k,  k];

% Natural frequencies (undamped)
omega1 = sqrt(k/M);  % ≈ 1.0 rad/s
omega2 = sqrt(k/m);  % ≈ 2.236 rad/s

% freq range
omega = linspace(0.01, 4*omega2, 1000); % avoid omega = 0

% amp vectors
X1_mag = zeros(size(omega));
X2_mag = zeros(size(omega));

% freq response
for i = 1:length(omega)
    w = omega(i);
    A = -w^2 * MM + 1i * w * CC + KK;
    F = [0; F0];
    X = A \ F;
    X1_mag(i) = abs(X(1));
    X2_mag(i) = abs(X(2));
end

% plots
figure;
plot(omega, X1_mag, 'r', 'LineWidth', 1.5); hold on;
plot(omega, X2_mag, 'b', 'LineWidth', 1.5);
xlabel('\omega (rad/s)');
ylabel('Steady-State Displacement Amplitude');
title('Frequency Response of Platform (M) and Block (m)');
legend('|X_1(\omega)| (Platform)', '|X_2(\omega)| (Block)');