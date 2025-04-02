clc; clear; close all;

% Given open-loop transfer function
K = 100000;  % Gain determined for steady-state error requirement
G = tf(K, [1 205 1000 0]); % G(s) = K / (s(s+5)(s+200))

% Bode plot of the uncompensated system
figure;
margin(G);
grid on;
title('Bode Plot of Uncompensated System');

% Lead compensator design
alpha = 1/3; % Chosen to provide 30° phase lead
omega_m = 40; % Desired frequency for maximum phase lead
T = 1 / (omega_m^2 * alpha); % Compute T

% Compensator transfer function
C = tf([T 1], [alpha*T 1]);

% Compensated open-loop system
G_comp = C * G;

% Bode plot of the compensated system
figure;
margin(G_comp);
grid on;
title('Bode Plot of Compensated System');

% Closed-loop system
T_cl = feedback(G_comp, 1);

% Step response to check transient behavior
figure;
step(T_cl);
grid on;
title('Step Response of Compensated System');

% Check closed-loop poles
poles = pole(T_cl);
disp('Closed-loop poles:');
disp(poles);

% Damping ratio and natural frequency
[zeta, wn] = damp(T_cl);
disp('Damping Ratios:');
disp(zeta);
disp('Natural Frequencies:');
disp(wn);
