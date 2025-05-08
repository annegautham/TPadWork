clear; close all; clc

T = 50;
dt = 0.01;
t = 0:dt:T;

V_all = zeros(3, length(t));
labels = {'Case 1: Const gain, \omega_1 = \omega_2', ...
          'Case 2: Const gain, \omega_1 \neq \omega_2', ...
          'Case 3: PI control'};
%params
theta1_0 = 2*pi*rand;
theta2_0 = 2*pi*rand;

% Random gains and frequencies
a1 = 1 + rand;
a2 = 1 + rand;
K1 = 1 + rand;
K2 = 1 + rand;
L1 = 1 + rand;
L2 = 1 + rand;

% phi = 0
omega1 = 2*pi*10;   % 10 Hz
omega2 = omega1;

theta1 = theta1_0;
theta2 = theta2_0;

for k = 1:length(t)
    V_all(1, k) = 1 - cos(theta1 - theta2);

    dtheta1 = omega1 + a1*K1*sin(theta2 - theta1);
    dtheta2 = omega2 + a2*K2*sin(theta1 - theta2);

    theta1 = theta1 + dt*dtheta1;
    theta2 = theta2 + dt*dtheta2;
end

% Constant gain, phase diff not equal to 0
omega1 = 2*pi*10;
omega2 = 2*pi*(10 + 1.5*rand);  % Slight offset

theta1 = theta1_0;
theta2 = theta2_0;

for k = 1:length(t)
    V_all(2, k) = 1 - cos(theta1 - theta2);

    dtheta1 = omega1 + a1*K1*sin(theta2 - theta1);
    dtheta2 = omega2 + a2*K2*sin(theta1 - theta2);

    theta1 = theta1 + dt*dtheta1;
    theta2 = theta2 + dt*dtheta2;
end

%PI Control
omega1 = 2*pi*10;
omega2 = 2*pi*(10 + 1.5*rand);

theta1 = theta1_0;
theta2 = theta2_0;
z1 = randn;
z2 = randn;

for k = 1:length(t)
    V_all(3, k) = 1 - cos(theta1 - theta2);
    
    sin_diff = sin(theta2 - theta1);
    
    dtheta1 = omega1 + a1*K1*sin_diff + a1*L1*z1;
    dtheta2 = omega2 + a2*K2*sin(-sin_diff) + a2*L2*z2;
    
    dz1 = sin_diff;
    dz2 = -sin_diff;

    theta1 = theta1 + dt*dtheta1;
    theta2 = theta2 + dt*dtheta2;
    z1 = z1 + dt*dz1;
    z2 = z2 + dt*dz2;
end

%% plot
figure;
semilogy(t, V_all(1,:), 'LineWidth', 1.5); hold on;
semilogy(t, V_all(2,:), 'LineWidth', 1.5);
semilogy(t, V_all(3,:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('V(\theta_1, \theta_2)');
title('Synchronization Metric V(t) for Three PLL Cases');
legend(labels, 'Location', 'northeast');

