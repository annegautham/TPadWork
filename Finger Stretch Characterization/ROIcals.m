% Given parameters
m1 = 1; m2 = 1;
b = 1;
a = b/4;
k1 = 2; k2 = 2;
k = 1;
c = 0.1;



% MOI
I1 = m1 * a^2 * (1/12);
I2 = m2 * a^2 * (1/12);

% Mass matrix M
M = [I1, 0; 0, I2];

% Stiffness matrix K
K = [k1 * a^2 + k * b^2, -k * b^2; -k * b^2, k2 * a^2 + k * b^2];

% natural freq
[eigvecs, eigvals] = eig(inv(M) * K);
natural_frequencies = sqrt(diag(eigvals)); % rad/s
C = c * b^2 * [1 -1; -1 1];


ode_func = @(t, x) [x(3); x(4); ...
                    -inv(M) * (C * [x(3); x(4)] + K * [x(1); x(2)])];


t_span = [0, 10];


initial_conditions_case1 = [1; 1; 0; 0]; % Case 1: Initial displacement 1 rad, velocity 0
initial_conditions_case2 = [-1; 1; 0; 0]; % Case 2: Initial displacement -1 rad and 1 rad, velocity 0
initial_conditions_case3 = [1; 0; 0; 0]; % Case 3: Initial displacement 1 rad, velocity 0

% Solve the ODE for each case
[t1, x1] = ode45(ode_func, t_span, initial_conditions_case1);
[t2, x2] = ode45(ode_func, t_span, initial_conditions_case2);
[t3, x3] = ode45(ode_func, t_span, initial_conditions_case3);

% plotting
figure;
subplot(3,1,1);
plot(t1, x1(:,1), 'ro');
hold on;
plot(t1, x1(:,2), 'b');
title('Case 1: Angular Displacements');
xlabel('Time (s)');
ylabel('Angular Displacement (rad)');
legend('Theta_1', 'Theta_2');

subplot(3,1,2);
plot(t2, x2(:,1), 'ro');
hold on;
plot(t2, x2(:,2), 'b');
title('Case 2: Angular Displacements');
xlabel('Time (s)');
ylabel('Ang Disp (rad)');
legend('Theta_1', 'Theta_2');

subplot(3,1,3);
plot(t3, x3(:,1), 'ro');
hold on;
plot(t3, x3(:,2), 'b', 'LineWidth',3);
title('Case 3: Angular Displacements');
xlabel('Time (s)');
ylabel('Ang Disp (rad)');
legend('Theta_1', 'Theta_2');


% Forced and damped problem (Cases 4 and 5)

F1 = 1; % Force on rod 1
F2 = 1; % Force on rod 2

omega1 = natural_frequencies(1); % First mode frequency
omega2 = natural_frequencies(2); % Second mode frequency

% Modify the ODE function to include external forces
ode_func_forced = @(t, x, omega) [x(3); x(4); ...
                                   -inv(M) * (C * [x(3); x(4)] + K * [x(1); x(2)]) + ...
                                   [F1 * cos(omega * t); F2 * cos(omega * t)]];
                               
% Solve
[t4, x4] = ode45(@(t, x) ode_func_forced(t, x, omega1), t_span, [0; 0; 0; 0]);
[t5, x5] = ode45(@(t, x) ode_func_forced(t, x, omega2), t_span, [0; 0; 0; 0]);

% Plotting
figure;
subplot(2,1,1);
plot(t4, x4(:,1), 'ro');
hold on;
plot(t4, x4(:,2), 'b', 'LineWidth', 3);
title('Case 4: Forced and Damped (omega1)');
xlabel('Time (s)');
ylabel('Angular Displacement (rad)');
legend('Theta_1', 'Theta_2');

subplot(2,1,2);
plot(t5, x5(:,1), 'ro');
hold on;
plot(t5, x5(:,2), 'b', 'LineWidth', 3);
title('Case 5: Forced and Damped (omega2)');
xlabel('Time (s)');
ylabel('Angular Displacement (rad)');
legend('Theta_1', 'Theta_2');
