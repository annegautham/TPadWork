% params
m = 1;
L = 1;
b = 1;
theta_r = 1;

% sampling
g = rand() * 10;                   % Uniform in [0, 10]
phi_g = (2 * rand() - 1) * pi;     % Uniform in [-pi, pi]
theta0 = (2 * rand() - 1) * pi;    % Initial angle
thetadot0 = (2 * rand() - 1) * pi; % Initial angular velocity

% ode function
f = @(t, x) [...
    x(2); % dx1/dt = x2
    - (b / (m * L^2)) * x(2) - (g / L) * sin(x(1) - theta_r) % dx2/dt
];

% solutoin
tspan = [0 20];
x0 = [theta0; thetadot0];
[t, x] = ode45(f, tspan, x0);

% plot over time
figure;
plot(t, x(:,1));
hold on;
yline(theta_r, '--r');
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('\theta(t)', '\theta_r');

% params
m = 1; L = 1; b = 1;
theta_r = 1;

% sampling
g = 10 * rand();
phi_g = -pi + 2*pi*rand();
theta0 = -pi + 2*pi*rand();
theta_dot0 = -pi + 2*pi*rand();
x0 = [theta0; theta_dot0];

% control gain
k = 15;
tspan = [0 10];

% ode function
f = @(t, x) [
    x(2);
    -(b/(m*L^2))*x(2) ...
    - (g/L)*sin(x(1) - phi_g) ...
    - (k/(m*L^2)) * sin(x(1) - theta_r)
];

[t, x] = ode45(f, tspan, x0);
theta_wrapped = mod(x(:,1) + pi, 2*pi) - pi;

% Plot
figure;
plot(t, theta_wrapped);
hold on;
yline(mod(theta_r + pi, 2*pi) - pi, '--r', 'θ_r');
xlabel('Time (s)');
ylabel('θ (rad)');
legend('θ(t)', 'θ_r');

% params
m = 1; L = 1; b = 1;
theta_r = 1;         % Desired angle

% sampling
g = 10 * rand();
phi_g = -pi + 2*pi*rand();
theta0 = -pi + 2*pi*rand();
omega0 = -pi + 2*pi*rand();
z0 = 0;
x0 = [theta0; omega0; z0];

% control gains
k_p = 10;
k_i = 5;
tspan = [0 15];

% ode function
f = @(t, x) [
    x(2);
    -(b / (m*L^2)) * x(2) - (g / L) * sin(x(1) - phi_g) ...
        + (1 / (m*L^2)) * (-k_p * sin(x(1) - theta_r) - k_i * x(3)); % domega/dt
    sin(x(1) - theta_r)
];

% solution
[t, x] = ode45(f, tspan, x0);
theta_wrapped = mod(x(:,1) + pi, 2*pi) - pi;

% Plot
figure;
plot(t, theta_wrapped);
hold on;
yline(mod(theta_r + pi, 2*pi) - pi, '--r', 'θ_r');
xlabel('Time (s)');
ylabel('θ (rad)');
legend('\theta(t)', '\theta_r');
