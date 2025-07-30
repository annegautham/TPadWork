

%System matrix
A = [-75/16, 1.5;
      0.75, -3/8];

% 2nd to first order system
% x = [theta1; theta2; theta1_dot; theta2_dot]
function dXdt = state_space(~, X)
    A = [-75/16, 1.5; 0.75, -3/8];
    theta = X(1:2);
    theta_dot = X(3:4);
    theta_ddot = A * theta;
    dXdt = [theta_dot; theta_ddot];
end

% time span
tspan = [0 20];
t_eval = linspace(0, 20, 1000);

% init conditions
ICs = {
    [0.1; 0.1; 0; 0], 'Case I: \theta_1(0)=0.1, \theta_2(0)=0.1'
    [0.1; -0.1; 0; 0], 'Case II: \theta_1(0)=0.1, \theta_2(0)=-0.1'
    [0.1; 0; 0; 0], 'Case III: \theta_1(0)=0.1, \theta_2(0)=0'
};

% plotting
figure;
for i = 1:3
    ic = ICs{i, 1};
    label = ICs{i, 2};
    [~, X] = ode45(@state_space, t_eval, ic);
    subplot(3,1,i);
    plot(t_eval, X(:,1), 'b', 'DisplayName', '\theta_1(t)');
    hold on;
    plot(t_eval, X(:,2), 'r', 'DisplayName', '\theta_2(t)');
    title(label);
    ylabel('Angle (rad)');
    legend;
end
xlabel('Time (s)');
