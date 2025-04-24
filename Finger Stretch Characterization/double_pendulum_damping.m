function double_pendulum_damping
    % params of system
    m1 = 1; m2 = 1;
    L1 = 0.5; L2 = 0.5;
    g = 9.81;
    Tmax = 0.02; % Max torque (for actuatotr saturation)
    k = 10000;     % control gain

    % sim settin    gs
    tspan = [0 1000];
    x0 = [pi/2; pi/2; 0; 0]; % init conds [theta1; theta2; dtheta1; dtheta2]

    % u1 system
    [t1, x1, E1] = simulate(@control_u1, x0, tspan, m1, m2, L1, L2, g, Tmax, k);

    % u2 system
    [t2, x2, E2] = simulate(@control_u2, x0, tspan, m1, m2, L1, L2, g, Tmax, k);

    % plot energy!!
    figure;
    semilogy(t1, E1);
    hold on;
    semilogy(t2, E2);
    xlabel('Time (s)');
    ylabel('Total Energy (J)');
    legend('u1 system', 'u2 system');
    title('Energy vs Time');
end

% u1 controller
function u = control_u1(dtheta1, dtheta2, k, Tmax)
    u1 = -k * dtheta1;
    u = [saturate(u1, Tmax); 0];
end

% u2 controller
function u = control_u2(dtheta1, dtheta2, k, Tmax)
    u2 = -k * dtheta2;
    u = [0; saturate(u2, Tmax)];
end

% sat fnction
function y = saturate(x, limit)
    y = max(-limit, min(limit, x));
end

% sim wrapper function
function [t, x, E] = simulate(controller, x0, tspan, m1, m2, L1, L2, g, Tmax, k)
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t, x] = ode45(@(t, x) dynamics(t, x, controller, m1, m2, L1, L2, g, Tmax, k), tspan, x0, opts);
    % energy calculation
    E = zeros(length(t),1);
    for i = 1:length(t)
        E(i) = total_energy(x(i,:)', m1, m2, L1, L2, g);
    end
end

% double pendulum dynamics
function dx = dynamics(~, x, controller, m1, m2, L1, L2, g, Tmax, k)
    theta1 = x(1); theta2 = x(2);
    dtheta1 = x(3); dtheta2 = x(4);
    M11 = (m1 + m2)*L1^2 + m2*L2^2 + 2*m2*L1*L2*cos(theta2);
    M12 = m2*L2^2 + m2*L1*L2*cos(theta2);
    M22 = m2*L2^2;
    M = [M11 M12; M12 M22];
    C = m2*L1*L2*sin(theta2)*[-dtheta2 - dtheta1, -dtheta2; dtheta1, 0];

    % grad of U
    dU1 = (m1 + m2)*g*L1*sin(theta1) + m2*g*L2*sin(theta1 + theta2);
    dU2 = m2*g*L2*sin(theta1 + theta2);
    dU = [dU1; dU2];

    % controler input
    u = controller(dtheta1, dtheta2, k, Tmax);

    % Compute acceleration
    ddtheta = M \ (u - C*[dtheta1; dtheta2] - dU);

    dx = [dtheta1; dtheta2; ddtheta];
end

% total energy
function E = total_energy(x, m1, m2, L1, L2, g)
    theta1 = x(1); theta2 = x(2);
    dtheta = x(3:4);
    M11 = (m1 + m2)*L1^2 + m2*L2^2 + 2*m2*L1*L2*cos(theta2);
    M12 = m2*L2^2 + m2*L1*L2*cos(theta2);
    M22 = m2*L2^2;
    M = [M11 M12; M12 M22];

    % KE
    KE = 0.5 * dtheta' * M * dtheta;

    % PE
    PE = (m1 + m2)*g*L1*(1 - cos(theta1)) + m2*g*L2*(1 - cos(theta1 + theta2));

    E = KE + PE;
end
