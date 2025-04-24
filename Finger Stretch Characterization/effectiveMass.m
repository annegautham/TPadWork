clc; clear;

% -------------------------------
% Beam and Material Parameters
% -------------------------------
L = 0.1;                 % Beam length (m) = 100 mm
x_meas = 0.06;           % Measurement point (m) = 60 mm
rho = 1600;              % Density (kg/m^3)
r = 0.5e-3;              % Radius (m)
A = pi * r^2;            % Cross-sectional area (m^2)
m_beam = rho * A * L;    % Total mass of the beam (kg)

% -------------------------------
% Mode Shape Parameters (1st Mode)
% -------------------------------
betaL = 1.87510407;               % First eigenvalue for cantilever beam
beta = betaL / L;
gamma = (cos(betaL) + cosh(betaL)) / ...
        (sin(betaL) + sinh(betaL));  % Mode shape coefficient

% -------------------------------
% Define Mode Shape Function
% -------------------------------
% Mode shape function for 1st mode
phi_raw = @(x) (sinh(beta * x) - sin(beta * x)) / ...
               (sinh(beta * L) - sin(beta * L));

% -------------------------------
% Compute Integral of Mode Shape Squared
% -------------------------------
% Compute the mode shape at each point
x_vals = linspace(0, L, 1000);
phi_vals = phi_raw(x_vals);

% Integral of mode shape squared
denominator = trapz(x_vals, phi_vals.^2);  % ∫φ(x)^2 dx

% -------------------------------
% Compute Effective Mass
% -------------------------------
% Mode shape at measurement point
phi_x = phi_raw(x_meas);                % Mode shape at measurement point
mass_participation = phi_x^2 / denominator;     % Mass participation factor
m_eff = mass_participation * m_beam;            % Effective mass in kg

% -------------------------------
% Display Results
% -------------------------------
fprintf('--- Effective Mass Calculation ---\n');
fprintf('Beam length            = %.1f mm\n', L*1000);
fprintf('Measurement location   = %.1f mm\n', x_meas*1000);
fprintf('Beam radius            = %.3f mm\n', r*1000);
fprintf('Beam density           = %.1f kg/m³\n', rho);
fprintf('Total beam mass        = %.6f kg (%.3f g)\n', m_beam, m_beam*1000);
fprintf('phi(x)                 = %.4f\n', phi_x);
fprintf('Integral of phi^2 dx   = %.4f\n', denominator);
fprintf('Effective mass ratio   = %.4f\n', mass_participation);
fprintf('Effective mass at %.1f mm = %.6f kg (%.3f g)\n', ...
        x_meas*1000, m_eff, m_eff*1000);
