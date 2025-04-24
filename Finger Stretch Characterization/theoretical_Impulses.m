% Beam + system parameters
E = 230e9;
r = 0.5e-3;
I = pi*r^4/4;
rho = 1600;
L = 0.12; % 100 x
x_meas = 0.055; % 55 mm
m_point = 0.0244; % 20 g

% beam mass + stiffness at ldv point
A = pi*r^2;
m_beam = rho * A * L;
m_eff = 0.21 * m_beam + m_point;
k_eff = 3*E*I / x_meas^3;

% damping
zeta = 0.06;
omega_n = sqrt(k_eff/m_eff);
omega_d = omega_n * sqrt(1 - zeta^2);
c_eff = 2*zeta*sqrt(k_eff*m_eff);

% Time vector: 10 seconds, 10k points (1 kHz sample rate)
t = linspace(0, 2.5, 10000);

% Impulse response (unit impulse at t=0)
h = (1/(m_eff*omega_d)) * exp(-zeta*omega_n*t) .* sin(omega_d*t);

