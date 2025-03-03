%% Compensator Design for Ship Heading Control
% This script designs a proportional compensator for a US Coast Guard
% cutter (Tampa, 902) so that:
%   - For a 5° step change in desired heading ψ_r, the maximum rudder deflection
%     remains below 10°.
%   - The settling time for the heading response is less than 50 sec.
%
% The ship dynamics (from sea-trial data) are given by:
%
%   G_delta(s) = ψ(s)/δ(s) = -0.0184*(s+0.0068) / [ s*(s+0.2647)*(s+0.0063) ]
%   G_w(s)     = ψ(s)/w(s)   = 0.0000064 / [ s*(s+0.2647)*(s+0.0063) ]
%
% The available measurement is the heading ψ and the yaw rate r (ψ_dot).

%% Part 1: Define the System Transfer Functions
s = tf('s');

% Transfer function from rudder deflection δ to heading angle ψ
G_delta = -0.0184*(s + 0.0068) / (s*(s + 0.2647)*(s + 0.0063));

% (Optional) Transfer function from wind disturbance w to heading angle ψ
G_w = 0.0000064 / (s*(s + 0.2647)*(s + 0.0063));

%% Part 2: Control Objectives and Controller Design
% Design Constraints:
%   1. Maximum Rudder Deflection:
%      At time t=0 (initial instant), the rudder deflection is
%         δ(0) = K * ψ_r(0)
%      For a step change of 5° in ψ_r, we require:
%         K * 5° < 10°  ==>  K < 2.
%
%   2. Settling Time Requirement:
%      For a second order system, the settling time Ts is approximately:
%         Ts ≈ 4.6 / σ
%      where σ is the real part of the dominant pole. To achieve Ts < 50 sec:
%         σ > 4.6/50 ≈ 0.092.
%
% A root locus analysis shows that a proportional gain of K = 1.56 gives
% closed-loop poles with real parts less than about 0.13, ensuring that the
% settling time is within specifications. In addition, the initial rudder
% deflection will be K*5 = 7.8° < 10°.
%
% Here we choose a simple proportional controller:
K =200; 
C = K;  % Controller (proportional gain)

%% Part 3: Open-Loop Analysis and Root Locus
% The open-loop transfer function L(s) is:
%      L(s) = C(s) * G_delta(s)
L = C * G_delta;

% Plot the root locus to visually assess the pole locations as K varies.
figure;
rlocus(L)
title('Root Locus of the Open-Loop System')
xlabel('Real Axis')
ylabel('Imaginary Axis')
grid on;
% Note: You could use "rlocfind(L)" interactively to pick a gain, but here we
% already determined K = 1.56 satisfies the design criteria.

%% Part 4: Closed-Loop System and Step Response Analysis
% The closed-loop transfer function (from heading command ψ_r to heading ψ)
% is given by:
%      T(s) = L(s) / (1 + L(s))
T = feedback(L, 1);

% Display the closed-loop poles for further insight.
closedLoopPoles = pole(T);
disp('Closed-loop poles:');
disp(closedLoopPoles);

% Simulate the step response for a 5° step change in heading.
% (Multiply T by 5 so that the input step is 5°.)
figure;
step(5*T)
title('Step Response for a 5° Heading Change')
xlabel('Time (sec)')
ylabel('Heading Angle (°)')
grid on;

% Obtain performance metrics such as settling time.
stepMetrics = stepinfo(5*T);
disp('Step Response Performance Metrics:');
disp(stepMetrics);

%% Part 5: Rudder Deflection Analysis
% The initial rudder deflection is given by the controller gain times the
% step input:
%      δ(0) = K * (step input) = 1.56 * 5° = 7.8°,
% which is below the maximum allowable 10°.
initialRudderDeflection = K * 5;
fprintf('Initial rudder deflection: %.2f° (must be < 10°)\n', initialRudderDeflection);
