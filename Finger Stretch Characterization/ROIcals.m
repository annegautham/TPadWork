clear; clf; clc;

% equilibrium
xe = [2; 6];

% grid
x1 = linspace(-10,10,600);
x2 = linspace(-10,10,600);
[X1,X2] = meshgrid(x1,x2);
DX1 = X1 - xe(1);
DX2 = X2 - xe(2);

% Lyapunov function V = (x - xe)' * I * (x - xe)
Vgrid = DX1.^2 + DX2.^2;
dVdx1 = 2*DX1;
dVdx2 = 2*DX2;

bestArea = 0;
bestK = 0;

% try various k values
ks = linspace(-10, 10, 200);  % sweep over feedbacks
for k = ks
    % control law u = -k*(x2 - 6)
    U = -k * (X2 - xe(2));

    % dynamics with feedback
    F1 = X1 - X1.^3 + X2;
    F2 = 3*X1 - X2 + U;

    % Lie derivative
    LfVgrid = dVdx1 .* F1 + dVdx2 .* F2;

    % largest c with LfV < 0 inside V <= c
    mask = (LfVgrid < 0);
    if any(mask(:))
        c = max(Vgrid(mask));
        area = pi * c;  % since P = I
        if area > bestArea
            bestArea = area;
            bestK = k;
            bestC = c;
            bestLfV = LfVgrid;
        end
    end
end

fprintf('Max area ≈ %.1f at k = %.3f\n', bestArea, bestK);

% final plot
figure(1); hold on; axis equal
colormap([0.8 0.95 1; 1 1 1; 1 0.8 0.8]);
%plotting boundaries and roA

contourf(X1, X2, bestLfV, [-1e6 0]);
contourf(X1, X2, Vgrid, [0 bestC]);

contour(X1, X2, bestLfV, [0 0]);
contour(X1, X2, Vgrid, [bestC bestC]);
plot(xe(1), xe(2), 'ko', 'MarkerFaceColor', 'k');

xlabel('x_1'); ylabel('x_2');
title(sprintf('Best \\Omega_c with feedback k = %.3f (area ≈ %.0f)', bestK, bestArea));
legend({'Vdot<0','\Omega_c','Vdot=0','V=c','equilibrium'}, 'Location','northwest');
xlim([-10 10]); ylim([-10 10]);
