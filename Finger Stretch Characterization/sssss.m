s = tf('s');

% G_r(s) = -0.0184*(s + 0.0068)/[(s+0.2647)*(s+0.0063)]
Gr = -0.0184*(s+0.0068)/((s+0.2647)*(s+0.0063));

final_value = dcgain(Gr);
fprintf('Final Value: %.5f\n', final_value);

figure;
step(Gr);
grid on;
title('Step resppomse of G_r(s)');

info = stepinfo(Gr, 'SettlingTimeThreshold', 0.01);
fprintf('1%% Settling Time: %.2f sec\n', info.SettlingTime);
