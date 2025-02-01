clear; close all; clc;

K_values = [0.5, 1.0, 2.0, 6.5e-2];
t = 0:0.01:150;
numPlots = length(K_values);

figure;
sgtitle('Step Response for Different K Values');

for i = 1:numPlots
    K = K_values(i);
    wn = sqrt(0.2 * K);
   
    num = wn^2;
    den = [1, 0.104, wn^2];
    sys = tf(num, den);
    
    % Compute step response
    [y, t_out] = step(sys, t);
    
    % Compute overshoot
    zeta = 0.104 / (2 * wn);
    if zeta < 1
        Mp = (max(y) - 1) * 100;
        overshootText = sprintf('Max Overshoot: %.2f%%', Mp);
    else
        overshootText = 'No Overshoot';
    end
    
    % rise time (from 10% to 90%)
    y_final = y(end);
    idx_10 = find(y >= 0.1 * y_final, 1);
    idx_90 = find(y >= 0.9 * y_final, 1);
    tr = t_out(idx_90) - t_out(idx_10);
    risetimeText = sprintf('Rise Time: %.2f sec', tr);
    
    %plot
    subplot(2, 2, i);
    plot(t_out, y, 'b', 'LineWidth', 1.5);
    grid on;
    title(sprintf('K = %.4f', K));
    xlabel('Time (s)');
    ylabel('Response');
    
    %show overshoot
    text(0.05 * max(t), 0.85 * max(y), overshootText, 'FontSize', 10);
    text(0.05 * max(t), 0.75 * max(y), risetimeText, 'FontSize', 10);
end
