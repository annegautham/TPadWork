load('fingerpad_session_current.mat');
cond = sessionData.conditions(2);
trials = cond.trials;

fs = trials(1).fs;
dur = trials(1).duration;
t_uniform = (0:1/fs:dur-1/fs)';


vel_mat = zeros(length(t_uniform), 4);
force_mat = zeros(length(t_uniform), 4);

forceFactor = -2.7;

for k = 1:4
    b = fir1(300, 20/(fs/2), 'high');

    vel = (trials(k).ldvSignal(:) - mean(trials(k).ldvSignal(:)))*125;
    vel = lowpass(vel, 1000, fs);
    vel = filter(b, 1, vel);

    curr = interp1(trials(k).inaTime(:), trials(k).inaCurrent(:), t_uniform, 'linear', 'extrap');
    curr = curr - mean(curr);
    force = curr * forceFactor;
    force = filter(b, 1, force);

    force_mat(:,k) = force;
    vel_mat(:,k) = vel;
end



figure;
plot(t_uniform, vel_mat, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Velocity (a.u.)');
title('LDV Trials After Mean Removal + Interpolation');
legend('Trial 1','Trial 2','Trial 3','Trial 4');

figure;
plot(t_uniform, force_mat, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Current (A, mean removed)');
title('INA219 Current Trials After Mean Removal + Interpolation');
legend('Trial 1','Trial 2','Trial 3','Trial 4');


window = 1024;
noverlap = 512;
nfft = 2048;

for k = 1:4
    % rev = flip(force_mat(:,k));
    % ir = conv(rev, vel_mat(:,k));
    % figure;
    % hold on;
    % plot(ir);
    vel = vel_mat(:,k);
    force = force_mat(:,k);
    N = length(vel);

    figure;
    spectrogram(vel, window, noverlap, nfft, fs, 'yaxis');
    title(sprintf('LDV Trial %d – Spectrogram', k));
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');
    colorbar;

    V = fft(vel);
    F = fft(force);   % or coil current * force constant
    f = (0:N-1)' * (fs/N);
    half = 1:floor(N/2);
    fmask = (f >= 30 & f <= 350);
    f_sel = f(fmask);

    V = V(fmask);
    F = F(fmask);
    H = V./ F;
    figure;
    subplot(2,1,1);
    plot(f_sel, abs(V));
    xlabel('Frequency (Hz)');
    ylabel('|V(ω)|');
    title('Velocity Spectrum');
    
    subplot(2,1,2);
    plot(f_sel, abs(F));
    xlabel('Frequency (Hz)');
    ylabel('|F(ω)|');
    title('Force / Current Spectrum');

    figure;
    plot(f_sel, abs(H));
    xlabel('Frequency (Hz)');
    ylabel('|H(ω)| = |V/F|');
    title('Measured FRF (Admittance) – Trial 1');

    H_mat = zeros(length(f_sel), 4);
    H_mat(:,k) = H;

end

H_mean = mean(H_mat, 2);
figure;
subplot(2,1,1);
plot(f_sel, abs(H_mean), 'b', 'LineWidth', 1.4);
xlabel('Frequency (Hz)');
ylabel('|H(\omega)|');
title('Mean FRF Magnitude (4 trials)');
grid on;

subplot(2,1,2);
plot(f_sel, unwrap(angle(H_mean)), 'r', 'LineWidth', 1.4);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
title('Mean FRF Phase (4 trials)');
grid on;


% %%
% 
% f = f_sel(:);                 % Hz
% omega = 2*pi*f;               % rad/s
% H = H_mean(:);                % admittance: (mm/s)/N
% 
% % --- Compute impedance ------------------------------------
% Z = 1 ./ H;                   % N·s/mm (complex)
% 
% % --- Extract damping and stiffness -------------------------
% b_w = real(Z);                % damping b(ω) [N·s/mm]
% k_w = - omega .* imag(Z);     % stiffness k(ω) [N/mm]
% 
% 
% %% --- 1. Constant Kelvin-Voigt (KV) Model Fit -----------------
% % Model: Z_KV(jω) = b_KV - j(k_KV/ω) (assuming m=0)
% 
% % Damping (b_KV): Average of the real part
% b_KV = mean(real(Z));
% 
% % Stiffness (k_KV): Linear fit of Im{Z} = k_KV * (-1/ω)
% X_KV = -1 ./ omega;
% Y_KV = imag(Z);
% % polyfit returns coefficients [slope, intercept]
% % Here, the slope is the stiffness k_KV, and the intercept (p_KV(2)) should be zero.
% p_KV = polyfit(X_KV, Y_KV, 1);
% k_KV = p_KV(1); % N/mm
% 
% Z_KV_model = b_KV + 1j * (-k_KV ./ omega);
% 
% fprintf('--- Constant Kelvin-Voigt (m=0) Fit Results -------------------\n');
% fprintf('Model: Z(jω) = b - j(k/ω)\n');
% fprintf('Damping (b_KV): %.4f N*s/mm\n', b_KV);
% fprintf('Stiffness (k_KV): %.4f N/mm\n', k_KV);
% fprintf('R-squared (Approximate Linear Fit): %.2f\n', 1 - sum((Y_KV - (k_KV*X_KV)).^2) / sum((Y_KV - mean(Y_KV)).^2));
% fprintf('\n');
% 
% 
% %% --- 2. Wiertlewski & Hayward (MSD) Model Fit -----------------
% % Model: Z_MSD(jω) = b + j(m*ω - k/ω)
% % This is the two-step method used in the reference paper.
% 
% % Step 1: Damping (b_MSD) is the average of the real part
% b_MSD = mean(real(Z));
% 
% % Step 2: Stiffness (k_MSD) and Mass (m_MSD) from non-linear least squares on the imaginary part.
% % x(1) = m_MSD, x(2) = k_MSD
% 
% % Objective function for least-squares (residuals vector of the imaginary part)
% ImagZ_fit_function = @(x) imag(Z) - (x(1) * omega - x(2) ./ omega);
% 
% % Initial Guess [m_init, k_init]
% % Wiertlewski & Hayward reported mass m in the range of 110-230 mg (~1e-7 N*s^2/mm).
% % Stiffness k in range of 0.6-2.1 N/mm.
% x0 = [1e-7, 1];
% 
% % Use lsqnonlin (requires Optimization Toolbox)
% options = optimset('Display', 'off');
% [x_fit, resnorm] = lsqnonlin(ImagZ_fit_function, x0, [], [], options);
% 
% m_MSD_N_s2_mm = x_fit(1); % N*s^c2/mm
% k_MSD = x_fit(2);         % N/mm
% 
% % Conversion for reporting (1 N*s^2/mm = 1000 kg = 1,000,000 g)
% m_MSD_g = m_MSD_N_s2_mm * 1e6;
% m_MSD_mg = m_MSD_g * 1000;
% 
% Z_MSD_model = b_MSD + 1j * (m_MSD_N_s2_mm * omega - k_MSD ./ omega);
% 
% fprintf('--- Wiertlewski & Hayward (MSD) Fit Results ------------------\n');
% fprintf('Model: Z(jω) = b + j(mω - k/ω)\n');
% fprintf('Damping (b_MSD): %.4f N*s/mm\n', b_MSD);
% fprintf('Stiffness (k_MSD): %.4f N/mm\n', k_MSD);
% fprintf('Mass (m_MSD): %.8f N*s^2/mm (%.2f g / %.2f mg)\n', m_MSD_N_s2_mm, m_MSD_g, m_MSD_mg);
% fprintf('\n');
% 
% 
% %% --- Plots of Fitted Models --------------------------------
% figure;
% subplot(2,1,1);
% plot(f, real(Z), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Measured Re\{Z\}');
% hold on;
% plot(f, b_MSD * ones(size(f)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Fitted');
