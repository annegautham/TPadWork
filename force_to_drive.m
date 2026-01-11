load('interpolatedImpedance_data.mat');
[~, ~, i_t] = chirpWithEnvelope(30, 350, 0.8, 5, 10000, ...
                                freq_interp, ohm_interp, true);
figure;

BL = 2.7;
in = BL*i_t;
%% 

load('18g_30-350.mat');


Fs = 10000;
Ts = 1/Fs;

raw = outputVec{1};
raw = raw(:)';
vel = (raw - mean(raw))*125;
vel = lowpass(vel, 1000, Fs);
b = fir1(500, 25/(Fs/2), 'high');
vel = filter(b, 1, vel);
pos = cumsum(vel)*Ts;
pos = pos - mean(pos);
[cc, lags] = xcorr(pos, in, 'none');
[~, max_idx] = max(abs(cc));
lag = lags(max_idx);
if lag >= 0
    pos = pos(lag+1:end);
end


N_sync = min(length(in), length(pos));
in = in(1:N_sync);
pos = pos(1:N_sync);
inv = fliplr(in');
ir = conv(inv, pos);
ir = ir(50000:end);
ir_trimmed = ir(1:5000);

data = iddata(ir_trimmed(:), [1; zeros(length(ir_trimmed)-1,1)], Ts);
sys_est = tfest(data, 2);
[ir_est, t_ir] = impulse(sys_est, length(ir_trimmed)*Ts);
plot(t_ir, ir_trimmed(1:length(t_ir)), 'LineWidth', 1.5)
hold on
plot(t_ir, ir_est, '--', 'LineWidth', 1.5)
legend('Measured IR', 'Fitted IR')
