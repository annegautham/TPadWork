% === PARAMETERS ===
Fs = 10000;
Ts = 1 / Fs;  % 10 kHz sampling time



% % === CREATE IDDATA OBJECT ===
% data = iddata(pos_sync, in_sync, Ts);
% 
% % === ESTIMATE IMPULSE RESPONSE ===
% model_ir = impulseest(data);
% 
% figure;
% impulseplot(model_ir);
% title('Estimated Impulse Response');
% 
% % === OPTIONAL: FIT 2ND-ORDER TRANSFER FUNCTION ===
% model_tf = tfest(data, 2, 0);  % 2 poles, 0 zeros
% figure;
% impulseplot(model_tf);
% title('2nd-Order Transfer Function Fit');
