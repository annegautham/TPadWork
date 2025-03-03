G = tf(10, [1 11 10 0]);
rlocus(G);
plead = 22.5;
zlead = 0.9;
Klead = 114;
Dlead = Klead * tf([1 zlead], [1 plead]);
plag = 1000;
zlag = 1;
Klag = 1000.74;
Dlag = Klag * tf([1 zlag], [1 plag]);

L = Dlead*G*Dlag;

rlocus(L);

T = feedback(L,1);

figure;
step(T)
info = stepinfo(T, 'SettlingTimeThreshold', 0.05);
fprintf('1%% Settling Time: %.2f sec\n', info.SettlingTime);