s = tf('s');
K= 600;
Kt = 1;
M0 = 0.6;
E = tf([K], [1 10]);
G = tf([1 3], [1 4 5]);
P = tf([M0], [1 0]);

theta = E* G* 1/s / (1+(1+Kt*s)*E*G*1/s);

newktsutfiio = E*G/ (2+E*G*1/s);
rlocus(newktsutfiio);