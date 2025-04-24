% 5-2 Part b
x_spring = 4;
f_spring = 5;
f_damper = -f_spring;
v_damper = f_damper/(3/2);
v_spring = v_damper;
fprintf('Displacement of spring (x_spring): %.2f mm\n', x_spring);
fprintf('Force in spring (f_spring): %.2f N\n', f_spring);
fprintf('Force in damper (f_damper): %.2f N\n', f_damper);
fprintf('Velocity in damper (v_damper): %.2f mm/s\n', v_damper);
fprintf('Velocity in spring (v_spring): %.2f mm/s\n', v_spring);
% 5-2 part c
xs = -6.4:0.1:8;
fs = zeros(size(xs));
for i = 1:length(xs)
x = xs(i);
if x >= -4 && x <= 4
fs(i) = 1.25 * xs(i);
elseif x > 4
fs(i) = 0.5*(xs(i) - 4) + 5;
else
fs(i) = NaN;
end
end
vd = -16/3:0.1:16;
fd = zeros(size(vd));
for i = 1:length(vd)
v = vd(i);
if v <=0
fd(i) = (3/2)*(v);
else
fd(i) = 0.5*(v);
end
end
figure;
subplot(1,2,1);
plot(xs, fs, "b", 'LineWidth', 2);
xlabel('x_s (mm)');
ylabel('f_s (N)');
title('Spring Force vs Displacement');
grid on;
subplot(1,2,2);
plot(vd, fd, 'r', 'LineWidth', 2);
xlabel('v_d (mm/s)');
ylabel('f_d (N)');
title('Damper Force vs Velocity');
grid on;
% 5-2 part d
F_range = -8:1:8;
x_s = zeros(size(F_range));
v_s = zeros(size(F_range));
for i = 1:length(F_range)
F = F_range(i);
if F >= -5 && F <= 5
x = F/1.25;
elseif F > 5
x = (F - 3) / 0.5;
else
x = NaN;
end
x_s(i) = x;
f_spring = F;
f_damper = -f_spring;
if f_damper >= 0
b_d = 0.5;
else
b_d = 1.5;
end
v = f_damper / b_d;
v_s(i) = v;
end
figure;
plot(x_s, v_s, 'LineWidth', 2);
xlabel('x_s (mm)');
ylabel("v_s (mm/s)");
title("Spring Velocity vs Displacement");
grid on;