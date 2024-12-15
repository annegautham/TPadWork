data = readtable('TEAX19C01-8_impedance', 'Format', '%f %f %f', 'HeaderLines', 1); 

freq = data{:, 1};
ohm = data{:, 2};
phase = data{:,3};
cuttoff = 2000;
range = freq <= cuttoff;
freq = freq(range);
ohm = ohm(range);
phase = phase(range);

freq_interp = linspace(min(freq), max(freq), 10^3);
ohm_interp = interp1(freq, ohm, freq_interp, 'spline');
phase_interp = interp1(freq, phase, freq_interp, 'spline');
save('interpolatedImpedance_data.mat', 'freq_interp', 'ohm_interp');

figure;
plot(freq, ohm, 'o', freq_interp, ohm_interp, '-');
xlabel('Frequency (Hz)');
ylabel('Impedance (Ohms)');
title('Impedance vs. Frequency');
legend('OG Data', 'Interp');
xticks(0:100:cuttoff); 
yticks(min(ohm):10:max(ohm));

figure;
plot(freq, phase, 'o', freq_interp, phase_interp, '-');
xlabel('Frequency (Hz)');
ylabel('Phase (Rad)');
title('Phase vs. Frequency');
legend('OG Data', 'Interp');
xticks(0:100:cuttoff); 
yticks(min(phase):10:max(phase));