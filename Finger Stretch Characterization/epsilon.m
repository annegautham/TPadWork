epsilon_min = 1e-6;
    epsilon_max = 1e3;
    epsilonVec = zeros(1, length(freq_lin));
    for j = 1:length(freq_lin)
        f_val = freq(j);
        if f_val < 20
            epsilonVec(j) = epsilon_max;
        elseif f_val >= 20 && f_val < 40
            epsilonVec(j) = epsilon_max - (f_val-20)/(40-20)*(epsilon_max - epsilon_min);
        elseif f_val >= 40 && f_val <= 340
            epsilonVec(j) = epsilon_min;
        elseif f_val > 340 && f_val <= 360
            epsilonVec(j) = epsilon_min + (f_val-340)/(360-340)*(epsilon_max - epsilon_min);
        else
            epsilonVec(j) = epsilon_max;
        end
    end
