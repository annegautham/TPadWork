mass = [0.018; 0.0244; 0.027; 0.0284];
theoretical_fn = [21.77; 18.7; 17.77; 17.33];
predicted_fn = [22.673; 19.472; 18.799; 18.126];
r = [0.9951; 0.9812; 0.9256; 0.9657];

% Create table
resultsTable = table(mass, theoretical_fn, predicted_fn, r,  ...
    'VariableNames', {'Mass_kg', 'Theoretical_fn_Hz', 'Predicted_fn_Hz', 'r'});

% Display table
disp(resultsTable);