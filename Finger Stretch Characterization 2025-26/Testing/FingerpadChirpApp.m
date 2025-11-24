classdef FingerpadChirpApp < matlab.apps.AppBase
    % FingerpadChirpApp
    % 175 conditions GUI. Each condition has 4 trials.
    % Teensy generates chirp, INA219 current is returned and stored.

    properties (Access = public)
        UIFigure matlab.ui.Figure
        TrialTable matlab.ui.control.Table

        PortLabel matlab.ui.control.Label
        PortField matlab.ui.control.EditField
        ConnectButton matlab.ui.control.Button

        RunTrialButton matlab.ui.control.Button
        AcceptButton matlab.ui.control.Button
        RetakeButton matlab.ui.control.Button

        AutoAdvanceCheckbox matlab.ui.control.CheckBox

        StatusLabel matlab.ui.control.Label
        SelectedConditionLabel matlab.ui.control.Label
        TrialInfoLabel matlab.ui.control.Label

        ChirpAxes matlab.ui.control.UIAxes
        INAaxes matlab.ui.control.UIAxes

        F1Label matlab.ui.control.Label
        F2Label matlab.ui.control.Label
        DurLabel matlab.ui.control.Label
        EnvLabel matlab.ui.control.Label
        FsLabel matlab.ui.control.Label

        F1Field matlab.ui.control.NumericEditField
        F2Field matlab.ui.control.NumericEditField
        DurField matlab.ui.control.NumericEditField
        EnvField matlab.ui.control.NumericEditField
        FsField matlab.ui.control.NumericEditField
    end

    properties (Access = private)
        Teensy              % serialport object
        Conditions          % table of 175 conditions
        SessionData         % struct storing all conditions + trials
        SessionFileName     % .mat session file
        SelectedRowIdx      % which condition row is selected
        CurrentTrialIdx     % trial 1..4 for the selected condition

        % --- INA219 data from last run ---
        LastINATime    double = []  % seconds
        LastINACurrent double = []  % amps
        LastHXTime     double = []
        LastHXValue    double = []
    end

    %======================================================================
    % PRIVATE METHODS
    %======================================================================
    methods (Access = private)

        function startupFcn(app)
            % Build the 175-condition grid: 7×5×5
            orientations = [0 15 30 45 60 75 90];
            indentations = [0.4 0.6 0.8 1.0 1.2];
            apertures    = [20 40 80 120 160];

            rows = [];
            for o = orientations
                for d = indentations
                    for a = apertures
                        rows = [rows; o, d, a]; %#ok<AGROW>
                    end
                end
            end

            nCond = size(rows,1);  % 175

            app.Conditions = table( ...
                rows(:,1), rows(:,2), rows(:,3), ...
                zeros(nCond,1), ...  % TrialsCompleted
                'VariableNames', {'Orientation','Indentation','Aperture','TrialsCompleted'});

            app.TrialTable.Data = app.Conditions;
            app.SelectedRowIdx  = [];
            app.CurrentTrialIdx = 1;

            % Initialize sessionData struct
            app.SessionData = struct();
            app.SessionData.meta.created     = datestr(now);
            app.SessionData.meta.description = 'Fingerpad dynamic experiments';
            app.SessionData.conditions       = struct([]);

            % Pre-allocate conditions struct
            app.SessionData.conditions = repmat( ...
                struct('orientation',[], 'indentation',[], 'aperture',[], ...
                'trialsCompleted',0, 'trials',struct([])), ...
                nCond, 1);

            for i = 1:nCond
                app.SessionData.conditions(i).orientation = app.Conditions.Orientation(i);
                app.SessionData.conditions(i).indentation = app.Conditions.Indentation(i);
                app.SessionData.conditions(i).aperture    = app.Conditions.Aperture(i);
            end

            % Session file name
            timestamp = datestr(now,'yyyymmdd_HHMMSS');
            app.SessionFileName = sprintf('fingerpad_session_%s.mat', timestamp);
            app.StatusLabel.Text = sprintf('Session file: %s', app.SessionFileName);
        end

        function createComponents(app)

            % Main figure
            app.UIFigure = uifigure('Name','Fingerpad Chirp Experiment');
            app.UIFigure.Position = [100 100 1350 700];
            app.UIFigure.Color = [0.98 0.98 0.98];

            % ============================================================
            % LEFT PANEL — Trial Table
            % ============================================================
            leftPanel = uipanel(app.UIFigure);
            leftPanel.Title = 'Conditions (175)';
            leftPanel.FontWeight = 'bold';
            leftPanel.Position = [20 80 520 600];

            % Trial table inside left panel
            app.TrialTable = uitable(leftPanel, ...
                'ColumnName', {'Orientation','Indentation','Aperture','TrialsCompleted'}, ...
                'ColumnEditable', [false false false false], ...
                'ColumnFormat', {'numeric','numeric','numeric','numeric'});
            app.TrialTable.Position = [10 10 500 550];
            app.TrialTable.SelectionType  = 'row';
            app.TrialTable.FontSize = 13;
            app.TrialTable.CellSelectionCallback = @(src,evt) app.onConditionSelected(evt);

            % Make columns sortable (R2021b+)
            if isprop(app.TrialTable,'ColumnSortable')
                app.TrialTable.ColumnSortable = true(1,4);
            end

            % ============================================================
            % RIGHT-TOP PANEL — Connection + Chirp Parameters
            % ============================================================
            topPanel = uipanel(app.UIFigure);
            topPanel.Title = 'Connection + Chirp Parameters';
            topPanel.FontWeight = 'bold';
            topPanel.Position = [560 500 760 180];

            % ---- Port input
            app.PortLabel = uilabel(topPanel);
            app.PortLabel.Position = [20 130 40 22];
            app.PortLabel.Text = 'Port:';

            app.PortField = uieditfield(topPanel,'text');
            app.PortField.Position = [60 128 100 25];
            app.PortField.Value = 'COM4';

            app.ConnectButton = uibutton(topPanel,'push');
            app.ConnectButton.Position = [170 126 130 30];
            app.ConnectButton.Text = 'Connect Teensy';
            app.ConnectButton.FontWeight = 'bold';
            app.ConnectButton.ButtonPushedFcn = @(btn,evt) app.onConnectTeensy();

            % ---- Chirp param labels
            xpos = 20; ypos = 80;

            app.F1Label = uilabel(topPanel,'Text','f1 (Hz)');
            app.F1Label.Position = [xpos ypos 60 22];

            app.F1Field = uieditfield(topPanel,'numeric');
            app.F1Field.Position = [xpos+70 ypos 70 24];
            app.F1Field.Value = 30;

            app.F2Label = uilabel(topPanel,'Text','f2 (Hz)');
            app.F2Label.Position = [xpos+160 ypos 60 22];

            app.F2Field = uieditfield(topPanel,'numeric');
            app.F2Field.Position = [xpos+230 ypos 70 24];
            app.F2Field.Value = 350;

            app.DurLabel = uilabel(topPanel,'Text','dur (s)');
            app.DurLabel.Position = [xpos+330 ypos 60 22];

            app.DurField = uieditfield(topPanel,'numeric');
            app.DurField.Position = [xpos+400 ypos 70 24];
            app.DurField.Value = 5;

            app.EnvLabel = uilabel(topPanel,'Text','env (0-1)');
            app.EnvLabel.Position = [xpos ypos-40 60 22];

            app.EnvField = uieditfield(topPanel,'numeric');
            app.EnvField.Position = [xpos+70 ypos-40 70 24];
            app.EnvField.Value = 0.05;

            app.FsLabel = uilabel(topPanel,'Text','fs (Hz)');
            app.FsLabel.Position = [xpos+160 ypos-40 60 22];

            app.FsField = uieditfield(topPanel,'numeric');
            app.FsField.Position = [xpos+230 ypos-40 90 24];
            app.FsField.Value = 10000;

            % Auto-advance checkbox
            app.AutoAdvanceCheckbox = uicheckbox(topPanel);
            app.AutoAdvanceCheckbox.Position = [xpos+350 ypos-40 220 24];
            app.AutoAdvanceCheckbox.Text = 'Auto-advance to next trial';
            app.AutoAdvanceCheckbox.Value = true;
            app.AutoAdvanceCheckbox.FontSize = 13;

            % ============================================================
            % RIGHT-MIDDLE PANEL — Trial Control
            % ============================================================
            midPanel = uipanel(app.UIFigure);
            midPanel.Title = 'Trial Control';
            midPanel.FontWeight = 'bold';
            midPanel.Position = [560 420 760 80];

            app.SelectedConditionLabel = uilabel(midPanel);
            app.SelectedConditionLabel.Position = [20 40 650 22];
            app.SelectedConditionLabel.Text = 'Selected condition: none';

            app.TrialInfoLabel = uilabel(midPanel);
            app.TrialInfoLabel.Position = [20 15 300 22];
            app.TrialInfoLabel.Text = 'Trial: -';

            % Buttons
            app.RunTrialButton = uibutton(midPanel,'push');
            app.RunTrialButton.Position = [350 20 120 35];
            app.RunTrialButton.Text = 'Run Trial';
            app.RunTrialButton.FontWeight = 'bold';
            app.RunTrialButton.ButtonPushedFcn = @(btn,evt) app.onRunTrial();

            app.AcceptButton = uibutton(midPanel,'push');
            app.AcceptButton.Position = [480 20 100 35];
            app.AcceptButton.Text = 'Accept';
            app.AcceptButton.Enable = 'off';
            app.AcceptButton.ButtonPushedFcn = @(btn,evt) app.onAcceptTrial();

            app.RetakeButton = uibutton(midPanel,'push');
            app.RetakeButton.Position = [590 20 100 35];
            app.RetakeButton.Text = 'Retake';
            app.RetakeButton.Enable = 'off';
            app.RetakeButton.ButtonPushedFcn = @(btn,evt) app.onRetakeTrial();

            % ============================================================
            % RIGHT-BOTTOM PANEL — Plot
            % ============================================================
            plotPanel = uipanel(app.UIFigure);
            plotPanel.Title = 'Chirp + INA219 Current';
            plotPanel.FontWeight = 'bold';
            plotPanel.Position = [560 80 760 330];

            % Top axes — expected chirp
            app.ChirpAxes = uiaxes(plotPanel);
            app.ChirpAxes.Position = [20 170 720 150];
            title(app.ChirpAxes, 'Expected Chirp');
            xlabel(app.ChirpAxes,'Time (s)');
            ylabel(app.ChirpAxes,'Norm amplitude');

            % Bottom axes — INA current
            app.INAaxes = uiaxes(plotPanel);
            app.INAaxes.Position = [20 20 720 130];
            title(app.INAaxes, 'INA219 Current');
            xlabel(app.INAaxes,'Time (s)');
            ylabel(app.INAaxes,'Current (A)');


            % ============================================================
            % Status bar
            % ============================================================
            app.StatusLabel = uilabel(app.UIFigure);
            app.StatusLabel.Position = [20 20 1200 40];
            app.StatusLabel.Text = 'Status: Not connected';
        end

        % ---------- GUI callback: connect Teensy ----------
        function onConnectTeensy(app)
            try
                if ~isempty(app.Teensy) && isvalid(app.Teensy)
                    clear app.Teensy;
                end
                port = app.PortField.Value;
                app.Teensy = serialport(port, 115200);
                pause(0.3);
                app.StatusLabel.Text = sprintf('Connected to Teensy on %s', port);
            catch ME
                app.StatusLabel.Text = sprintf('Error connecting: %s', ME.message);
            end
        end

        % ---------- GUI callback: select a condition ----------
        function onConditionSelected(app, evt)
            if isempty(evt.Indices)
                return;
            end

            row = evt.Indices(1);
            app.SelectedRowIdx = row;

            C   = app.Conditions;
            ori = C.Orientation(row);
            ind = C.Indentation(row);
            ap  = C.Aperture(row);
            tc  = C.TrialsCompleted(row);

            app.SelectedConditionLabel.Text = sprintf( ...
                'Selected condition: Ori=%g°, Ind=%g mm, Ap=%g mm^2 (Trials %d/4)', ...
                ori, ind, ap, tc);

            % Next trial index is trialsCompleted+1 (clamped to 4)
            app.CurrentTrialIdx = min(tc + 1, 4);
            app.TrialInfoLabel.Text = sprintf('Trial: %d of 4', app.CurrentTrialIdx);

            app.AcceptButton.Enable = 'off';
            app.RetakeButton.Enable = 'off';
        end

        % ---------- GUI callback: run trial ----------
        function onRunTrial(app)
            if isempty(app.SelectedRowIdx)
                app.StatusLabel.Text = 'Select a condition first.';
                return;
            end
            if isempty(app.Teensy) || ~isvalid(app.Teensy)
                app.StatusLabel.Text = 'Connect to Teensy first.';
                return;
            end

            row = app.SelectedRowIdx;
            tc  = app.Conditions.TrialsCompleted(row);
            if tc >= 4
                app.StatusLabel.Text = 'This condition already has 4 accepted trials.';
                return;
            end

            % Ensure current trial index is consistent
            app.CurrentTrialIdx = tc + 1;
            app.TrialInfoLabel.Text = sprintf('Trial: %d of 4', app.CurrentTrialIdx);

            % Get chirp parameters
            f1  = app.F1Field.Value;
            f2  = app.F2Field.Value;
            dur = app.DurField.Value;
            env = app.EnvField.Value;
            fs  = app.FsField.Value;

            cmd = sprintf('CHIRP %.4f %.4f %.4f %.4f %.4f\n', f1, f2, dur, env, fs);
            app.StatusLabel.Text = sprintf('Sending: %s', strtrim(cmd));

            try
                write(app.Teensy, cmd, "char");
            catch ME
                app.StatusLabel.Text = sprintf('Serial write error: %s', ME.message);
                return;
            end

            % ---------- listen for STARTED / INA data / DONE ----------
            % ---------- robust Teensy reader (full INA capture) ----------
            rawLines  = strings(0,1);
            N_expected = -1;
            startedSeen = false;

            tStart = tic;

            % -------- STEP 1 — wait for "STARTED" and "INA N" --------
            while toc(tStart) < (dur + 5)
                if app.Teensy.NumBytesAvailable > 0
                    line = strtrim(readline(app.Teensy));
                    sline = string(line);
                    rawLines(end+1) = sline;   %#ok<AGROW>

                    % STARTED arrives first
                    if strcmp(line,"STARTED")
                        app.StatusLabel.Text = "Chirp started...";
                        startedSeen = true;
                    end

                    % INA header: “INA N”
                    if startsWith(line,"INA ")
                        N_expected = str2double(extractAfter(line,"INA "));
                        break;   % stop scanning; go to fixed-count read
                    end
                end
                pause(0.005);
            end

            if N_expected < 1
                app.StatusLabel.Text = "ERROR: never received INA header.";
                return;
            end

            % -------- STEP 2 — read EXACTLY N_expected data lines --------
            for k = 1:N_expected
                line = strtrim(readline(app.Teensy));
                rawLines(end+1) = string(line); %#ok<AGROW>
            end

            % -------- STEP 3 — continue reading until DONE --------
            doneSeen = false;
            while toc(tStart) < (dur + 5)
                if app.Teensy.NumBytesAvailable > 0
                    line = strtrim(readline(app.Teensy));
                    rawLines(end+1) = string(line); %#ok<AGROW>
                    if strcmp(line,"DONE")
                        doneSeen = true;
                        break;
                    end
                end
                pause(0.005);
            end

            if ~doneSeen
                app.StatusLabel.Text = "WARNING: did not see DONE (but data captured).";
            else
                app.StatusLabel.Text = sprintf("Chirp finished. INA samples: %d.", N_expected);
            end


            % ---------- Parse INA219 data from rawLines ----------
            % =====================================================
            % PARSE RETURNED DATA BLOCKS (INA + HX)
            % =====================================================

            app.LastINATime    = [];
            app.LastINACurrent = [];
            app.LastHXTime     = [];
            app.LastHXValue    = [];

            if isempty(rawLines)
                app.StatusLabel.Text = 'Chirp finished. No data returned.';
            else

                % ----------- Locate headers -----------
                idxINA  = find(startsWith(rawLines,"INA "),1,'first');
                idxHX   = find(startsWith(rawLines,"HX "),1,'first');
                idxDONE = find(strcmp(rawLines,"DONE"),1,'first');

                % ----------------- INA DATA -----------------
                if ~isempty(idxINA)
                    N = str2double(extractAfter(rawLines(idxINA),"INA "));
                    if isnan(N), N = 0; end
                    first = idxINA + 1;
                    last  = min(first + N - 1, numel(rawLines));

                    t = nan(N,1);  i = nan(N,1);
                    for k = 1:(last-first+1)
                        vals = sscanf(rawLines(first+k-1), "%f,%f");
                        if numel(vals)==2
                            t(k) = vals(1);  i(k)=vals(2);
                        end
                    end
                    mask = ~isnan(t) & ~isnan(i);
                    app.LastINATime    = t(mask);
                    app.LastINACurrent = i(mask);
                end

                % ----------------- HX DATA -----------------
                if ~isempty(idxHX)
                    M = str2double(extractAfter(rawLines(idxHX),"HX "));
                    if isnan(M), M = 0; end
                    first = idxHX + 1;
                    last  = min(first + M - 1, numel(rawLines));

                    t = nan(M,1);  v = nan(M,1);
                    for k = 1:(last-first+1)
                        vals = sscanf(rawLines(first+k-1), "%f,%f");
                        if numel(vals)==2
                            t(k) = vals(1);  v(k)=vals(2);
                        end
                    end
                    mask = ~isnan(t) & ~isnan(v);
                    app.LastHXTime  = t(mask);
                    app.LastHXValue = v(mask);
                end

                app.StatusLabel.Text = sprintf("INA=%d   HX=%d samples", ...
                    numel(app.LastINATime), numel(app.LastHXTime));

            end


            % Plot expected chirp (for now; DAQ later)
            [tChirp, yNorm] = app.generateExpectedChirp(f1, f2, dur, env, fs);
            plot(app.ChirpAxes, tChirp, yNorm);
            xlabel(app.ChirpAxes,'Time (s)');
            ylabel(app.ChirpAxes,'Norm amplitude');
            title(app.ChirpAxes, sprintf('Expected chirp (Trial %d of 4)', app.CurrentTrialIdx));

            app.AcceptButton.Enable = 'on';
            app.RetakeButton.Enable = 'on';

            % ========== Plot INA219 data ==========
            if ~isempty(app.LastINATime) && ~isempty(app.LastINACurrent)
                plot(app.INAaxes, app.LastINATime, app.LastINACurrent, 'b-');
                xlabel(app.INAaxes, 'Time (s)');
                ylabel(app.INAaxes, 'Current (A)');
                title(app.INAaxes, sprintf('INA219 Current (%d samples)', numel(app.LastINATime)));
            else
                cla(app.INAaxes);
                title(app.INAaxes, 'INA219 Current (no data)');
            end

            % ========== Plot HX711 data (if available) ==========
            if ~isempty(app.LastHXTime)
                hold(app.INAaxes,'on');
                plot(app.INAaxes, app.LastHXTime, ...
                    app.LastHXValue / max(abs(app.LastHXValue)), 'r-');
                legend(app.INAaxes, {'INA Current','HX Load'}, 'Location','best');
                hold(app.INAaxes,'off');
            end


        end



        % ---------- Generate expected chirp ----------
        function [t, yNorm] = generateExpectedChirp(app, f1, f2, dur, env, fs) %#ok<INUSL>
            N  = round(dur * fs);
            t  = (0:N-1)' / fs;
            k  = (f2 - f1) / dur;
            ph = 2*pi * (f1*t + 0.5*k.*t.^2);
            s  = env * sin(ph);           % [-env, env]
            if max(abs(s)) > 0
                yNorm = s / max(abs(s));
            else
                yNorm = s;
            end
        end

        % ---------- GUI callback: accept trial ----------
        function onAcceptTrial(app)
            if isempty(app.SelectedRowIdx)
                return;
            end
            row = app.SelectedRowIdx;
            condIdx = row;  % 1:1 mapping

            tc = app.Conditions.TrialsCompleted(row);
            if app.CurrentTrialIdx ~= tc + 1 && tc < 4
                % Ensure consistent logic
                app.CurrentTrialIdx = tc + 1;
            end

            % Grab condition info
            C   = app.Conditions;
            ori = C.Orientation(row);
            ind = C.Indentation(row);
            ap  = C.Aperture(row);

            % Chirp params & expected waveform
            f1  = app.F1Field.Value;
            f2  = app.F2Field.Value;
            dur = app.DurField.Value;
            env = app.EnvField.Value;
            fs  = app.FsField.Value;
            [t, yNorm] = app.generateExpectedChirp(f1, f2, dur, env, fs);

            % Build trial struct
            trialStruct = struct();
            trialStruct.timestamp   = datestr(now);
            trialStruct.trialIndex  = app.CurrentTrialIdx;
            trialStruct.f1          = f1;
            trialStruct.f2          = f2;
            trialStruct.duration    = dur;
            trialStruct.env         = env;
            trialStruct.fs          = fs;
            trialStruct.time        = t;
            trialStruct.chirpNorm   = yNorm;

            % --- INA219 data from last run ---
            trialStruct.inaTime    = app.LastINATime;
            trialStruct.inaCurrent = app.LastINACurrent;

            % FUTURE: add measured LDV, drive, sync, load cell, etc.

            condStruct = app.SessionData.conditions(condIdx);
            if isempty(condStruct.trials)
                condStruct.trials = trialStruct;
            else
                idx = app.CurrentTrialIdx;
                if length(condStruct.trials) < idx
                    condStruct.trials(idx) = trialStruct;
                else
                    condStruct.trials(idx) = trialStruct;
                end
            end
            condStruct.trialsCompleted = app.CurrentTrialIdx;
            app.SessionData.conditions(condIdx) = condStruct;

            % Update table status (TrialsCompleted)
            app.Conditions.TrialsCompleted(row) = app.CurrentTrialIdx;
            app.TrialTable.Data = app.Conditions;

            % Save sessionData to disk
            try
                sessionData = app.SessionData; %#ok<NASGU>
                save(app.SessionFileName, 'sessionData', '-v7.3');
                app.StatusLabel.Text = sprintf( ...
                    'Accepted trial %d/4 for condition row %d. Session saved.', ...
                    app.CurrentTrialIdx, row);
            catch ME
                app.StatusLabel.Text = sprintf('Error saving session: %s', ME.message);
            end

            % If we just completed 4/4, auto-select next incomplete condition
            if app.CurrentTrialIdx >= 4
                app.AcceptButton.Enable = 'off';
                app.RetakeButton.Enable = 'off';

                % Find next condition with TrialsCompleted < 4
                nextIdx = find(app.Conditions.TrialsCompleted < 4, 1, 'first');
                if ~isempty(nextIdx)
                    app.SelectedRowIdx = nextIdx;
                    app.TrialTable.Selection = [nextIdx 1];

                    C2 = app.Conditions;
                    ori2 = C2.Orientation(nextIdx);
                    ind2 = C2.Indentation(nextIdx);
                    ap2  = C2.Aperture(nextIdx);
                    tc2  = C2.TrialsCompleted(nextIdx);

                    app.SelectedConditionLabel.Text = sprintf( ...
                        'Selected condition: Ori=%g°, Ind=%g mm, Ap=%g mm^2 (Trials %d/4)', ...
                        ori2, ind2, ap2, tc2);
                    app.CurrentTrialIdx = tc2 + 1;
                    app.TrialInfoLabel.Text = sprintf('Trial: %d of 4', app.CurrentTrialIdx);
                else
                    app.StatusLabel.Text = [app.StatusLabel.Text, ' All conditions complete.'];
                end

            else
                % Not yet 4/4: optionally auto-run next trial
                app.CurrentTrialIdx = app.CurrentTrialIdx + 1;
                app.TrialInfoLabel.Text = sprintf('Trial: %d of 4', app.CurrentTrialIdx);

                if app.AutoAdvanceCheckbox.Value
                    % Auto-run next trial
                    app.AcceptButton.Enable = 'off';
                    app.RetakeButton.Enable = 'off';
                    drawnow;
                    app.onRunTrial();
                    return;
                else
                    % Wait for user to press Run Trial again
                    app.AcceptButton.Enable = 'off';
                    app.RetakeButton.Enable = 'off';
                end
            end
        end

        % ---------- GUI callback: retake trial ----------
        function onRetakeTrial(app)
            if isempty(app.SelectedRowIdx)
                return;
            end
            % Retake just re-runs same trial number without incrementing.
            app.StatusLabel.Text = sprintf('Retaking trial %d...', app.CurrentTrialIdx);
            app.AcceptButton.Enable = 'off';
            app.RetakeButton.Enable = 'off';
            app.onRunTrial();
        end
    end

    %======================================================================
    % PUBLIC METHODS
    %======================================================================
    methods (Access = public)
        function app = FingerpadChirpApp
            createComponents(app);
            registerApp(app, app.UIFigure);
            runStartupFcn(app, @startupFcn);
        end

        function delete(app)
            % Try to clean up Teensy serial
            try
                if ~isempty(app.Teensy) && isvalid(app.Teensy)
                    clear app.Teensy;
                end
            catch
            end
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end
end
