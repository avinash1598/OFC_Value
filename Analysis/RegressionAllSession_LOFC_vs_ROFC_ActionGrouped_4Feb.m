
clear all

bhvBaseDir = 'C:\Users\Ranjan\Desktop\Data\bhv'; % Base directory
recdataBaseDir = 'C:\Users\Ranjan\Desktop\Data\recdata'; % Base directory
folders = dir(fullfile(bhvBaseDir, 'George*')); % Get all folders starting with 'George'
folders = folders([folders.isdir]); % Filter only directories (ignore files)

% Process each folder
for i = 1:length(folders)
    % Get the full path of the folder
    session = folders(i).name;
    bhvFolderPath = fullfile(bhvBaseDir, folders(i).name);
    recdataFolderPath = fullfile(recdataBaseDir, folders(i).name);
    
    disp(session)

    % Prepare the behavioral data needed for regression analysis
    bhvData = loadBehavioralData(bhvFolderPath); % Load behavioral data

    numSaccades = bhvData.trialinfo.saccN;
    trialNumbers = bhvData.trialinfo.TrialNumber;
    leverResponse = bhvData.trialinfo.lever; % -1: Left, 1: right (for each trial)
    stimVal = bhvData.trialinfo.valbin_lin;  % value information for left and right pic shown as part of stimuli
    
    chosenValPerTrial = NaN(size(leverResponse)); % For each trial, extract the value of chosen picture
    chosenValPerTrial(leverResponse == -1) = stimVal(leverResponse == -1, 1);  % For -1 (left picture), value in 1st column
    chosenValPerTrial(leverResponse == 1) = stimVal(leverResponse == 1, 2);    % For 1 (right picture), value in 2nd column
    
    % TODO: Perform this separately for left and right action
    for action=[-1, 1] % -1 for left and 1 for right
        % Consider trials with follwing criteria:
        % 1. only single saccade
        % 2. leverResponse is either 1 or -1
        % 3. Chosen value is not NaN i.e. in 1, 2, 3, 4
        fltValidTrialIndex = (numSaccades == 1) & (leverResponse == action) & (~isnan(chosenValPerTrial));
        
        flt_chosenValPerTrial = chosenValPerTrial(fltValidTrialIndex); % Filtered data
        flt_leverResponse = leverResponse(fltValidTrialIndex); % Filtered data
        flt_trialNumbers = trialNumbers(fltValidTrialIndex); % Filtered data
        
        % Load data from LOFC (from stim onset) and perform regression analysis
        nrnDataLOFC = loadNeuronsData(recdataFolderPath, 'LOFC', 'stimOnset');
        
        zScoredSpkDataLOFC = nrnDataLOFC.SPKfrnorm_units;
        timeBinsLOFC = nrnDataLOFC.t_mids;
        flt_zScoredSpkDataLOFC = zScoredSpkDataLOFC(fltValidTrialIndex, :, :); % Filtered data
        
        regressionAnalysisDataLOFC = PerformRegression( ...
            flt_zScoredSpkDataLOFC, flt_chosenValPerTrial, flt_leverResponse, flt_trialNumbers);
        
        % Load data from ROFC (from stim onset) and perform regression analysis
        nrnDataROFC = loadNeuronsData(recdataFolderPath, 'ROFC', 'stimOnset');
        
        zScoredSpkDataROFC = nrnDataROFC.SPKfrnorm_units;
        timeBinsROFC = nrnDataROFC.t_mids;
        flt_zScoredSpkDataROFC = zScoredSpkDataROFC(fltValidTrialIndex, :, :); % Filtered data
        
        regressionAnalysisDataROFC = PerformRegression( ...
            flt_zScoredSpkDataROFC, flt_chosenValPerTrial, flt_leverResponse, flt_trialNumbers);
        
        fields = {'pValOverall', 'pValTrialNumbers', 'pValTrialValue', 'pValActionPerformed', 'betaCoeffTrialNumbers', 'betaCoeffTrialValue', 'betaCoeffActionPerformed'};

        % Save data to file
        if action == -1 % Left

            % LOFC
            for idx = 1:numel(fields)
                dataToSave.left.LOFC.(fields{idx}) = regressionAnalysisDataLOFC.(fields{idx});
            end
            dataToSave.left.LOFC.timeBins = timeBinsLOFC;
        
            % ROFC
            for idx = 1:numel(fields)
                dataToSave.left.ROFC.(fields{idx}) = regressionAnalysisDataROFC.(fields{idx});
            end
            dataToSave.left.ROFC.timeBins = timeBinsROFC;
            
        elseif action == 1 % Right
            % LOFC
            for idx = 1:numel(fields)
                dataToSave.right.LOFC.(fields{idx}) = regressionAnalysisDataLOFC.(fields{idx});
            end
            dataToSave.right.LOFC.timeBins = timeBinsLOFC;
        
            % ROFC
            for idx = 1:numel(fields)
                dataToSave.right.ROFC.(fields{idx}) = regressionAnalysisDataROFC.(fields{idx});
            end
            dataToSave.right.ROFC.timeBins = timeBinsROFC;
        end
        
    end
    
    sessionPath = fullfile('C:/Users/Ranjan/Documents/MATLAB/OFC_Value/Data', session);
    
    % Check if the directory exists
    if ~exist(sessionPath, 'dir')
        mkdir(sessionPath); % Create the directory if it doesn't exist
    end
    
    sessionPath = fullfile(sessionPath, "pvalue_LOFC_vs_ROFC_ActionGrouped.mat");
    save(sessionPath, 'dataToSave')
end

print("done")

%%
pval_Value.left.LOFC = [];
pval_Value.left.ROFC = [];
pval_Value.right.LOFC = [];
pval_Value.right.ROFC = [];

timeBins = [];

dataBaseDIR = "C:/Users/Ranjan/Documents/MATLAB/OFC_Value/Data";
folders = dir(dataBaseDIR);
folders = folders([folders.isdir]); % Filter only directories (ignore files)

for i = 1:length(folders)
    if startsWith(folders(i).name, '.')
        continue;  % Skip the current iteration if file is a '._' file
    end

    session = folders(i).name;
    sessionFilePath = fullfile(dataBaseDIR, session, "pvalue_LOFC_vs_ROFC_ActionGrouped.mat");

    % Check if all the time bins are same
    regressionData = load(sessionFilePath); 
    
    % Left
    pval_ = regressionData.dataToSave.left.LOFC.pValTrialValue; % TODO: change
    pval_Value.left.LOFC = [pval_Value.left.LOFC; pval_]; % Append rows from this session

    pval_ = regressionData.dataToSave.left.ROFC.pValTrialValue; % TODO: change
    pval_Value.left.ROFC = [pval_Value.left.ROFC; pval_]; % Append rows from this session

    % Right
    pval_ = regressionData.dataToSave.right.LOFC.pValTrialValue; % TODO: change
    pval_Value.right.LOFC = [pval_Value.right.LOFC; pval_]; % Append rows from this session

    pval_ = regressionData.dataToSave.right.ROFC.pValTrialValue; % TODO: change
    pval_Value.right.ROFC = [pval_Value.right.ROFC; pval_]; % Append rows from this session

    % Assuming all regions have same timebins
    timeBins = [timeBins; regressionData.dataToSave.left.LOFC.timeBins];
    
end

% Classify a neurons as value neuron if p-value 
% is significant for 4 consecutive windows 500 ms after stim onset

figure;

subplot(2, 2, 1)

% Plot 1: All significant neurons
% If p-value is significant for any neuron for 4 consecutive timebin within
% 500 ms of stim onset then the neuron is significant
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.left.LOFC, tbins);
significantMaskLOFC = retData.significantMask;

significantNeuronsProp1 = sum(significantMaskLOFC, 1) / size(pval_Value.left.LOFC, 1);
plot(tbins, significantNeuronsProp1, 'LineWidth', 2, 'DisplayName', 'LOFC');

hold on;

% Plot 1: All significant neurons
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.left.ROFC, tbins);
significantMaskROFC = retData.significantMask;

significantNeuronsProp2 = sum(significantMaskROFC, 1) / size(pval_Value.left.ROFC, 1);
plot(tbins, significantNeuronsProp2, 'LineWidth', 2, 'DisplayName', 'ROFC');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from stim onset (ms)');
ylabel('Proportion of significant neurons');
title('Left action');

% Perform Chi-squared test
chiTestData.prop1 = significantNeuronsProp1;
chiTestData.prop2 = significantNeuronsProp2;
chiTestData.totCount1 = size(pval_Value.left.LOFC, 1);
chiTestData.totCount2 = size(pval_Value.left.ROFC, 1);

retData = PerformChiSquaredTest(chiTestData);

% Overlay significance markers (e.g., small bars at the bottom)
significantBins = retData.pvals < 0.05;
bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

hold off;


subplot(2, 2, 2)

% Plot 1: All significant neurons
% If p-value is significant for any neuron for 4 consecutive timebin within
% 500 ms of stim onset then the neuron is significant
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.right.LOFC, tbins);
significantMaskLOFC = retData.significantMask;

significantNeuronsProp1 = sum(significantMaskLOFC, 1) / size(pval_Value.right.LOFC, 1);
plot(tbins, significantNeuronsProp1, 'LineWidth', 2, 'DisplayName', 'LOFC');

hold on;

% Plot 1: All significant neurons
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.right.ROFC, tbins);
significantMaskROFC = retData.significantMask;

significantNeuronsProp2 = sum(significantMaskROFC, 1) / size(pval_Value.right.ROFC, 1);
plot(tbins, significantNeuronsProp2, 'LineWidth', 2, 'DisplayName', 'ROFC');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from stim onset (ms)');
ylabel('Proportion of significant neurons');
title('Right action');

% Perform Chi-squared test
chiTestData.prop1 = significantNeuronsProp1;
chiTestData.prop2 = significantNeuronsProp2;
chiTestData.totCount1 = size(pval_Value.right.LOFC, 1);
chiTestData.totCount2 = size(pval_Value.right.ROFC, 1);

retData = PerformChiSquaredTest(chiTestData);

% Overlay significance markers (e.g., small bars at the bottom)
significantBins = retData.pvals < 0.05;
bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

hold off;


% Another way of plotting the same thing
subplot(2, 2, 3)

% Plot 1: All significant neurons
% If p-value is significant for any neuron for 4 consecutive timebin within
% 500 ms of stim onset then the neuron is significant
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.left.LOFC, tbins);
significantMask = retData.significantMask;

significantNeuronsProp1 = sum(significantMask, 1) / size(pval_Value.left.LOFC, 1);
plot(tbins, significantNeuronsProp1, 'LineWidth', 2, 'DisplayName', 'Left');

hold on;

% Plot 1: All significant neurons
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.right.LOFC, tbins);
significantMask = retData.significantMask;

significantNeuronsProp2 = sum(significantMask, 1) / size(pval_Value.right.LOFC, 1);
plot(tbins, significantNeuronsProp2, 'LineWidth', 2, 'DisplayName', 'Right');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from stim onset (ms)');
ylabel('Proportion of significant neurons');
title('LOFC');

% Perform Chi-squared test
chiTestData.prop1 = significantNeuronsProp1;
chiTestData.prop2 = significantNeuronsProp2;
chiTestData.totCount1 = size(pval_Value.left.LOFC, 1);
chiTestData.totCount2 = size(pval_Value.right.LOFC, 1);

retData = PerformChiSquaredTest(chiTestData);

% Overlay significance markers (e.g., small bars at the bottom)
significantBins = retData.pvals < 0.05;
bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

hold off;


subplot(2, 2, 4)

% Plot 1: All significant neurons
% If p-value is significant for any neuron for 4 consecutive timebin within
% 500 ms of stim onset then the neuron is significant
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.left.ROFC, tbins);
significantMask = retData.significantMask;

significantNeuronsProp1 = sum(significantMask, 1) / size(pval_Value.left.ROFC, 1);
plot(tbins, significantNeuronsProp1, 'LineWidth', 2, 'DisplayName', 'Left');

hold on;

% Plot 1: All significant neurons
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.right.ROFC, tbins);
significantMask = retData.significantMask;

significantNeuronsProp2 = sum(significantMask, 1) / size(pval_Value.right.ROFC, 1);
plot(tbins, significantNeuronsProp2, 'LineWidth', 2, 'DisplayName', 'Right');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from stim onset (ms)');
ylabel('Proportion of significant neurons');
title('ROFC');

% Perform Chi-squared test
chiTestData.prop1 = significantNeuronsProp1;
chiTestData.prop2 = significantNeuronsProp2;
chiTestData.totCount1 = size(pval_Value.left.ROFC, 1);
chiTestData.totCount2 = size(pval_Value.right.ROFC, 1);

retData = PerformChiSquaredTest(chiTestData);

% Overlay significance markers (e.g., small bars at the bottom)
significantBins = retData.pvals < 0.05;
bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

hold off;

%%
% Ipsilateral and contralateral (Isn't this same as left and right action)
figure;

% Plot 1: All significant neurons
% If p-value is significant for any neuron for 4 consecutive timebin within
% 500 ms of stim onset then the neuron is significant
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.left.LOFC, tbins);
significantMask = retData.significantMask;
retData2 = GetSignificantTimeBinsFromPVal(pval_Value.right.ROFC, tbins);
significantMask2 = retData2.significantMask;

significantNeuronsProp1 = (sum(significantMask, 1) + sum(significantMask2, 1)) / (size(pval_Value.left.LOFC, 1) + size(pval_Value.right.ROFC, 1));
plot(tbins, significantNeuronsProp1, 'LineWidth', 2, 'DisplayName', 'Ipsi');

hold on;

% Plot 1: All significant neurons
tbins = mean(timeBins, 1); 
retData = GetSignificantTimeBinsFromPVal(pval_Value.right.LOFC, tbins);
significantMask = retData.significantMask;
retData2 = GetSignificantTimeBinsFromPVal(pval_Value.left.ROFC, tbins);
significantMask2 = retData2.significantMask;

significantNeuronsProp2 = (sum(significantMask, 1) + sum(significantMask2, 1)) / (size(pval_Value.left.ROFC, 1) + size(pval_Value.right.LOFC, 1));
plot(tbins, significantNeuronsProp2, 'LineWidth', 2, 'DisplayName', 'Contra');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from stim onset (ms)');
ylabel('Proportion of significant neurons');
title('Ipsi/Contra value dynamics');

% Perform Chi-squared test
chiTestData.prop1 = significantNeuronsProp1;
chiTestData.prop2 = significantNeuronsProp2;
chiTestData.totCount1 = (size(pval_Value.left.LOFC, 1) + size(pval_Value.right.ROFC, 1));
chiTestData.totCount2 = (size(pval_Value.left.ROFC, 1) + size(pval_Value.right.LOFC, 1));

retData = PerformChiSquaredTest(chiTestData);

% Overlay significance markers (e.g., small bars at the bottom)
significantBins = retData.pvals < 0.05;
bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

hold off;