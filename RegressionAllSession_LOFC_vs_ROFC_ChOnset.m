
clear allShohini Ghose

bhvBaseDir = '/Volumes/LaCie/bhv/'; % Base directory
recdataBaseDir = '/Volumes/LaCie/recdata/'; % Base directory
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
    
    % Consider trials with follwing criteria:
    % 1. only single saccade
    % 2. leverResponse is either 1 or -1
    % 3. Chosen value is not NaN i.e. in 1, 2, 3, 4
    fltValidTrialIndex = (numSaccades == 1) & ( (leverResponse == 1) | (leverResponse == -1) ) & (~isnan(chosenValPerTrial));

    flt_chosenValPerTrial = chosenValPerTrial(fltValidTrialIndex); % Filtered data
    flt_leverResponse = leverResponse(fltValidTrialIndex); % Filtered data
    flt_trialNumbers = trialNumbers(fltValidTrialIndex); % Filtered data
    
    % Load data from LOFC (from stim onset) and perform regression analysis
    nrnDataLOFC = loadNeuronsData(recdataFolderPath, 'LOFC', 'choiceOnset');
    
    zScoredSpkDataLOFC = nrnDataLOFC.SPKfrnorm_units;
    timeBinsLOFC = nrnDataLOFC.t_mids;
    flt_zScoredSpkDataLOFC = zScoredSpkDataLOFC(fltValidTrialIndex, :, :); % Filtered data
    
    regressionAnalysisDataLOFC = PerformRegression( ...
        flt_zScoredSpkDataLOFC, flt_chosenValPerTrial, flt_leverResponse, flt_trialNumbers);
    
    % Load data from ROFC (from stim onset) and perform regression analysis
    nrnDataROFC = loadNeuronsData(recdataFolderPath, 'ROFC', 'choiceOnset');
    
    zScoredSpkDataROFC = nrnDataROFC.SPKfrnorm_units;
    timeBinsROFC = nrnDataROFC.t_mids;
    flt_zScoredSpkDataROFC = zScoredSpkDataROFC(fltValidTrialIndex, :, :); % Filtered data

    regressionAnalysisDataROFC = PerformRegression( ...
        flt_zScoredSpkDataROFC, flt_chosenValPerTrial, flt_leverResponse, flt_trialNumbers);
    
    % Save data to file
    dataToSave.LOFC.pValOverall = regressionAnalysisDataLOFC.pValOverall;
    dataToSave.LOFC.pValTrialNumbers = regressionAnalysisDataLOFC.pValTrialNumbers;
    dataToSave.LOFC.pValTrialValue = regressionAnalysisDataLOFC.pValTrialValue;
    dataToSave.LOFC.pValActionPerformed = regressionAnalysisDataLOFC.pValActionPerformed;
    dataToSave.LOFC.timeBins = timeBinsLOFC;

    dataToSave.ROFC.pValOverall = regressionAnalysisDataROFC.pValOverall;
    dataToSave.ROFC.pValTrialNumbers = regressionAnalysisDataROFC.pValTrialNumbers;
    dataToSave.ROFC.pValTrialValue = regressionAnalysisDataROFC.pValTrialValue;
    dataToSave.ROFC.pValActionPerformed = regressionAnalysisDataROFC.pValActionPerformed;
    dataToSave.ROFC.timeBins = timeBinsROFC;
    
    sessionPath = fullfile('/Users/avinashranjan/Documents/MATLAB/ElstonLab/OFC_Value/Data/', session);
    
    % Check if the directory exists
    if ~exist(sessionPath, 'dir')
        mkdir(sessionPath); % Create the directory if it doesn't exist
    end
    
    sessionPath = fullfile(sessionPath, "pvalue_LOFC_vs_ROFC_ChOnset.mat");
    save(sessionPath, 'dataToSave')
end

%%
pval_ValueLOFC = [];
pval_ValueROFC = [];
timeBinsLOFC = [];
timeBinsROFC = [];

dataBaseDIR = '/Users/avinashranjan/Documents/MATLAB/ElstonLab/OFC_Value/Data';
folders = dir(dataBaseDIR);
folders = folders([folders.isdir]); % Filter only directories (ignore files)

for i = 1:length(folders)
    if startsWith(folders(i).name, '.')
        continue;  % Skip the current iteration if file is a '._' file
    end

    session = folders(i).name;
    sessionFilePath = fullfile(dataBaseDIR, session, "pvalue_LOFC_vs_ROFC_ChOnset.mat");

    % Check if all the time bins are same
    regressionData = load(sessionFilePath); 
    
    pval_ValueLOFC_ = regressionData.dataToSave.LOFC.pValActionPerformed;
    pval_ValueLOFC = [pval_ValueLOFC; pval_ValueLOFC_]; % Append rows from this session
    timeBinsLOFC = [timeBinsLOFC; regressionData.dataToSave.LOFC.timeBins];

    pval_ValueROFC_ = regressionData.dataToSave.ROFC.pValActionPerformed;
    pval_ValueROFC = [pval_ValueROFC; pval_ValueROFC_]; % Append rows from this session
    timeBinsROFC = [timeBinsROFC; regressionData.dataToSave.ROFC.timeBins];
    
end

figure;

% Plot 1: All significant neurons
significantNeuronsProp = sum(pval_ValueLOFC < 0.01, 1) / size(pval_ValueLOFC, 1);
plot(mean(timeBinsLOFC, 1), significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'LOFC');

hold on;

% Plot 1: All significant neurons
significantNeuronsProp = sum(pval_ValueROFC < 0.01, 1) / size(pval_ValueROFC, 1);
plot(mean(timeBinsROFC, 1), significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'ROFC');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from choice onset (ms)');
ylabel('Proportion of significant neurons');
title('Proportion of Significant Action Neurons');

hold off;
