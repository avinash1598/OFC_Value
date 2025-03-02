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
    
    % Perform this for all trials (left and right)

    % for action=[-1, 1] % -1 for left and 1 for right
    % Consider trials with follwing criteria:
    % 1. only single saccade
    % 2. leverResponse is either 1 or -1
    % 3. Chosen value is not NaN i.e. in 1, 2, 3, 4
    fltValidTrialIndex = (numSaccades == 1) & (~isnan(chosenValPerTrial));
    
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
    % LOFC
    for idx = 1:numel(fields)
        dataToSave.LOFC.(fields{idx}) = regressionAnalysisDataLOFC.(fields{idx});
    end
    dataToSave.LOFC.timeBins = timeBinsLOFC;

    % ROFC
    for idx = 1:numel(fields)
        dataToSave.ROFC.(fields{idx}) = regressionAnalysisDataROFC.(fields{idx});
    end
    dataToSave.ROFC.timeBins = timeBinsROFC;
    
    sessionPath = fullfile('C:\Users\Ranjan\Documents\MATLAB\OFC_Value\Data', session);
    
    % Check if the directory exists
    if ~exist(sessionPath, 'dir')
        mkdir(sessionPath); % Create the directory if it doesn't exist
    end
    
    sessionPath = fullfile(sessionPath, "pvalue_regrerssion_all_factors.mat");
    save(sessionPath, 'dataToSave')
end

disp("done")