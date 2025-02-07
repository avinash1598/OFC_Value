
clear all

bhvBaseDir = '/Volumes/LaCie/bhv'; % Base directory
recdataBaseDir = '/Volumes/LaCie/recdata'; % Base directory
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
    
    sessionPath = fullfile('/Users/avinashranjan/Documents/MATLAB/ElstonLab/OFC_Value/Data', session);
    
    % Check if the directory exists
    if ~exist(sessionPath, 'dir')
        mkdir(sessionPath); % Create the directory if it doesn't exist
    end
    
    sessionPath = fullfile(sessionPath, "pvalue_LOFC_vs_ROFC_ActionGrouped.mat");
    save(sessionPath, 'dataToSave')
end

disp("done")

%%
dataBaseDIR = "/Users/avinashranjan/Documents/MATLAB/ElstonLab/OFC_Value/Data";
folders = dir(dataBaseDIR);
folders = folders([folders.isdir]); % Filter only directories (ignore files)

% Initialize the list
actionDirs = {'left', 'right'};   % Cell array of character vectors
ofcSides = {'LOFC', 'ROFC'};      % Cell array of character vectors
brainSides = {'IPSI', 'CONTRA'};
    
for aIdx = 1:length(actionDirs)
    for bIdx = 1:length(ofcSides)
        actionDir = actionDirs{aIdx};   % Extract string from cell
        ofcSide = ofcSides{bIdx};       % Extract string from cell
        
        % Determine IPSI or CONTRA based on the mapping rule
        if (strcmp(actionDir, 'left') && strcmp(ofcSide, 'LOFC')) || ...
           (strcmp(actionDir, 'right') && strcmp(ofcSide, 'ROFC'))
            brainSide = 'IPSI';
        else
            brainSide = 'CONTRA';
        end

        % TODO: create an array for factor as well
        cumDataAllSession.pval.Value.(ofcSide).(brainSide) = [];
        cumDataAllSession.betaCoeff.Value.(ofcSide).(brainSide) = [];
    end
end

for i = 1:length(folders) % Iterate through each session file
    if startsWith(folders(i).name, '.')
        continue;  % Skip the current iteration if file is a '._' file
    end

    session = folders(i).name;
    sessionFilePath = fullfile(dataBaseDIR, session, "pvalue_LOFC_vs_ROFC_ActionGrouped.mat");

    % TODO: Check if all the time bins are same
    regressionData = load(sessionFilePath); 

    for aIdx = 1:length(actionDirs)
        for bIdx = 1:length(ofcSides)
            actionDir = actionDirs{aIdx};   % Extract string from cell
            ofcSide = ofcSides{bIdx};       % Extract string from cell
            
            % Determine IPSI or CONTRA based on the mapping rule
            if (strcmp(actionDir, 'left') && strcmp(ofcSide, 'LOFC')) || ...
               (strcmp(actionDir, 'right') && strcmp(ofcSide, 'ROFC'))
                brainSide = 'IPSI';
            else
                brainSide = 'CONTRA';
            end

            disp([actionDir, ' ', ofcSide, ' ', brainSide])

            % Value information
            pval_ = regressionData.dataToSave.(actionDir).(ofcSide).pValTrialValue; % TODO: change
            cumDataAllSession.pval.Value.(ofcSide).(brainSide) = [cumDataAllSession.pval.Value.(ofcSide).(brainSide); pval_]; % Append rows from this session

            % Beta coefficient
            betaCoeff_ = regressionData.dataToSave.(actionDir).(ofcSide).betaCoeffTrialValue; % TODO: change
            cumDataAllSession.betaCoeff.Value.(ofcSide).(brainSide) = [cumDataAllSession.betaCoeff.Value.(ofcSide).(brainSide); betaCoeff_]; % Append rows from this session

        end
    end

    % Assuming all the session have same bins - safe assumption
    cumDataAllSession.timeBins = regressionData.dataToSave.left.LOFC.timeBins;
end

% Get index of neuron which are significant only within 0 to 500 ms after
% stim onset

valueDynamicsDataStruct = [];

% Store value dynamics in a structure first
for idx1=1:length(ofcSides) % 
    for idx2=1:length(brainSides)
        
        tbins = cumDataAllSession.timeBins;
        pvalData = cumDataAllSession.pval.Value.(ofcSides{idx1}).(brainSides{idx2});
        betaCoeffData = cumDataAllSession.betaCoeff.Value.(ofcSides{idx1}).(brainSides{idx2});
        
        % Get significant time bins
        retData1 = GetSignificantTimeBinsFromPVal(pvalData, tbins); 
        significantMask = retData1.significantMask;
        
        % Get significant neuron IDx (significant timebin within 0 - 500 ms interval)
        retData2 = GetSignificantNeuronIDx(significantMask, tbins);
        significantNrnIDx = retData2.sigNrnIDx;
        
        % Get first significant timebin post stim onset
        % And get this only for significant neuron
        fristSigBinPostStim = GetFirstSigPostStim(significantMask(significantNrnIDx, :), tbins);
        
        % value dynamics
        % Get significant proportion at each time bin for significant neuron only
        significantNeuronsProp = sum(significantMask(significantNrnIDx, :), 1) / size(pvalData, 1);
        % significantNeuronsProp = sum(significantMask, 1) / size(data, 1);

        % beta dynamics
        betaDynamics = mean(abs(betaCoeffData(significantNrnIDx, :)), 1);

        % Mean beta value (0-500 ms) of all neurons (not just significant)
        tIDx = (tbins > 0) & (tbins <= 500);
        meanBeta = mean(betaCoeffData(:,tIDx), 2);

        valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).(brainSides{idx2}) = significantNeuronsProp;
        valueDynamicsDataStruct.betaDynamics.Value.(ofcSides{idx1}).(brainSides{idx2}) = betaDynamics;
        valueDynamicsDataStruct.fristSigBinPostStim.Value.(ofcSides{idx1}).(brainSides{idx2}) = fristSigBinPostStim;
        valueDynamicsDataStruct.meanBeta.Value.(ofcSides{idx1}).(brainSides{idx2}) = meanBeta;
        valueDynamicsDataStruct.significantNrnMask.Value.(ofcSides{idx1}).(brainSides{idx2}) = significantNrnIDx;
        
    end

end

% Perform and store Chi-squared test
for idx1=1:length(ofcSides)
    % Chi-squared test for significantNeuronsProp
    chiTestData.prop1 = valueDynamicsDataStruct.significantNeuronsProp.Value.LOFC.IPSI;
    chiTestData.prop2 = valueDynamicsDataStruct.significantNeuronsProp.Value.LOFC.CONTRA;
    chiTestData.totCount1 = size(cumDataAllSession.pval.Value.(ofcSides{idx1}).IPSI, 1);
    chiTestData.totCount2 = size(cumDataAllSession.pval.Value.(ofcSides{idx1}).CONTRA, 1);
    
    retData = PerformChiSquaredTest(chiTestData);
    
    valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).ChiTestResult = retData.pvals;
end


figure
hold on

% plotType = ['pval', 'betaCoeff'];
% factor = ['Value', 'Action'];
% ofcSides = {'LOFC', 'ROFC'};
% brainSides = {'CONTRA', 'IPSI'};

for idx1=1:length(ofcSides) % 
    for idx2=1:length(brainSides)
        
        tbins = cumDataAllSession.timeBins;
        significantNeuronsProp = valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).(brainSides{idx2});
        betaDynamics = valueDynamicsDataStruct.betaDynamics.Value.(ofcSides{idx1}).(brainSides{idx2});
        fristSigBinPostStim = valueDynamicsDataStruct.fristSigBinPostStim.Value.(ofcSides{idx1}).(brainSides{idx2});
        
        % Plot value dynamics
        subplot(3, 2, idx1)
        hold on
        plot(tbins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', brainSides{idx2});

        title(ofcSides{idx1})
        legend('Location', 'best');
        xlabel('Time - from stim onset (ms)');
        ylabel('Proportion of significant value neurons');
        
        % Plot beta dynamics
        
        subplot(3, 2, 2 + idx1)
        hold on
        plot(tbins, betaDynamics, 'LineWidth', 2, 'DisplayName', brainSides{idx2});
        
        title(ofcSides{idx1})
        legend('Location', 'best');
        xlabel('Time - from stim onset (ms)');
        ylabel('Mean value regression coeff');

        % Histogram of first significant time bin post stim onset
        subplot(3, 2, 4 + idx1)
        hold on
        h = histogram(fristSigBinPostStim, 'NumBins', 17);  % Create the histogram
        h.DisplayName = brainSides{idx2};  % Assign the DisplayName property

        title(ofcSides{idx1})
        legend('Location', 'best');
        xlabel('Frist significant timebin post stim-onset (ms)');
        ylabel('Count');
    end

    % Overlay significance markers (e.g., small bars at the bottom)
    pvals = valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).ChiTestResult;
    significantBins = pvals < 0.05;

    subplot(3, 2, idx1)
    bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

    % ttest result on first sig bin post stim-onset
    d1 = valueDynamicsDataStruct.fristSigBinPostStim.Value.(ofcSides{idx1}).IPSI;
    d2 = valueDynamicsDataStruct.fristSigBinPostStim.Value.(ofcSides{idx1}).CONTRA;
    [h, p, ci, stats] =  ttest2(d1, d2);
    disp(p)
end

hold off


%%
% P-values IPSI vs CONTRA
figure

for idx1=1:length(ofcSides) % 
        
    tbins = cumDataAllSession.timeBins;
    meanBetaIPSI = valueDynamicsDataStruct.meanBeta.Value.(ofcSides{idx1}).IPSI;
    meanBetaCONTRA = valueDynamicsDataStruct.meanBeta.Value.(ofcSides{idx1}).CONTRA;

    significantNrnMaskIPSI = valueDynamicsDataStruct.significantNrnMask.Value.(ofcSides{idx1}).IPSI;
    significantNrnMaskCONTRA = valueDynamicsDataStruct.significantNrnMask.Value.(ofcSides{idx1}).CONTRA;
    
    nrnSigInBothIPSI_CONTRA = significantNrnMaskIPSI & significantNrnMaskCONTRA;

    % Define colors
    colors = repmat([0 0 1], length(nrnSigInBothIPSI_CONTRA), 1);
    colors(nrnSigInBothIPSI_CONTRA, :) = repmat([1 0 0], sum(nrnSigInBothIPSI_CONTRA), 1); % Set significant points to red
    
    subplot(1, 2, idx1)
    hold on
    scatter(meanBetaIPSI, meanBetaCONTRA, 50, colors, 'filled', 'HandleVisibility', 'off'); % Scatter plot

    % Add invisible scatter plots for legend
    scatter(NaN, NaN, 50, 'b', 'filled', 'DisplayName', 'Non-significant'); % Blue legend
    scatter(NaN, NaN, 50, 'r', 'filled', 'DisplayName', 'Significant in both IPSI/CONTRA'); % Red legend

    % Set labels and title
    title(ofcSides{idx1});
    xlabel('Beta value IPSI');
    ylabel('Beta value CONTRA');
    
    % Show legend
    legend('Location', 'best');
    hold off

end



