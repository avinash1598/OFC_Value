
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
    
    sessionPath = fullfile('C:\Users\Ranjan\Documents\MATLAB\OFC_Value\Data', session);
    
    % Check if the directory exists
    if ~exist(sessionPath, 'dir')
        mkdir(sessionPath); % Create the directory if it doesn't exist
    end
    
    sessionPath = fullfile(sessionPath, "pvalue_LOFC_vs_ROFC_ActionGrouped.mat");
    save(sessionPath, 'dataToSave')
end

disp("done")

%%
dataBaseDIR = "C:\Users\Ranjan\Documents\MATLAB\OFC_Value\Data";
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

% Get all the value neurons considering both left and right action
valueNrnMask_LR_Action = GetAllValueCells();

valueDynamicsDataStruct = [];

% Store value dynamics in a structure first
for idx1=1:length(ofcSides) % 
    % Get significant value neurons found considering both left and right
    % action
    significantNrnIDx = valueNrnMask_LR_Action.significantNeuronsMask.(ofcSides{idx1});
    
    for idx2=1:length(brainSides)
        
        tbins = cumDataAllSession.timeBins;
        pvalData = cumDataAllSession.pval.Value.(ofcSides{idx1}).(brainSides{idx2});
        betaCoeffData = cumDataAllSession.betaCoeff.Value.(ofcSides{idx1}).(brainSides{idx2});
        
        % Get significant time bins
        retData1 = GetSignificantTimeBinsFromPVal(pvalData, tbins); 
        significantMask = retData1.significantMask;
        
        % Get significant neuron IDx (significant timebin within 0 - 500 ms interval)
        % Instead of this maybe get all the value neurons as mask.
        % retData2 = GetSignificantNeuronIDx(significantMask, tbins);
        % significantNrnIDx = retData2.sigNrnIDx; % Don't do this

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

% Aggregate from LOFC and ROFC for IPSI and CONTRA
for idx2=1:length(brainSides)
    significantNrnIDx_L = valueNrnMask_LR_Action.significantNeuronsMask.LOFC;
    significantNrnIDx_R = valueNrnMask_LR_Action.significantNeuronsMask.ROFC;

    tbins = cumDataAllSession.timeBins;
    pvalData_L = cumDataAllSession.pval.Value.LOFC.(brainSides{idx2});
    betaCoeffData_L = cumDataAllSession.betaCoeff.Value.LOFC.(brainSides{idx2});

    pvalData_R = cumDataAllSession.pval.Value.ROFC.(brainSides{idx2});
    betaCoeffData_R = cumDataAllSession.betaCoeff.Value.ROFC.(brainSides{idx2});

    % Concatinate data from LOFC and ROFC
    pvalData = [pvalData_L; pvalData_R];
    betaCoeffData = [betaCoeffData_L; betaCoeffData_R];
    significantNrnIDx = [significantNrnIDx_L; significantNrnIDx_R];
    
    % Get significant time bins
    retData1 = GetSignificantTimeBinsFromPVal(pvalData, tbins); 
    significantMask = retData1.significantMask;
    
    % Get first significant timebin post stim onset
    % And get this only for significant neuron
    fristSigBinPostStim = GetFirstSigPostStim(significantMask(significantNrnIDx, :), tbins);
    
    % value dynamics
    % Get significant proportion at each time bin for significant neuron only
    significantNeuronsProp = sum(significantMask(significantNrnIDx, :), 1) / size(pvalData, 1);
    % significantNeuronsProp = sum(significantMask, 1) / size(data, 1);
    
    % beta dynamics
    betaDynamics = mean(abs(betaCoeffData(significantNrnIDx, :)), 1);

    valueDynamicsDataStruct.significantNeuronsProp.Value.(brainSides{idx2}) = significantNeuronsProp;
    valueDynamicsDataStruct.betaDynamics.Value.(brainSides{idx2}) = betaDynamics;
    valueDynamicsDataStruct.fristSigBinPostStim.Value.(brainSides{idx2}) = fristSigBinPostStim;
    valueDynamicsDataStruct.significantNrnMask.Value.(brainSides{idx2}) = significantNrnIDx;
    valueDynamicsDataStruct.betaCoeffData.Value.(brainSides{idx2}) = abs(betaCoeffData(significantNrnIDx, :)); % significant neuron only
end


% Perform and store Chi-squared test on proportion neurons
for idx1=1:length(ofcSides)
    % Chi-squared test for significantNeuronsProp
    chiTestData.prop1 = valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).IPSI;
    chiTestData.prop2 = valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).CONTRA;
    chiTestData.totCount1 = size(cumDataAllSession.pval.Value.(ofcSides{idx1}).IPSI, 1);
    chiTestData.totCount2 = size(cumDataAllSession.pval.Value.(ofcSides{idx1}).CONTRA, 1);
    
    retData = PerformChiSquaredTest(chiTestData);
    valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).ChiTestResult = retData.pvals;
end

chiTestData.prop1 = valueDynamicsDataStruct.significantNeuronsProp.Value.IPSI;
chiTestData.prop2 = valueDynamicsDataStruct.significantNeuronsProp.Value.CONTRA;
chiTestData.totCount1 = length(valueDynamicsDataStruct.significantNrnMask.Value.IPSI);
chiTestData.totCount2 = length(valueDynamicsDataStruct.significantNrnMask.Value.CONTRA);

retData = PerformChiSquaredTest(chiTestData);
valueDynamicsDataStruct.significantNeuronsProp.Value.ChiTestResult = retData.pvals;

%
% Perform Cluster-Based Permutation Test (Better for Multiple Comparisons with data that has correlation over time)
% Compute empirical t-statistics
% Cluster based permutation test better for correlated time points data
FR_cond1 = valueDynamicsDataStruct.betaCoeffData.Value.IPSI;
FR_cond2 = valueDynamicsDataStruct.betaCoeffData.Value.CONTRA;

alpha = 0.05; % Cluster-forming threshold
n_perms = 1000;
T = size(FR_cond1, 2); % Number of time points

% Step 1: Compute actual t-values
[~, ~, ~, stats] = ttest(FR_cond1, FR_cond2);
t_real = stats.tstat;

% Step 2: Find clusters in real data
sig_clusters = bwconncomp(abs(t_real) > tinv(1-alpha, size(FR_cond1,1)-1)); 

% Step 3: Permutation test
max_cluster_sums = zeros(n_perms,1);

for i = 1:n_perms
    % Randomly flip conditions
    flip_sign = (rand(size(FR_cond1,1),1) > 0.5);
    FR_perm1 = FR_cond1;
    FR_perm2 = FR_cond2;
    FR_perm1(flip_sign,:) = FR_cond2(flip_sign,:);
    FR_perm2(flip_sign,:) = FR_cond1(flip_sign,:);
    
    % Compute permuted t-values
    [~, ~, ~, perm_stats] = ttest(FR_perm1, FR_perm2);
    t_perm = perm_stats.tstat;
    
    % Find max cluster sum in permuted data
    perm_clusters = bwconncomp(abs(t_perm) > tinv(1-alpha, size(FR_cond1,1)-1));
    cluster_sums = cellfun(@(x) sum(abs(t_perm(x))), perm_clusters.PixelIdxList);
    max_cluster_sums(i) = max([cluster_sums, 0]); % Include 0 to avoid empty case
end

% Step 4: Find significance threshold
cluster_thresh = prctile(max_cluster_sums, 95); % 95th percentile

% Step 5: Identify significant clusters
real_cluster_sums = cellfun(@(x) sum(abs(t_real(x))), sig_clusters.PixelIdxList);
sig_clusters_idx = find(real_cluster_sums > cluster_thresh);

% Extract significant time points
sig_time_points = cell2mat(sig_clusters.PixelIdxList(sig_clusters_idx));



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
        subplot(3, 3, idx1)
        hold on
        plot(tbins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', brainSides{idx2});

        title(ofcSides{idx1})
        legend('Location', 'best');
        xlabel('Time - from stim onset (ms)');
        ylabel('Proportion of significant value neurons');
        xlim([-1000, 1000])
        
        % Plot beta dynamics
        
        subplot(3, 3, 3 + idx1)
        hold on
        plot(tbins, betaDynamics, 'LineWidth', 2, 'DisplayName', brainSides{idx2});
        
        title(ofcSides{idx1})
        legend('Location', 'best');
        xlabel('Time - from stim onset (ms)');
        ylabel('Mean value regression coeff');
        xlim([-1000, 1000])

        % Histogram of first significant time bin post stim onset
        subplot(3, 3, 6 + idx1)
        hold on
        edges = linspace(0, 1000, 18); % 17 bins between 0 and 1000
        h = histogram(fristSigBinPostStim, 'BinEdges', edges);
        % h = histogram(fristSigBinPostStim, 'NumBins', 17);  % Create the histogram
        h.DisplayName = brainSides{idx2};  % Assign the DisplayName property

        title(ofcSides{idx1})
        legend('Location', 'best');
        xlabel('Frist significant timebin post stim-onset (ms)');
        ylabel('Count');
    end

    % Overlay significance markers (e.g., small bars at the bottom)
    pvals = valueDynamicsDataStruct.significantNeuronsProp.Value.(ofcSides{idx1}).ChiTestResult;
    significantBins = pvals < 0.05;

    subplot(3, 3, idx1)
    bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

    % ttest result on first sig bin post stim-onset
    d1 = valueDynamicsDataStruct.fristSigBinPostStim.Value.(ofcSides{idx1}).IPSI;
    d2 = valueDynamicsDataStruct.fristSigBinPostStim.Value.(ofcSides{idx1}).CONTRA;
    [h, p, ci, stats] =  ttest2(d1, d2);
    disp(p)
end


for idx2=1:length(brainSides)
    
    tbins = cumDataAllSession.timeBins;
    significantNeuronsProp = valueDynamicsDataStruct.significantNeuronsProp.Value.(brainSides{idx2});
    betaDynamics = valueDynamicsDataStruct.betaDynamics.Value.(brainSides{idx2});
    fristSigBinPostStim = valueDynamicsDataStruct.fristSigBinPostStim.Value.(brainSides{idx2});
    
    % Plot value dynamics
    subplot(3, 3, 3)
    hold on
    plot(tbins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', brainSides{idx2});

    legend('Location', 'best');
    xlabel('Time - from stim onset (ms)');
    ylabel('Proportion of significant value neurons');
    xlim([-1000, 1000])
    
    % Plot beta dynamics
    subplot(3, 3, 3 + 3)
    hold on
    plot(tbins, betaDynamics, 'LineWidth', 2, 'DisplayName', brainSides{idx2});
    
    legend('Location', 'best');
    xlabel('Time - from stim onset (ms)');
    ylabel('Mean value regression coeff');
    xlim([-1000, 1000])

    % Histogram of first significant time bin post stim onset
    subplot(3, 3, 6 + 3)
    hold on
    edges = linspace(0, 1000, 18); % 17 bins between 0 and 1000
    tempArr = fristSigBinPostStim(fristSigBinPostStim < 600);
    h = histogram(tempArr, 'BinEdges', edges);
    % h = histogram(fristSigBinPostStim, 'NumBins', 17);  % Create the histogram
    h.DisplayName = brainSides{idx2};  % Assign the DisplayName property
    
    legend('Location', 'best');
    xlabel('Frist significant timebin post stim-onset (ms)');
    ylabel('Count');
    % xlim([])
end

% Overlay significance markers (e.g., small bars at the bottom)
% Significance test for cell proportion markers
pvals = valueDynamicsDataStruct.significantNeuronsProp.Value.ChiTestResult;
significantBins = pvals < 0.05;

subplot(3, 3, 3)
bar(tbins, significantBins * 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
  
% Beta regression coefficient t-test
subplot(3, 3, 6)
title("Cluster-based permutation test")
bar(tbins(sig_time_points), zeros(size(sig_time_points)) + 0.01, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

% ttest result on first sig bin post stim-onset
d1 = valueDynamicsDataStruct.fristSigBinPostStim.Value.IPSI;
d2 = valueDynamicsDataStruct.fristSigBinPostStim.Value.CONTRA;
[h, p, ci, stats] =  ttest2(d1, d2);
disp(p)


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



