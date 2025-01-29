
% Do wrt stimulus onset
sessionDIR = '/Volumes/LaCie/recdata/George00_rec29_04052021';
nrnData = loadNeuronsData(sessionDIR, 'LOFC', 'stimOnset');

bhvSessionDIR = '/Volumes/LaCie/bhv/George00_rec29_04052021';
bhvData = loadBehavioralData(bhvSessionDIR);

% Prepare the behavioral data needed for regression analysis
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
    
% Load all neurons z-scored spike data
zScoredSpkData = nrnData.SPKfrnorm_units;
timeBins = nrnData.t_mids;
flt_zScoredSpkData = zScoredSpkData(fltValidTrialIndex, :, :); % Filtered data

% Perform regresssion analysis
regressionAnalysisData = PerformRegression( ...
    flt_zScoredSpkData, flt_chosenValPerTrial, flt_leverResponse, flt_trialNumbers);


%% Plots
% Proportion of significant neuron with time
figure();

hold on;

% Plot 1: All significant neurons
significantNeuronsProp = sum(regressionAnalysisData.pValOverall < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'All significant');

% Plot 2: Significant neurons by Trial Number
significantNeuronsProp = sum(regressionAnalysisData.pValTrialNumbers < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'Trial Number');

% Plot 3: Significant neurons by Trial Value
significantNeuronsProp = sum(regressionAnalysisData.pValTrialValue < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'Trial Value');

% Plot 4: Significant neurons by Action Performed
significantNeuronsProp = sum(regressionAnalysisData.pValActionPerformed < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'Action');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from stim onset');
ylabel('Proportion of significant neurons');
title('Proportion of Significant Neurons Over Time');

hold off;
