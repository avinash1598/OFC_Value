% Neurons aligned to onset of picture and from left OFC
% Firing rates data binned in sliding window of 100 ms and slide of 25 ms
sessionDIR = '/Volumes/LaCie/recdata/George00_rec29_04052021';
nrnData = loadNeuronsData(sessionDIR, 'LOFC', 'choiceOnset');

bhvSessionDIR = '/Volumes/LaCie/bhv/George00_rec29_04052021';
bhvData = loadBehavioralData(bhvSessionDIR);

trialNumbers = bhvData.trialinfo.TrialNumber;
leverResponse = bhvData.trialinfo.lever; % -1: Left, 1: right (for each trial)
stimVal = bhvData.trialinfo.valbin_lin;  % value information for left and right pic shown as part of stimuli

% For each trial, extract the value of chosen picture
chosenValPerTrial = NaN(size(leverResponse));

chosenValPerTrial(leverResponse == -1) = stimVal(leverResponse == -1, 1);  % For -1 (left picture), value in 1st column
chosenValPerTrial(leverResponse == 1) = stimVal(leverResponse == 1, 2);    % For 1 (right picture), value in 2nd column

% For all the leverResponses of zero, set the action performed to NaN
actionPerformedPerTrial = leverResponse;
actionPerformedPerTrial(actionPerformedPerTrial == 0) = NaN;

% Load all neurons z-scored spike data
zScoredSpkData = nrnData.SPKfrnorm_units;
timeBins = nrnData.t_mids;

% Row: Neuron, Column: time bins
pValOfAllNeurons = zeros(size(zScoredSpkData, 3), size(zScoredSpkData, 2));
pValTrialNumbers = zeros(size(zScoredSpkData, 3), size(zScoredSpkData, 2));
pValTrialValue = zeros(size(zScoredSpkData, 3), size(zScoredSpkData, 2));
pValActionPerformed = zeros(size(zScoredSpkData, 3), size(zScoredSpkData, 2));

% Regression analysis for all the neurons
for nrnIDX=1:size(zScoredSpkData, 3)
    singleNrnSpkData = zScoredSpkData(:,:,nrnIDX);
    pValOfThisF_Neuron = zeros(size(timeBins));
    
    % For each time bin perform regression analysis
    for t=1:size(singleNrnSpkData, 2)
        firingRates_ = singleNrnSpkData(:,t);
        trialValues_ = chosenValPerTrial;
        
        % Remove rows where trialValues or firingRates is NaN
        validData_ = ~isnan(trialValues_) & ~isnan(firingRates_) & ~isnan(actionPerformedPerTrial);  % Logical index for valid (non-NaN) rows
        firingRates_ = firingRates_(validData_);  % Keep only valid firing rates
        trialValues_ = trialValues_(validData_);  % Keep only valid trial values
        actionPerformedPerTrial_ = actionPerformedPerTrial(validData_);
        trialNumbers_ = trialNumbers(validData_);
        
        % Create a table for fitting the model
        data = table(trialNumbers_, trialValues_, actionPerformedPerTrial_, firingRates_, ...
            'VariableNames', {'TrialNumbers', 'TrialValue', 'ActionPerformed', 'FiringRate'});
        
        % TODO: do this analysis using regress (other method) as well to
        % see if results remain the same (not TODO anymore)
        % Fit the linear regression model (firing rate as a function of trial value)
        mdl = fitlm(data, 'FiringRate ~ TrialNumbers + TrialValue + ActionPerformed');
        
        % Extract overall p-value
        anovaTable = anova(mdl, 'summary');
        pValue = anovaTable.pValue(2);  % Overall p-value
        
        % TODO: record p-values for individual variable
        pValOfThisF_Neuron(t) = pValue;
        pValTrialNumbers(nrnIDX, t) = mdl.Coefficients.pValue(2);    % P-val for TrialNumber predictor
        pValTrialValue(nrnIDX, t) = mdl.Coefficients.pValue(3);      % P-val for TrialValue predictor
        pValActionPerformed(nrnIDX, t) = mdl.Coefficients.pValue(4); % P-val for ActionPerfomed predictor
        
    end
    
    pValOfAllNeurons(nrnIDX,:) = pValOfThisF_Neuron;
end


%%
% figure()
% imagesc(pValOfAllNeurons);            % Create a heatmap
% colorbar;                             % Add color bar to indicate scale
% colormap jet;                         % Set the color map (optional)
% xlabel('time (from choice onset)');                       % Add label for x-axis
% ylabel('Neuron');                     % Add label for y-axis
% title('Heatmap of pValOfAllNeurons'); % Title of the plot
% 
% step = 5;
% xticks(1:step:size(pValOfAllNeurons, 2));  % Set x-tick positions
% xticklabels(timeBins(1:step:size(pValOfAllNeurons, 2)));

figure();

% significantMask = pValOfAllNeurons < 0.01;  % Size: [27 x 150]
% consecutiveCount = movsum(significantMask, 4, 2);  % Moving sum with a window of 4
% significantNeuronsCount = sum(consecutiveCount >= 4, 1);

significantNeuronsCount = sum(pValOfAllNeurons < 0.01, 1);
scatter(timeBins, significantNeuronsCount);

% Axis labels and title
xlabel('Time - from choice onset');
ylabel('Count of significant neurons');
title('Count of Significant Neurons Over Time');


%%

% Proportion of significant neuron with time
figure();

hold on;

% Plot 1: All significant neurons
significantNeuronsProp = sum(pValOfAllNeurons < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'All significant');

% Plot 2: Significant neurons by Trial Number
significantNeuronsProp = sum(pValTrialNumbers < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'Trial Number');

% Plot 3: Significant neurons by Trial Value
significantNeuronsProp = sum(pValTrialValue < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'Trial Value');

% Plot 4: Significant neurons by Action Performed
significantNeuronsProp = sum(pValActionPerformed < 0.01, 1) / size(zScoredSpkData, 3);
plot(timeBins, significantNeuronsProp, 'LineWidth', 2, 'DisplayName', 'Action');

% Add legend
legend('Location', 'best');

% Axis labels and title
xlabel('Time - from choice onset');
ylabel('Proportion of significant neurons');
title('Proportion of Significant Neurons Over Time');

hold off;


%%
% Get significant neurons index in 500 ms window preceding choice onset
windowSize = 500; % ms
tIDx = (timeBins < 0) & (timeBins > -windowSize);
% tIDx = (timeBins > 0) & (timeBins < windowSize);
significantNeuronsIDx = find(sum(pValTrialValue < 0.01, 2)); % P-val only for value modulated neurons

% For each neuron plot firing rate for trials grouped by different values
idxNrn = significantNeuronsIDx(3); % IDX of actual value neuron: 3, 5, 15, 19, 23, Inverted: 6, 7
groupedTrialsByValue = GetGroupedTrials(bhvData, 'byValue', true);
uniqValues = keys(groupedTrialsByValue);

figure()
hold on

legendLabels = [];

for idx=1:size(groupedTrialsByValue.keys, 2)
    val_ = uniqValues{idx};
    trialNumbersForThisVal_ = groupedTrialsByValue(val_);
    
    % Cell data averaged over the filtered trial indexes for this value
    cellData = zScoredSpkData(trialNumbersForThisVal_,:,idxNrn);
    avgTrialData = nanmean(cellData, 1);
    
    plot(timeBins, avgTrialData, LineWidth=1.5)
    legendLabels = [legendLabels, "Value " + val_];
end

xlabel("Time bins (ms) - from choice onset")
ylabel("Firing rate")
title("Value dependent Firing rate")
legend(legendLabels)
hold off

%%
% Single neuron firing rate map across all the trials

reactionTimes = bhvData.trialinfo.rt; 
fltTrialIDx = ~isnan(reactionTimes);
fltreactionTimes = reactionTimes(fltTrialIDx);

nrnIDx = significantNeuronsIDx(23);
trialData = zScoredSpkData(fltTrialIDx,:,nrnIDx); % Average over all the trials for each neuron
trialData = squeeze(trialData);

imagesc(timeBins, 1:size(trialData, 1), trialData);  % Use timeBins for x-axis
hold on
scatter(fltreactionTimes, 1:size(fltreactionTimes, 1), 1, 'white', 'filled')

% Customize the plot
colormap('jet');          % Set the colormap (e.g., 'jet', 'parula', 'hot', etc.)
colorbar;                 % Add a colorbar
xlabel('Time (ms) from stim onset');      % Label the x-axis (modify units as needed)
ylabel('Trials');        % Label the y-axis
title('Average Trial Data Heatmap for single neuron');  % Add a title

xlim([min(timeBins), max(timeBins)]);  % Set x-axis limits based on timeBins
hold off

%%
% Firing rate of all neurons 

avgTrialData = nanmean(zScoredSpkData(:,:,significantNeuronsIDx), 1); % Average over all the trials for each neuron
avgTrialData = squeeze(avgTrialData)';

% Plot the color map
imagesc(timeBins, 1:size(avgTrialData, 1), avgTrialData);  % Use timeBins for x-axis

% Customize the plot
colormap('jet');          % Set the colormap (e.g., 'jet', 'parula', 'hot', etc.)
colorbar;                 % Add a colorbar
xlabel('Time (ms) from choice onset');      % Label the x-axis (modify units as needed)
ylabel('Neurons');        % Label the y-axis
title('Average Trial Data Heatmap');  % Add a title
