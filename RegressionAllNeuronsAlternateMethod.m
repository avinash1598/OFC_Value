% Neurons aligned to onset of picture and from left OFC
% Firing rates data binned in sliding window of 100 ms and slide of 25 ms
sessionDIR = '/Volumes/LaCie/recdata/George00_rec29_04052021';
nrnData = loadNeuronsData(sessionDIR, 'LOFC', 'stimOnset');

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
        
        % Prepare the design matrix (X) and the response vector (y)
        X_ = [ones(length(trialNumbers_), 1), trialNumbers_, trialValues_, actionPerformedPerTrial_];  % Add a column of ones for the intercept
        Y_ = firingRates_;
        
        % Perform linear regression using the regress function
        [b, bint, r, rint, stats] = regress(Y_, X_);
        
        % Extract the p-value for the TrialValue predictor
        pValue = stats(3); 
       
        % TODO: record p-values for individual variable
        pValOfThisF_Neuron(t) = pValue;
        % pValTrialNumbers(nrnIDX, t) = mdl.Coefficients.pValue(2);
        % pValTrialValue(nrnIDX, t) = mdl.Coefficients.pValue(3);
        % pValActionPerformed(nrnIDX, t) = mdl.Coefficients.pValue(4);
        
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


% Number of significant neuron with time
significantNeuronsCount = sum(pValOfAllNeurons < 0.01, 1);
figure()
scatter(timeBins, significantNeuronsCount)
xlabel("time - from stim onset")
ylabel("Count of significant neurons")


%%
% Get significant neurons index in 500 ms window preceding choice onset
windowSize = 500; % ms
tIDx = (timeBins < 0) & (timeBins > -windowSize);
significantNeuronsIDx = find(sum(pValOfAllNeurons < 0.01, 2));

