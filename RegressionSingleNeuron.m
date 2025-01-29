% Neurons aligned to onset of picture and from left OFC
% Firing rates data binned in sliding window of 100 ms and slide of 25 ms
% sessionDIR = '/Volumes/LaCie/recdata/George00_rec29_04052021';
% nrnData = loadNeuronsData(sessionDIR, 'LOFC', 'stimulusOnset');
% 
% bhvSessionDIR = '/Volumes/LaCie/bhv/George00_rec29_04052021';
% bhvData = loadBehavioralData(bhvSessionDIR);

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

% Regression analysis for single neuron
nrnIDX = 22;
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

    %     % Convert trial values to categorical variable
    %     trialValues_ = categorical(trialValues_);
    %     
    %     % Convert action performed to categorical variable
    %     actionPerformedPerTrial_ = categorical(actionPerformedPerTrial_);
    
    % Create a table for fitting the model
    data = table(trialNumbers_, trialValues_, actionPerformedPerTrial_, firingRates_, ...
        'VariableNames', {'TrialNumbers', 'TrialValue', 'ActionPerformed', 'FiringRate'});
    
    % Fit the linear regression model (firing rate as a function of trial value)
    mdl = fitlm(data, 'FiringRate ~ TrialNumbers + TrialValue + ActionPerformed');
    
    % Extract the p-value for the TrialValue predictor
    pValue = mdl.Coefficients.pValue(2);  % 2nd row corresponds to TrialValue
   
    pValOfThisF_Neuron(t) = pValue;
end

%%
figure()
plot(timeBins, pValOfThisF_Neuron)
xlabel("time")
ylabel("P-val")
ylim([0, 0.05])