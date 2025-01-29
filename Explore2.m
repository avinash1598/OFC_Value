% Neurons aligned to onset of picture and from left OFC
% Firing rates data binned in sliding window of 100 ms and slide of 25 ms
sessionDIR = '/Volumes/LaCie/recdata/George00_rec29_04052021';
nrnData = loadNeuronsData(sessionDIR, 'LOFC', 'choiceOnset');

bhvSessionDIR = '/Volumes/LaCie/bhv/George00_rec29_04052021';
bhvData = loadBehavioralData(bhvSessionDIR);

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

% For single neuron, plot responses for each value
uniqueVals = unique(chosenValPerTrial(~isnan(chosenValPerTrial)));
nrnIDX = 3;
singleNrnSpkData = zScoredSpkData(:,:,nrnIDX);

% figure()
% hold on
% 
% % For each uniqVals, the avg trial response for given neuron
% for i=1:length(uniqueVals)
%     val_ = uniqueVals(i);
%     fltTrialIDX_ = chosenValPerTrial == val_;
%     nrnSpkRateFromGivenTrial_ = singleNrnSpkData(fltTrialIDX_, :);
%     avgTrialData_ = nanmean(nrnSpkRateFromGivenTrial_, 1);
% 
%     plot(timeBins, avgTrialData_, LineWidth=2)
% end
% 
% xlabel("Time bin centers (unit??) - w.r.t pic onset")
% ylabel("Firing rate")
% title("Z-scored spk rate")
% legend(arrayfun(@(v) ['value ' num2str(v)], uniqueVals, 'UniformOutput', false));
% hold off

zScoredSpkData = nrnData.SPKfrnorm_units; % Number of trials, Number of time bins, Number of cells
timeBins = nrnData.t_mids;

figure()
hold on

% Plot trial average of each neuron
for idx=1:size(zScoredSpkData, 3) % Iterate through each cell
    % Get average of response of all the trials for each cell
    cellData = zScoredSpkData(:,:,idx);
    avgTrialData = nanmean(cellData, 1);

    plot(timeBins, avgTrialData, LineWidth=1.5)
end

xlabel("Time bin centers (unit??)")
ylabel("Firing rate")
title("Z-scored spk rate")
hold off



