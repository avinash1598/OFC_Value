bhvSessionDIR = '/Volumes/LaCie/bhv/George00_rec29_04052021';
bhvData = loadBehavioralData(bhvSessionDIR);

trialNumbers = bhvData.trialinfo.TrialNumber;
numSaccades = bhvData.trialinfo.saccN;
leverResponse = bhvData.trialinfo.lever; % -1: Left, 1: right (for each trial)
stimVal = bhvData.trialinfo.valbin_lin;  % value information for left and right pic shown as part of stimuli
reactionTimes = bhvData.trialinfo.rt;    % Reaction times of each trial

% For each trial, extract the value of chosen picture
chosenValPerTrial = NaN(size(leverResponse));

chosenValPerTrial(leverResponse == -1) = stimVal(leverResponse == -1, 1);  % For -1 (left picture), value in 1st column
chosenValPerTrial(leverResponse == 1) = stimVal(leverResponse == 1, 2);    % For 1 (right picture), value in 2nd column

% For all the leverResponses of zero, set the action performed to NaN
actionPerformedPerTrial = leverResponse;
actionPerformedPerTrial(actionPerformedPerTrial == 0) = NaN;

% Set chosen value and action performed to NaN for trials for which number 
% of saccades is greater than 1
nanIDx = (numSaccades > 1) | isnan(numSaccades);
chosenValPerTrial(nanIDx) = NaN;
actionPerformedPerTrial(nanIDx) = NaN;
reactionTimes(nanIDx) = NaN;

% Left choice RTs and Right choice RTs
leftIDx = actionPerformedPerTrial == -1;
leftRTs = reactionTimes(leftIDx);

rightIDx = actionPerformedPerTrial == 1;
rightRTs = reactionTimes(rightIDx);

% Plot distribution of left RTs and right RTs
figure;
hold on;

histogram(leftRTs, 'BinWidth', 10, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
histogram(rightRTs, 'BinWidth', 10, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

xlabel('Reaction Times (ms)');
ylabel('Frequency');
title('Histogram of Left and Right Reaction Times');
legend({'Left RTs', 'Right RTs'});
hold off;

% Plot distribution for each value
figure;

uniqueVals = unique(chosenValPerTrial(~isnan(chosenValPerTrial)));

for i=1:length(uniqueVals)
    val_ = uniqueVals(i);
    fltIDx_ = chosenValPerTrial == val_;

    fltActions_ = actionPerformedPerTrial(fltIDx_);
    fltRTs_ = reactionTimes(fltIDx_);

    leftIDx_ = fltActions_ == -1;
    leftRTs_ = fltRTs_(leftIDx_);
    
    rightIDx_ = fltActions_ == 1;
    rightRTs_ = fltRTs_(rightIDx_);

    % Plot for each value
    subplot(2, 2, i)
    hold on

    histogram(leftRTs_, 'BinWidth', 10, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    histogram(rightRTs_, 'BinWidth', 10, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    title(sprintf('Value: %d', val_));
    xlabel('Reaction Times (ms)');
    ylabel('Frequency');
    legend({'Left RTs', 'Right RTs'});
    xlim([0, 1000])
    hold off
    
end

sgtitle('Histogram of Left and Right Reaction Times');
