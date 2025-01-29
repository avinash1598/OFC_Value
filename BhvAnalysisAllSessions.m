% For each session get reaction times for left and right choices

% Strucutures containing left and right RTs
bhvDataAllSession.leftRTs = [];
bhvDataAllSession.rightRTs = [];
bhvDataAllSession.leftRTsbyValue = containers.Map('KeyType', 'double', 'ValueType', 'any');
bhvDataAllSession.rightRTsbyValue = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Base directory
baseDir = '/Volumes/LaCie/bhv/';

% Get all folders starting with 'George'
folders = dir(fullfile(baseDir, 'George*'));

% Filter only directories (ignore files)
folders = folders([folders.isdir]);

% Process each folder
for i = 1:length(folders)
    % Get the full path of the folder
    folderPath = fullfile(baseDir, folders(i).name);
    
    % Display the folder being processed
    fprintf('Processing folder: %s\n', folderPath);
   
    % Get reaction times for this session
    data_ = GetReactionTimes(folderPath);
    
    if isempty(bhvDataAllSession.leftRTs)
        bhvDataAllSession.leftRTs = data_.leftRTs;
    else
        bhvDataAllSession.leftRTs = [bhvDataAllSession.leftRTs; data_.leftRTs];
    end
    
    if isempty(bhvDataAllSession.rightRTs)
        bhvDataAllSession.rightRTs = data_.rightRTs;
    else
        bhvDataAllSession.rightRTs = [bhvDataAllSession.rightRTs; data_.rightRTs];
    end
    
    % RTs by value - left choices
    uniqValues = keys(data_.leftRTsbyValue);

    for idx=1:size(data_.leftRTsbyValue.keys, 2)
        val_ = uniqValues{idx};

        if ~isKey(bhvDataAllSession.leftRTsbyValue, val_)
            bhvDataAllSession.leftRTsbyValue(val_) = [];
        end

        if isempty(bhvDataAllSession.leftRTsbyValue(val_))
            bhvDataAllSession.leftRTsbyValue(val_) = data_.leftRTsbyValue(val_);
        else
            bhvDataAllSession.leftRTsbyValue(val_) = [bhvDataAllSession.leftRTsbyValue(val_); data_.leftRTsbyValue(val_)];
        end
    end

    % RTs by value - right choices
    uniqValues = keys(data_.rightRTsbyValue);

    for idx=1:size(data_.rightRTsbyValue.keys, 2)
        val_ = uniqValues{idx};

        if ~isKey(bhvDataAllSession.rightRTsbyValue, val_)
            bhvDataAllSession.rightRTsbyValue(val_) = [];
        end

        if isempty(bhvDataAllSession.rightRTsbyValue(val_))
            bhvDataAllSession.rightRTsbyValue(val_) = data_.rightRTsbyValue(val_);
        else
            bhvDataAllSession.rightRTsbyValue(val_) = [bhvDataAllSession.rightRTsbyValue(val_); data_.rightRTsbyValue(val_)];
        end
    end
end

%% Plot distribution of left RTs and right RTs
figure;
hold on;

histogram(bhvDataAllSession.leftRTs, 'BinWidth', 10, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
histogram(bhvDataAllSession.rightRTs, 'BinWidth', 10, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

xlabel('Reaction Times (ms)');
ylabel('Frequency');
title('Histogram of Left and Right Reaction Times');
legend({'Left RTs', 'Right RTs'});
hold off;


figure;

uniqValues = keys(bhvDataAllSession.rightRTsbyValue);
for idx=1:size(bhvDataAllSession.rightRTsbyValue.keys, 2)
    val_ = uniqValues{idx};

    subplot(2, 2, idx)
    hold on

    histogram(bhvDataAllSession.leftRTsbyValue(val_), 'BinWidth', 10, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    histogram(bhvDataAllSession.rightRTsbyValue(val_), 'BinWidth', 10, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    title(sprintf('Value: %d', val_));
    xlabel('Reaction Times (ms)');
    ylabel('Frequency');
    legend({'Left RTs', 'Right RTs'});
    xlim([0, 1000])
    hold off

end

sgtitle('Histogram of Left and Right Reaction Times');

%%

% Utility function
function [returnData] = GetReactionTimes(bhvSessionDIR)
    
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

    % Retrun reaction times for left and right choices
    returnData.leftRTs = leftRTs;
    returnData.rightRTs = rightRTs;

    returnData.leftRTsbyValue = containers.Map('KeyType', 'double', 'ValueType', 'any');
    returnData.rightRTsbyValue = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % RTs by value
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

        returnData.leftRTsbyValue(val_) = leftRTs_;
        returnData.rightRTsbyValue(val_) = rightRTs_;
    end

end
