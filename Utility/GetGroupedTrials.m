function [returnData] = GetGroupedTrials(bhvData, varargin)
    % Define default values
    byValue = false;  % Default for 'byValue'
    
    % Parse optional arguments
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'byValue'
                byValue = varargin{i+1};
            otherwise
                error('Unknown parameter name: %s', varargin{i});
        end
    end
    
    trialNumbers = bhvData.trialinfo.TrialNumber;
    leverResponse = bhvData.trialinfo.lever; % -1: Left, 1: right (for each trial)
    stimVal = bhvData.trialinfo.valbin_lin;  % value information for left and right pic shown as part of stimuli
    
    if byValue == true
        % For each trial, extract the value of chosen picture
        chosenValPerTrial = NaN(size(leverResponse));
    
        chosenValPerTrial(leverResponse == -1) = stimVal(leverResponse == -1, 1);  % For -1 (left picture), value in 1st column
        chosenValPerTrial(leverResponse == 1) = stimVal(leverResponse == 1, 2);    % For 1 (right picture), value in 2nd column
        
        % Get unique possible values for the stimuli and get trial indexes for
        % each of those unique values
        uniqueVals = unique(chosenValPerTrial(~isnan(chosenValPerTrial)));
        
        % Create a dictionary (containers.Map)
        returnData = containers.Map('KeyType', 'double', 'ValueType', 'any');
        
        for i=1:length(uniqueVals)
            val = uniqueVals(i);
            returnData(val) = trialNumbers(chosenValPerTrial == val);
        end
    else
        disp("WTF")
    end
    
    % No else case handling for now
end