function [retData] = GetSignificantNeuronIDx(pValMask, timeBins)
    % pValMask - Format of matrix should be same as structure returned by 
    % the function GetSignificantTimeBinsFromPVal.

    % Get index of timebins from 0 to 500 ms
    idx = (timeBins > 0) & (timeBins <= 500);
    sigNrnIDx = ( sum(pValMask(:,idx), 2) >= 1);

    retData.sigNrnIDx = sigNrnIDx;
end