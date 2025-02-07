function [retData] = GetSignificantTimeBinsFromPVal(pValueAllNeurons, timeBins)
    warning("Make sure that p-value across time for all " + ...
        "neurons is aligned to stim onset");

    significantMask = false(size(pValueAllNeurons)); 
    
    for n = 1:size(pValueAllNeurons, 1)
        sigBins = pValueAllNeurons(n, :) < 0.01;
%         sigBins = pValueAllNeurons(n, :) < 0.05;
        count = 0;
    
        for t = 1:length(sigBins)
            if sigBins(t)
                count = count + 1;
            else
                count = 0;
            end
            
            % Mark time bins as significant if there are 4 or more consecutive significant bins
            if count >= 4
                significantMask(n, t-count+1:t) = true;
            end
        end
    end
    
    retData.significantMask = significantMask;
end






%     significantMask = pValueAllNeurons < 0.01;  
%     consecutiveCount = movsum(significantMask, 4, 2);  % Moving sum (along each row) with a window of 4
%     mask1 = consecutiveCount >= 4;
% 
%     % 500 ms window after stim onset
%     tIDX = (timeBins > 0) & (timeBins <= 500);
%     mask2 = mask1 & tIDX;
% 
%     % Index of significant neuron
%     idxSignificantNeuron = find(sum(mask2, 2) >= 1);
