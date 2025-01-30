function [retData] = GetFirstSigPostStim(sigMaskFromPValues, tbins)
    % Find the time bin for each neuron when it first becomes significant
    % after the stim onset
    dataMatrix = sigMaskFromPValues & (tbins > 0); % Consider only time bins after stim onset
    first_significant_bin_idx = arrayfun(@(r) find(dataMatrix(r, :), 1, 'first'), 1:size(dataMatrix, 1), 'UniformOutput', false);

    % Replace empty arrays with NaN to handle cases where no '1' is found
    % firstSigBin(cellfun(@isempty, num2cell(firstSigBin))) = NaN;

    % Convert cell array to a numeric array
    first_significant_bin_idx_array = cell2mat(first_significant_bin_idx);
    
    % Replace zeros with NaN
    first_significant_bin_idx_array(first_significant_bin_idx_array == 0) = NaN;

    % retData = tbins(first_significant_bin_idx);
    retData = tbins(first_significant_bin_idx_array);
end