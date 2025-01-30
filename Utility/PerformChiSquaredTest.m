function [retData] = PerformChiSquaredTest(inputData)
    prop1 = inputData.prop1;
    prop2 = inputData.prop2;
    totCount1 = inputData.totCount1;
    totCount2 = inputData.totCount2;
    
    % if (totCount1 ~= totCount2) || (size(prop1) ~= size(prop2))
    %     disp("WTF total cell count is different for both?")
    %     return
    % end
    
    if size(prop1) ~= size(prop2)
        disp("WTF no of time bins are different for both?")
        return
    end
    
    pVals = nan(size(prop1)); % Preallocate p-values

    for t = 1:length(prop1)
        % Get proportions
        prop1_ = prop1(t);
        prop2_ = prop2(t);
    
        % Convert proportions into a contingency table using a common total (e.g., 100)
        active1_ = round(prop1_ * totCount1);
        inactive1_ = totCount1 - active1_;
    
        active2_ = round(prop2_ * totCount2);
        inactive2_ = totCount2 - active2_;
    
        % Construct 2x2 contingency table
        contingencyTable = [active1_, inactive1_; 
                            active2_, inactive2_];
        
        % Perform Chi-squared test
        % [~, pVals(t), ~] = chi2test(contingencyTable); % Correct test function
        observed = contingencyTable;
        expected = sum(observed, 2) * sum(observed, 1) / sum(observed, 'all');
        
        % Compute chi-square statistic
        chi2_stat = sum((observed - expected).^2 ./ expected, 'all');
        
        % Get p-value using chi-square CDF
        df = (size(contingencyTable, 1) - 1) * (size(contingencyTable, 2) - 1);
        pVals(t) = 1 - chi2cdf(chi2_stat, df);
    end
    
    retData.pvals = pVals;
end