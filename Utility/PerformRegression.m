function [retData] = PerformRegression(trialFiringRates, trialChosenValues, trialActionPerformed, trialNumbers)
    % trialFiringRates: num trial x time bins x num neurons
    % trialChosenValues: num trials
    % trialActionPerformed: num trials
    
    % Center these values to avoid multi-collinearity
    trialChosenValues_c = trialChosenValues - mean(trialChosenValues(~isnan(trialChosenValues)));
    trialActionPerformed_c = trialActionPerformed - mean(trialActionPerformed(~isnan(trialActionPerformed)));
    trialNumbers_c = trialNumbers - mean(trialNumbers(~isnan(trialNumbers)));
    
    % Row: Neuron, Column: time bins
    pValOverall = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2)); % For each neuron and for each trial
    pValTrialNumbers = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2));
    pValTrialValue = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2));
    pValActionPerformed = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2));
    
    % Regression analysis for all the neurons
    for nrnIDX=1:size(trialFiringRates, 3)
        singleNrnSpkData = trialFiringRates(:,:,nrnIDX);
        
        % For each time bin perform regression analysis
        for t=1:size(singleNrnSpkData, 2)
            firingRates_ = singleNrnSpkData(:,t);
            trialValues_ = trialChosenValues_c;
            actionPerformed_ = trialActionPerformed_c;
            
            % Remove rows where trialValues or firingRates is NaN
            validData_ = ~isnan(trialValues_) & ~isnan(firingRates_) & ~isnan(actionPerformed_);  % Logical index for valid (non-NaN) rows
            firingRates_ = firingRates_(validData_);  % Keep only valid firing rates
            trialValues_ = trialValues_(validData_);  % Keep only valid trial values
            actionPerformedPerTrial_ = actionPerformed_(validData_);
            trialNumbers_ = trialNumbers_c(validData_);
            
            % Create a table for fitting the model
            data = table(trialNumbers_, trialValues_, actionPerformedPerTrial_, firingRates_, ...
                'VariableNames', {'TrialNumbers', 'TrialValue', 'ActionPerformed', 'FiringRate'});
            
            % Fit the linear regression model (firing rate as a function of trial value)
            mdl = fitlm(data, 'FiringRate ~ TrialNumbers + TrialValue + ActionPerformed');
            
            % Extract overall p-value
            anovaTable = anova(mdl, 'summary');
            pValue = anovaTable.pValue(2);  % Overall p-value
           
            % TODO: record p-values for individual variable
            pValOverall(nrnIDX, t) = pValue;
            pValTrialNumbers(nrnIDX, t) = mdl.Coefficients.pValue(2);
            pValTrialValue(nrnIDX, t) = mdl.Coefficients.pValue(3);
            pValActionPerformed(nrnIDX, t) = mdl.Coefficients.pValue(4);
            
        end
    end
    
    retData.pValOverall = pValOverall;
    retData.pValTrialNumbers = pValTrialNumbers;
    retData.pValTrialValue = pValTrialValue;
    retData.pValActionPerformed = pValActionPerformed;
end