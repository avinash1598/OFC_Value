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
    
    % Estimated regression coefficient (beta)
    betaCoeffTrialNumbers = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2));
    betaCoeffTrialValue = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2));
    betaCoeffActionPerformed = zeros(size(trialFiringRates, 3), size(trialFiringRates, 2));
    
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
            
            % Only include predictors which have some variance in the data
            % Check variance of each predictor
            predictors = {'TrialNumbers', 'TrialValue', 'ActionPerformed'};
            data = table(trialNumbers_, trialValues_, actionPerformedPerTrial_, firingRates_, ...
                'VariableNames', {'TrialNumbers', 'TrialValue', 'ActionPerformed', 'FiringRate'});
            validPredictors = predictors(var(table2array(data(:, predictors))) > 0); % Keep only variables with non-zero variance
            
            % Construct the regression formula dynamically
            if isempty(validPredictors)
                continue; % Skip if no predictors remain
            end
            formula = ['FiringRate ~ ', strjoin(validPredictors, ' + ')];

            mdl = fitlm(data, formula);

            % % Create a table for fitting the model
            % data = table(trialNumbers_, trialValues_, actionPerformedPerTrial_, firingRates_, ...
            %     'VariableNames', {'TrialNumbers', 'TrialValue', 'ActionPerformed', 'FiringRate'});
            
            % Fit the linear regression model (firing rate as a function of trial value)
            % mdl = fitlm(data, 'FiringRate ~ TrialNumbers + TrialValue + ActionPerformed');
            
            % Extract overall p-value
            anovaTable = anova(mdl, 'summary');
            pValue = anovaTable.pValue(2);  % Overall p-value
           
            % TODO: record p-values for individual variable
            % pValOverall(nrnIDX, t) = pValue;
            % pValTrialNumbers(nrnIDX, t) = mdl.Coefficients.pValue(2);
            % pValTrialValue(nrnIDX, t) = mdl.Coefficients.pValue(3);
            % pValActionPerformed(nrnIDX, t) = mdl.Coefficients.pValue(4);
            % mdl.Coefficients.Estimate(1) => Intercept
            % mdl.Coefficients.Estimate(2) => Coefficeint for factor 1
            % mdl.Coefficients.Estimate(3) => Coefficeint for factor 2
            % mdl.Coefficients.Estimate(4) => Coefficeint for factor 3
            % Store results
            
            pValOverall(nrnIDX, t) = pValue;
            if any(strcmp(validPredictors, 'TrialNumbers'))
                pValTrialNumbers(nrnIDX, t) = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames, 'TrialNumbers'));
            end
            if any(strcmp(validPredictors, 'TrialValue'))
                pValTrialValue(nrnIDX, t) = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames, 'TrialValue'));
            end
            if any(strcmp(validPredictors, 'ActionPerformed'))
                pValActionPerformed(nrnIDX, t) = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames, 'ActionPerformed'));
            end
            
            % Estimated coefficients for regression factor
            if any(strcmp(validPredictors, 'TrialNumbers'))
                betaCoeffTrialNumbers(nrnIDX, t) = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, 'TrialNumbers'));
            end
            if any(strcmp(validPredictors, 'TrialValue'))
                betaCoeffTrialValue(nrnIDX, t) = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, 'TrialValue'));
            end
            if any(strcmp(validPredictors, 'ActionPerformed'))
                betaCoeffActionPerformed(nrnIDX, t) = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, 'ActionPerformed'));
            end
            
        end
    end

    disp(formula)
    retData.pValOverall = pValOverall;
    retData.pValTrialNumbers = pValTrialNumbers;
    retData.pValTrialValue = pValTrialValue;
    retData.pValActionPerformed = pValActionPerformed;

    retData.betaCoeffTrialNumbers = betaCoeffTrialNumbers;
    retData.betaCoeffTrialValue = betaCoeffTrialValue;
    retData.betaCoeffActionPerformed = betaCoeffActionPerformed;
end