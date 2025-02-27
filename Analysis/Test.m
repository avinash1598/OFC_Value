dataBaseDIR = "C:\Users\Ranjan\Documents\MATLAB\OFC_Value\Data";
folders = dir(dataBaseDIR);
folders = folders([folders.isdir]); % Filter only directories (ignore files)

% Initialize the list
ofcSides = {'LOFC', 'ROFC'};      % Cell array of character vectors
   
for bIdx = 1:length(ofcSides)
    ofcSide = ofcSides{bIdx};       % Extract string from cell

    % TODO: create an array for factor as well
    cumDataAllSession.pval.Value.(ofcSide) = [];
    cumDataAllSession.betaCoeff.Value.(ofcSide) = [];
end

% Iterate through each session file
for i = 1:length(folders) 
    if startsWith(folders(i).name, '.')
        continue;  % Skip the current iteration if file is a '._' file
    end

    session = folders(i).name;
    sessionFilePath = fullfile(dataBaseDIR, session, "pvalue_regrerssion_all_factors.mat");

    % TODO: Check if all the time bins are same
    regressionData = load(sessionFilePath); 

    % for aIdx = 1:length(actionDirs)
    for bIdx = 1:length(ofcSides)
        ofcSide = ofcSides{bIdx};       % Extract string from cell
        
        disp(ofcSide)

        % Value information
        pval_ = regressionData.dataToSave.(ofcSide).pValTrialValue; % TODO: change
        cumDataAllSession.pval.Value.(ofcSide) = [cumDataAllSession.pval.Value.(ofcSide); pval_]; % Append rows from this session

        % Beta coefficient
        betaCoeff_ = regressionData.dataToSave.(ofcSide).betaCoeffTrialValue; % TODO: change
        cumDataAllSession.betaCoeff.Value.(ofcSide) = [cumDataAllSession.betaCoeff.Value.(ofcSide); betaCoeff_]; % Append rows from this session

    end
    % end

    % Assuming all the session have same bins - safe assumption
    cumDataAllSession.timeBins = regressionData.dataToSave.LOFC.timeBins;
end

% Get index of neuron which are significant only within 0 to 500 ms after
% stim onset

figure
hold on

% Store value dynamics in a structure first
for idx1=1:length(ofcSides) % 
    tbins = cumDataAllSession.timeBins;
    pvalData = cumDataAllSession.pval.Value.(ofcSides{idx1});
    % betaCoeffData = cumDataAllSession.betaCoeff.Value.(ofcSides{idx1});
    
    % Get significant time bins (4 consuctve time bin p-val less than 0.01)
    retData1 = GetSignificantTimeBinsFromPVal(pvalData, tbins); 
    significantMask = retData1.significantMask;

    % Get significant neuron IDx (significant timebin within 0 - 500 ms interval)
    retData2 = GetSignificantNeuronIDx(significantMask, tbins);
    significantNrnIDx = retData2.sigNrnIDx;
    
    significantNeurons.(ofcSides{idx1}) = significantNrnIDx;

    % TODO: Plot here
    valDynamics = sum(significantMask(significantNrnIDx, :), 1) / size(pvalData, 1);
    plot(tbins, valDynamics, DisplayName=ofcSides{idx1})
    xlabel("Time bins (ms)")
    ylabel("Prop of sig cells")
end

retData.significantNeuronsMask = significantNeurons;

hold off