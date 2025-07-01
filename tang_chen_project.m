% predict_concentrations_with_adjusted_references.m
% This script will:
% 1. Ask the user to select the 4.9-standard file to train a regression model.
% 2. Then ask the user to select the 4.9-predict file to predict concentrations for the samples.
% 3. Validate the predicted concentrations against reference values with a 3–9% difference tolerance.

%% Step 1: Select and load the 4.9-standard file
[standardFile, standardPath] = uigetfile('*.csv', 'Select the 4.9-standard CSV file');
if isequal(standardFile, 0)
    error('No file selected for the standard data.');
end
standardFullPath = fullfile(standardPath, standardFile);
Tstandard = readtable(standardFullPath, 'VariableNamingRule', 'preserve');

% Extract wavenumbers, absorption data, and Total concentrations
wavenumbers = Tstandard.Wavenumber; % Assuming column name is 'Wavenumber'
totalConcentration = Tstandard.Total; % Assuming column name is 'Total'
absorptionData = Tstandard{:, 6:end}; % Assuming absorption data starts from the 6th column

% Average the replicates for each sample
meanAbsorptionStandard = mean(absorptionData, 2);

% Prepare the training data
X_train = [wavenumbers, meanAbsorptionStandard]; % Features: Wavenumber and averaged absorption
Y_train = totalConcentration;                    % Target: Total concentration

%% Step 2: Train the regression model
model = fitlm(X_train, Y_train);
disp('Trained Regression Model:');
disp(model);

%% Step 3: Select and load the 4.9-predict file
[predictFile, predictPath] = uigetfile('*.csv', 'Select the 4.9-predict CSV file');
if isequal(predictFile, 0)
    error('No file selected for the prediction data.');
end
predictFullPath = fullfile(predictPath, predictFile);
Tpredict = readtable(predictFullPath, 'VariableNamingRule', 'preserve');

% Extract absorption data from the predict file
absorptionDataPredict = Tpredict{:, 2:end}; % Assuming absorption data starts from the 2nd column

% Validate number of replicate columns and samples
if mod(size(absorptionDataPredict, 2), 3) ~= 0
    error('The number of replicate columns in prediction data is not a multiple of 3. Check your input file.');
end
numSamples = size(absorptionDataPredict, 2) / 3; % Each sample has 3 replicates
fprintf('Number of replicate columns in prediction data: %d\n', size(absorptionDataPredict, 2));
fprintf('Number of samples in prediction data: %d\n', numSamples);

%% Step 4: Adjust Reference Concentrations
% Assuming the pattern from the given reference concentrations:
% Reference values (adjusted to match the number of samples in the predict file)
referenceConcentrations = repmat([0.59; 2; 3.41; 1; 2; 3; 3; 2; 2; 2; 1], ceil(numSamples / 11), 1);
referenceConcentrations = referenceConcentrations(1:numSamples); % Trim to match the exact number of samples
fprintf('Number of reference concentrations provided: %d\n', length(referenceConcentrations));

% Initialize variables for storing predictions
predictedConcentrations = zeros(numSamples, 1);
sampleNames = cell(numSamples, 1);

% Loop through each sample (group of 3 replicates) and predict concentration
for sampleIdx = 1:numSamples
    % Get replicate columns for the current sample (3 replicates per sample)
    replicateStartCol = (sampleIdx - 1) * 3 + 1;
    replicateEndCol = replicateStartCol + 2;
    replicateData = absorptionDataPredict(:, replicateStartCol:replicateEndCol);
    
    % Average the replicates
    meanAbsorptionSample = mean(replicateData, 2);
    
    % Prepare prediction data for the current sample
    X_sample = [Tpredict.Wavenumber, meanAbsorptionSample];
    
    % Predict concentration for the sample
    predictedConcentration = mean(predict(model, X_sample));
    predictedConcentrations(sampleIdx) = predictedConcentration;
    
    % Store sample name
    sampleNames{sampleIdx} = sprintf('Sample-%d', sampleIdx);
    
    % Validate the prediction against the reference value
    refValue = referenceConcentrations(sampleIdx);
    lowerBound = refValue * 0.91; % 91% of the reference value
    upperBound = refValue * 1.03; % 103% of the reference value
    if predictedConcentration < lowerBound || predictedConcentration > upperBound
        fprintf('Sample %s: Predicted concentration %.2f is outside the 3–9%% range (%.2f–%.2f)\n', ...
            sampleNames{sampleIdx}, predictedConcentration, lowerBound, upperBound);
    end
end

%% Step 5: Display and optionally save the predictions
resultsTable = table(sampleNames, predictedConcentrations, ...
    'VariableNames', {'SampleName', 'PredictedConcentration'});

disp('Predicted Concentrations for Each Sample:');
disp(resultsTable);

% Save predictions to a file
[saveFile, savePath] = uiputfile('predicted_concentrations_with_validation.csv', 'Save Predicted Concentrations');
if ischar(saveFile)
    writetable(resultsTable, fullfile(savePath, saveFile));
    fprintf('Predictions saved to %s\n', fullfile(savePath, saveFile));
end  