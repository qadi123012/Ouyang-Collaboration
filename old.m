%% Load Standard Dataset (with concentrations)
[filename1, pathname1] = uigetfile('*.xlsx', 'Select the Standard Dataset (with concentrations)');
if isequal(filename1, 0)
    error('No file was selected for the standard dataset. Operation aborted.');
end
filepath1 = fullfile(pathname1, filename1);
data1 = readtable(filepath1, 'Sheet', 'Sheet1');

%% Load Test Dataset (without concentrations)
[filename2, pathname2] = uigetfile('*.xlsx', 'Select the Test Dataset (without concentrations)');
if isequal(filename2, 0)
    error('No file was selected for the test dataset. Operation aborted.');
end
filepath2 = fullfile(pathname2, filename2);
data2 = readtable(filepath2, 'Sheet', 'Sheet1');

%% Fix Column Names for Both Datasets
originalColNames1 = data1.Properties.VariableNames;
sanitizedColNames1 = matlab.lang.makeValidName(originalColNames1);
data1.Properties.VariableNames = sanitizedColNames1;

originalColNames2 = data2.Properties.VariableNames;
sanitizedColNames2 = matlab.lang.makeValidName(originalColNames2);
data2.Properties.VariableNames = sanitizedColNames2;

disp('Available variable names in standard dataset:');
disp(data1.Properties.VariableNames);
disp('Available variable names in test dataset:');
disp(data2.Properties.VariableNames);

%% Separate Metadata and Spectral Data for Standard Dataset
if ismember('name', data1.Properties.VariableNames)
    validRows1 = ~cellfun(@isempty, data1.name);
    metadata1 = data1(validRows1, :);
    spectralData1 = data1(~validRows1, :);
else
    error('The standard dataset does not contain a "name" column. Please check the dataset format.');
end

%% Process Spectral Data for Test Dataset
if ismember('Wavenumber', data2.Properties.VariableNames)
    spectralData2 = data2; % Assume the entire dataset is spectral data
else
    error('The test dataset does not contain a "Wavenumber" column. Please check the dataset format.');
end

%% Extract Wavenumbers and Sample Names
% For Standard Dataset
wavenumbers = spectralData1.Wavenumber;  % [nWaves x 1] vector
sampleNames1 = metadata1.name;           % e.g., {'R1', 'R2', ...}
sanitizedSampleNames1 = matlab.lang.makeValidName(sampleNames1);

% For Test Dataset
% No sample names in test dataset; process by columns directly

%% Build the Spectral Data Matrix (X) for Standard Dataset
specVarNames1 = setdiff(spectralData1.Properties.VariableNames, {'Wavenumber'});
nSamples1 = height(metadata1);   % number of samples as given in the metadata
nWaves1   = height(spectralData1); % number of wavelength points

X1 = zeros(nSamples1, nWaves1);
for i = 1:nSamples1
    pattern = ['^', sanitizedSampleNames1{i}, '[_-]'];
    idx = find(~cellfun(@isempty, regexp(specVarNames1, pattern)));
    if isempty(idx)
        error('No spectral data found for sample "%s".', sanitizedSampleNames1{i});
    end
    averagedSpectrum = mean(spectralData1{:, specVarNames1(idx)}, 2);  % [nWaves x 1]
    X1(i, :) = averagedSpectrum';  % assign as a row vector
end

%% Build the Spectral Data Matrix (X) for Test Dataset
specVarNames2 = setdiff(spectralData2.Properties.VariableNames, {'Wavenumber'});
nWaves2 = height(spectralData2); % number of wavelength points
nSamples2 = length(specVarNames2); % number of samples assumed from columns

X2 = zeros(nSamples2, nWaves2);
for i = 1:nSamples2
    X2(i, :) = spectralData2{:, specVarNames2{i}}';  % Assign each column as a row vector
end

%% Extract Concentrations from Metadata (Standard Dataset)
concentration_K2CO3 = metadata1.c_K2CO3;
concentration_KHCO3 = metadata1.c_KHCO3;
Y1 = [concentration_K2CO3, concentration_KHCO3];  % [nSamples x 2]

%% Preprocess Data (Mean Centering)
X1_mean = mean(X1, 1);
Y1_mean = mean(Y1, 1);
X1_centered = X1 - X1_mean;
Y1_centered = Y1 - Y1_mean;

%% Determine Maximum Allowable PLS Components
% Adjust the maximum allowable components based on the dataset
rnk1 = rank(X1_centered);
max_allowed_ncomp1 = min([size(X1_centered, 1) - 1, size(X1_centered, 2), rnk1]);
fprintf('Maximum allowable PLS components for these data: %d\n', max_allowed_ncomp1);

%% Cross-Validation to Determine Optimal PLS Components
mse1 = zeros(max_allowed_ncomp1, 1); % Use the dynamic max limit
k_folds1 = 5;
cv1 = cvpartition(nSamples1, 'KFold', k_folds1);
rng(0); % Set random seed for reproducibility

for ncomp = 1:max_allowed_ncomp1 % Ensure ncomp does not exceed the dynamic limit
    cv_mse = 0;
    for fold = 1:k_folds1
        train = cv1.training(fold);
        test  = cv1.test(fold);
        % Use ncomp within the allowable range
        [~, ~, ~, ~, beta_cv] = plsregress(X1_centered(train, :), Y1_centered(train, :), ncomp);
        Y_pred_cv = [ones(sum(test), 1), X1_centered(test, :)] * beta_cv;
        cv_mse = cv_mse + sum((Y_pred_cv - Y1_centered(test, :)).^2, 'all') / (k_folds1 * sum(test));
    end
    mse1(ncomp) = cv_mse;
end

[~, optimal_ncomp1] = min(mse1);
fprintf('Optimal PLS components: %d\n', optimal_ncomp1);

% Train Final PLS Model
[~, ~, ~, ~, beta1] = plsregress(X1_centered, Y1_centered, optimal_ncomp1);

%% Predict Concentrations for Test Dataset
% Preprocess Test Dataset
X2_centered = X2 - X1_mean;

% Predict Concentrations
Y2_pred_centered = [ones(nSamples2, 1), X2_centered] * beta1;
Y2_pred = Y2_pred_centered + Y1_mean;

%% Output Predicted Concentrations
predictedConcentrations = array2table(Y2_pred, ...
    'VariableNames', {'Predicted_c_K2CO3', 'Predicted_c_KHCO3'});

disp('Predicted concentrations for test dataset:');
disp(predictedConcentrations);

%% Save Predicted Concentrations to File
outputFile = fullfile(pathname2, 'Predicted_Concentrations.xlsx');
writetable(predictedConcentrations, outputFile, 'WriteRowNames', true);
fprintf('Predicted concentrations saved to: %s\n', outputFile);