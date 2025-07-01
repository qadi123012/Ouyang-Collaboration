%% Load Data
[filename, pathname] = uigetfile('*.xlsx', 'Select the Excel file');
filepath = fullfile(pathname, filename);
data = readtable(filepath, 'Sheet', 'Sheet1');

%% Fix Column Names
% MATLAB modifies column headers to be valid identifiers.
originalColNames = data.Properties.VariableNames;
sanitizedColNames = matlab.lang.makeValidName(originalColNames);
data.Properties.VariableNames = sanitizedColNames;
disp('Available variable names:');
disp(data.Properties.VariableNames);

%% Separate Metadata and Spectral Data
% Assume metadata rows have nonempty "name"; spectral-data rows have an empty "name".
validRows = ~cellfun(@isempty, data.name);
metadata = data(validRows, :);
spectralData = data(~validRows, :);

%% Extract Wavenumbers and Sample Names
% The spectralData table must contain a "Wavenumber" column.
wavenumbers = spectralData.Wavenumber;  % [nWaves x 1] vector
% The metadata table contains one row per sample.
sampleNames = metadata.name;            % e.g., {'R1', 'R2', ...}
sanitizedSampleNames = matlab.lang.makeValidName(sampleNames);

%% Build the Spectral Data Matrix (X) by Averaging Replicates
% In spectralData, aside from "Wavenumber", all columns are replicate absorption measurements.
% Their headers follow the pattern 'R<number>_1', 'R<number>_2', etc.
specVarNames = setdiff(spectralData.Properties.VariableNames, {'Wavenumber'});
nSamples = height(metadata);   % number of samples as given in the metadata
nWaves   = height(spectralData); % number of wavelength points

% Preallocate X so that each row corresponds to one sample’s averaged spectrum.
X = zeros(nSamples, nWaves);
for i = 1:nSamples
    % Build a pattern to match replicate columns for sample i.
    % For instance, if sample "R1" is given (after sanitization, "R1"),
    % then use a pattern that matches "R1_1", "R1_2", etc.
    pattern = ['^', sanitizedSampleNames{i}, '[_-]'];
    idx = find(~cellfun(@isempty, regexp(specVarNames, pattern)));
    if isempty(idx)
        error('No spectral data found for sample "%s".', sanitizedSampleNames{i});
    end
    % Average the replicate columns point‐by‐point.
    averagedSpectrum = mean(spectralData{:, specVarNames(idx)}, 2);  % [nWaves x 1]
    X(i, :) = averagedSpectrum';  % assign as a row vector
end

%% Get Concentrations (Response Variables)
% Extract concentration values from metadata.
concentration_K2CO3 = metadata.c_K2CO3;
concentration_KHCO3 = metadata.c_KHCO3;
Y = [concentration_K2CO3, concentration_KHCO3];  % [nSamples x 2]

%% Classical Least Squares (CLS) to Estimate Pure Component Spectra
% Regress X = Y * S for S (pure component spectra).
S = Y \ X;  % S has dimensions [2 x nWaves] (each row corresponds to a component)
pure_K2CO3 = S(1, :);
pure_KHCO3 = S(2, :);
X_reconstructed = Y * S;  % (Optional) Reconstructed spectra

%% Preprocess Data (Mean Centering)
X_mean = mean(X, 1);
Y_mean = mean(Y, 1);
X_centered = X - X_mean;
Y_centered = Y - Y_mean;

%% Determine Maximum Allowable PLS Components (Full Data)
% For the full dataset, a standard upper bound is:
%     min( nSamples - 1, nPredictors, rank(X_centered) )
rnk = rank(X_centered);
max_allowed_ncomp = min([size(X_centered,1)-1, size(X_centered,2), rnk]);
% Initially, use the minimum of 10 and this allowed number.
max_components = min(10, max_allowed_ncomp);
fprintf('Initial maximum allowable PLS components (full data): %d\n', max_components);

%% Adjust Maximum Components Based on CV Training Partitions
% When performing k-fold CV, some training sets may have fewer samples or lower rank.
k_folds = 5;
cv = cvpartition(nSamples, 'KFold', k_folds);
global_allowed = inf;
for fold = 1:k_folds
    train_idx = cv.training(fold);
    % For this training set, allowed latent variables is:
    allowed_fold = min(sum(train_idx)-1, rank(X_centered(train_idx, :)));
    global_allowed = min(global_allowed, allowed_fold);
end
% Adjust max_components for CV to be the smaller of the current setting and the CV-limit.
max_components = min(max_components, global_allowed);
fprintf('Adjusted maximum allowable PLS components for CV: %d\n', max_components);

%% Cross-Validation to Determine Optimal PLS Components
mse = zeros(max_components, 1);
hasStatsToolbox = license('test', 'statistics_toolbox');
rng(0); % Set random seed for reproducibility

for ncomp = 1:max_components
    cv_mse = 0;
    for fold = 1:k_folds
        train = cv.training(fold);
        test  = cv.test(fold);
        % Since we've adjusted max_components, ncomp is always within allowed limits.
        [~, ~, ~, ~, beta_cv] = plsregress(X_centered(train, :), Y_centered(train, :), ncomp);
        Y_pred_cv = [ones(sum(test),1), X_centered(test, :)] * beta_cv;
        cv_mse = cv_mse + sum((Y_pred_cv - Y_centered(test, :)).^2, 'all') / (k_folds * sum(test));
    end
    mse(ncomp) = cv_mse;
end

[~, optimal_ncomp] = min(mse);
fprintf('Optimal PLS components: %d\n', optimal_ncomp);

%% Train Final PLS Model
[~, ~, ~, ~, beta] = plsregress(X_centered, Y_centered, optimal_ncomp);
Y_pred_centered = [ones(nSamples, 1), X_centered] * beta;
Y_pred = Y_pred_centered + Y_mean;

%% Model Evaluation
r2_k2co3 = 1 - sum((Y(:,1) - Y_pred(:,1)).^2) / sum((Y(:,1) - mean(Y(:,1))).^2);
rmse_k2co3 = sqrt(mean((Y(:,1) - Y_pred(:,1)).^2));
r2_khco3 = 1 - sum((Y(:,2) - Y_pred(:,2)).^2) / sum((Y(:,2) - mean(Y(:,2))).^2);
rmse_khco3 = sqrt(mean((Y(:,2) - Y_pred(:,2)).^2));

fprintf('c_K2CO3 Performance:\nR² = %.4f, RMSE = %.4f\n', r2_k2co3, rmse_k2co3);
fprintf('c_KHCO3 Performance:\nR² = %.4f, RMSE = %.4f\n', r2_khco3, rmse_khco3);

%% Plot Pure Component Spectra (CLS)
figure;
plot(wavenumbers, pure_K2CO3, 'b-', 'LineWidth', 2);
hold on;
plot(wavenumbers, pure_KHCO3, 'r-', 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})');
ylabel('Absorbance');
title('Estimated Pure Component Spectra (CLS)');
legend('c\_K2CO3', 'c\_KHCO3', 'Location', 'best');
grid on;
hold off;

%% Plot Measured Spectra with Concentration Tooltips
% In this plot, the x-axis is wavenumbers and the y-axis is the averaged absorbance.
figure;
hold on;
colors = lines(nSamples);
for i = 1:nSamples
    % Each sample’s averaged absorption spectrum is a row of X.
    absorption = X(i, :)';  % [nWaves x 1]
    h = plot(wavenumbers, absorption, 'Color', colors(i,:), 'DisplayName', sampleNames{i});
    % Store the actual and predicted concentrations for use in data tips.
    set(h, 'UserData', {Y(i,:), Y_pred(i,:)});
end
xlabel('Wavenumber (cm^{-1})');
ylabel('Absorbance');
title('Spectra with Concentration Predictions');
legend('Location', 'bestoutside');
grid on;
hold off;

%% Configure Data Tips for the Spectral Plot
dcm = datacursormode(gcf);
set(dcm, 'UpdateFcn', @enhancedDatatip);

%% Plot Actual vs. Predicted Concentrations
figure;
% c_K2CO3 subplot
subplot(1,2,1);
scatter(Y(:,1), Y_pred(:,1), 40, 'filled', 'MarkerFaceColor', [0 0.447 0.741]);
hold on;
minVal = min([Y(:,1); Y_pred(:,1)]);
maxVal = max([Y(:,1); Y_pred(:,1)]);
plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 1.5);
xlabel('Actual c\_K2CO3 (mol/L)');
ylabel('Predicted c\_K2CO3 (mol/L)');
title(sprintf('c\_K2CO3\nR² = %.4f, RMSE = %.4f', r2_k2co3, rmse_k2co3));
grid on;
legend('Predictions', '1:1 Line', 'Location', 'southeast');
axis equal tight;

% c_KHCO3 subplot
subplot(1,2,2);
scatter(Y(:,2), Y_pred(:,2), 40, 'filled', 'MarkerFaceColor', [0.85 0.325 0.098]);
hold on;
minVal = min([Y(:,2); Y_pred(:,2)]);
maxVal = max([Y(:,2); Y_pred(:,2)]);
plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 1.5);
xlabel('Actual c\_KHCO3 (mol/L)');
ylabel('Predicted c\_KHCO3 (mol/L)');
title(sprintf('c\_KHCO3\nR² = %.4f, RMSE = %.4f', r2_khco3, rmse_khco3));
grid on;
legend('Predictions', '1:1 Line', 'Location', 'southeast');
axis equal tight;

%% Custom Data Tip Function
function txt = enhancedDatatip(~, event_obj)
    target = get(event_obj.Target, 'UserData');
    pos = get(event_obj, 'Position');
    actual = target{1};
    predicted = target{2};
    
    error1 = predicted(1) - actual(1);
    error2 = predicted(2) - actual(2);
    
    txt = { ...
        ['Wavenumber: ', num2str(pos(1))], ...
        ['Absorbance: ', num2str(pos(2))], ...
        '--- Actual ---', ...
        ['c_K2CO3: ', num2str(actual(1)), ' mol/L'], ...
        ['c_KHCO3: ', num2str(actual(2)), ' mol/L'], ...
        '--- Predicted ---', ...
        ['c_K2CO3: ', num2str(predicted(1)), ' mol/L (Δ: ', num2str(error1), ')'], ...
        ['c_KHCO3: ', num2str(predicted(2)), ' mol/L (Δ: ', num2str(error2), ')'] ...
    };
end