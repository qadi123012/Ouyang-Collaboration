%% Load Data
[filename, pathname] = uigetfile('*.xlsx', 'Select the Excel file');
filepath = fullfile(pathname, filename);

% Load Metadata Sheet
metadata = readtable(filepath, 'Sheet', 'Metadata');
metadata = metadata(:, ~cellfun(@isempty, metadata.Properties.VariableNames));

% Load Spectra Sheet
spectralData = readtable(filepath, 'Sheet', 'Spectra');
spectralData = spectralData(:, ~cellfun(@isempty, spectralData.Properties.VariableNames));

%% Fix Column Names
% Process metadata
metaVars = metadata.Properties.VariableNames;
sanitizedMetaVars = matlab.lang.makeValidName(metaVars);
metadata.Properties.VariableNames = sanitizedMetaVars;

% Process spectral data
specVars = spectralData.Properties.VariableNames;
sanitizedSpecVars = matlab.lang.makeValidName(specVars);
spectralData.Properties.VariableNames = sanitizedSpecVars;

%% Extract Metadata
% Get concentrations
k2co3_col = find(contains(sanitizedMetaVars, 'c_K2CO3'));
khco3_col = find(contains(sanitizedMetaVars, 'c_KHCO3'));

concentration_K2CO3 = metadata.(sanitizedMetaVars{k2co3_col});
concentration_KHCO3 = metadata.(sanitizedMetaVars{khco3_col});

% Extract sample names
sampleNames = metadata.name;

%% Extract Spectral Data (X)
% Get wavenumbers
wavenumbers = spectralData.Wavenumber;

% Match sample names between metadata and spectral data
validSamples = ismember(sanitizedSpecVars, sanitizedMetaVars);
spectralColumns = find(validSamples(2:end)) + 1; % Skip Wavenumber column

% Verify sample matching
if numel(spectralColumns) ~= numel(sampleNames)
    error('Sample name mismatch between metadata and spectral data');
end

X = spectralData{:, spectralColumns}'; % Transpose to [Samples x Wavenumbers]

%% Prepare Y Matrix
Y = [concentration_K2CO3, concentration_KHCO3];

%% Classical Least Squares (CLS) to Estimate Pure Component Spectra
S = Y \ X;  % Pure component spectra (rows: K2CO3, KHCO3)
pure_K2CO3 = S(1, :);
pure_KHCO3 = S(2, :);

%% Preprocess Data (Mean Centering)
X_mean = mean(X, 1);
Y_mean = mean(Y, 1);
X_centered = X - X_mean;
Y_centered = Y - Y_mean;

%% Cross-Validation to Determine Optimal PLS Components
max_components = 10;
k_folds = 5;
mse = zeros(max_components, 1);

hasStatsToolbox = license('test', 'statistics_toolbox');
rng(0); % Set random seed for reproducibility

for ncomp = 1:max_components
    cv_mse = 0;
    
    if hasStatsToolbox
        cv = cvpartition(size(X_centered, 1), 'KFold', k_folds);
        for fold = 1:k_folds
            train = cv.training(fold);
            test = cv.test(fold);
            
            [~, ~, ~, ~, beta_cv] = plsregress(X_centered(train, :), Y_centered(train, :), ncomp);
            
            Y_pred_cv = [ones(sum(test),1), X_centered(test, :)] * beta_cv;
            cv_mse = cv_mse + sum((Y_pred_cv - Y_centered(test, :)).^2, 'all') / (k_folds * sum(test));
        end
    else
        n_samples = size(X_centered, 1);
        indices = randperm(n_samples);
        fold_size = floor(n_samples / k_folds);
        
        for fold = 1:k_folds
            test_indices = indices((fold-1)*fold_size+1 : min(fold*fold_size, n_samples));
            test = false(n_samples, 1); test(test_indices) = true;
            train = ~test;
            
            [~, ~, ~, ~, beta_cv] = plsregress(X_centered(train, :), Y_centered(train, :), ncomp);
            
            Y_pred_cv = [ones(sum(test),1), X_centered(test, :)] * beta_cv;
            cv_mse = cv_mse + sum((Y_pred_cv - Y_centered(test, :)).^2, 'all') / (k_folds * sum(test));
        end
    end
    
    mse(ncomp) = cv_mse;
end

[~, optimal_ncomp] = min(mse);
fprintf('Optimal PLS components: %d\n', optimal_ncomp);

%% Train Final PLS Model
[~, ~, ~, ~, beta] = plsregress(X_centered, Y_centered, optimal_ncomp);
Y_pred_centered = [ones(size(X_centered,1),1), X_centered] * beta;
Y_pred = Y_pred_centered + Y_mean;

%% Model Evaluation
r2_k2co3 = 1 - sum((Y(:,1) - Y_pred(:,1)).^2) / sum((Y(:,1) - mean(Y(:,1))).^2);
rmse_k2co3 = sqrt(mean((Y(:,1) - Y_pred(:,1)).^2));

r2_khco3 = 1 - sum((Y(:,2) - Y_pred(:,2)).^2) / sum((Y(:,2) - mean(Y(:,2))).^2);
rmse_khco3 = sqrt(mean((Y(:,2) - Y_pred(:,2)).^2));

fprintf('K₂CO₃ Performance:\nR² = %.4f\nRMSE = %.4f\n', r2_k2co3, rmse_k2co3);
fprintf('KHCO₃ Performance:\nR² = %.4f\nRMSE = %.4f\n', r2_khco3, rmse_khco3);

%% Plot Pure Component Spectra
figure;
plot(wavenumbers, pure_K2CO3, 'b-', 'LineWidth', 2);
hold on;
plot(wavenumbers, pure_KHCO3, 'r-', 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})');
ylabel('Absorbance');
title('Estimated Pure Component Spectra (CLS)');
legend('K₂CO₃', 'KHCO₃', 'Location', 'best');
grid on;
hold off;

%% Plot Measured Spectra with Concentration Tooltips
figure;
hold on;
colors = lines(numel(sampleNames));

for i = 1:numel(sampleNames)
    absorption = X(i, :);
    h = plot(wavenumbers, absorption, ...
             'Color', colors(i,:), ...
             'DisplayName', sampleNames{i});
    
    % Store actual and predicted values for tooltips
    set(h, 'UserData', {[concentration_K2CO3(i), concentration_KHCO3(i)], Y_pred(i,:)});
end

xlabel('Wavenumber (cm^{-1})');
ylabel('Absorption');
title('Spectra with Concentration Predictions');
legend('Location', 'bestoutside');
grid on;

%% Configure Data Tips for Spectral Plot
dcm = datacursormode(gcf);
set(dcm, 'UpdateFcn', @enhancedDatatip);

%% Plot Actual vs. Predicted Concentrations
figure;

% K₂CO₃ subplot
subplot(1,2,1);
scatter(Y(:,1), Y_pred(:,1), 40, 'filled', 'MarkerFaceColor', [0 0.447 0.741]);
hold on;
minVal = min([Y(:,1); Y_pred(:,1)]);
maxVal = max([Y(:,1); Y_pred(:,1)]);
plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 1.5);
xlabel('Actual K₂CO₃ (mol/L)');
ylabel('Predicted K₂CO₃ (mol/L)');
title(sprintf('K₂CO₃\nR² = %.4f, RMSE = %.4f', r2_k2co3, rmse_k2co3));
grid on;
legend('Predictions', '1:1 Line', 'Location', 'southeast');
axis equal tight

% KHCO₃ subplot
subplot(1,2,2);
scatter(Y(:,2), Y_pred(:,2), 40, 'filled', 'MarkerFaceColor', [0.85 0.325 0.098]);
hold on;
minVal = min([Y(:,2); Y_pred(:,2)]);
maxVal = max([Y(:,2); Y_pred(:,2)]);
plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 1.5);
xlabel('Actual KHCO₃ (mol/L)');
ylabel('Predicted KHCO₃ (mol/L)');
title(sprintf('KHCO₃\nR² = %.4f, RMSE = %.4f', r2_khco3, rmse_khco3));
grid on;
legend('Predictions', '1:1 Line', 'Location', 'southeast');
axis equal tight

%% Custom Data Tip Function
function txt = enhancedDatatip(~, event_obj)
    target = get(event_obj.Target, 'UserData');
    pos = get(event_obj, 'Position');
    actual = target{1};
    predicted = target{2};
    
    error_k2co3 = predicted(1) - actual(1);
    error_khco3 = predicted(2) - actual(2);
    
    txt = { 
        ['Wavenumber: ', num2str(pos(1))],
        ['Absorption: ', num2str(pos(2))],
        '--- Actual ---',
        ['K₂CO₃: ', num2str(actual(1)), ' mol/L'],
        ['KHCO₃: ', num2str(actual(2)), ' mol/L'],
        '--- Predicted ---',
        ['K₂CO₃: ', num2str(predicted(1)), ' mol/L (Δ: ', num2str(error_k2co3), ')'],
        ['KHCO₃: ', num2str(predicted(2)), ' mol/L (Δ: ', num2str(error_khco3), ')']
    };
end