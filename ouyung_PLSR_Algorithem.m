%% Load Data
[filename, pathname] = uigetfile('*.xlsx', 'Select the Excel file');
filepath = fullfile(pathname, filename);
data = readtable(filepath, 'Sheet', 'Sheet1');

%% Fix Column Names
originalColNames = data.Properties.VariableNames;
sanitizedColNames = matlab.lang.makeValidName(originalColNames);
data.Properties.VariableNames = sanitizedColNames;   

%% Extract Metadata and Spectral Data      
validRows = ~cellfun(@isempty, data.name); 
metadata = data(validRows, :);
spectralData = data(~validRows, :);  

% Get concentrations
k2co3_col = find(contains(originalColNames, 'c_K2CO3'));
khco3_col = find(contains(originalColNames, 'c_KHCO3'));

concentration_K2CO3 = metadata.(sanitizedColNames{k2co3_col});
concentration_KHCO3 = metadata.(sanitizedColNames{khco3_col});

% Extract wavenumbers and sample names
wavenumbers = spectralData.Wavenumber;
sampleNames = metadata.name;

%% Correctly Extract Spectral Data (X)
sanitizedSampleNames = matlab.lang.makeValidName(sampleNames);
X = spectralData{:, sanitizedSampleNames}';  % Samples x Wavenumbers
Y = [concentration_K2CO3, concentration_KHCO3];

%% Preprocess Data (Mean Centering)
X_mean = mean(X, 1);
Y_mean = mean(Y, 1);
X_centered = X - X_mean;
Y_centered = Y - Y_mean;

%% Cross-Validation to Determine Optimal PLS Components
max_components = 10;
k_folds = 5;
mse = zeros(max_components, 1);

% Check for Statistics Toolbox
hasStatsToolbox = license('test', 'statistics_toolbox');

% Set random seed for reproducibility in cross-validation splits
rng(0); % <--- Critical fix: Ensures consistent random splits every run

for ncomp = 1:max_components
    cv_mse = 0;
    
    % K-fold implementation
    if hasStatsToolbox
        % Use cvpartition with fixed seed
        cv = cvpartition(size(X_centered, 1), 'KFold', k_folds);
        for fold = 1:k_folds
            train = cv.training(fold);
            test = cv.test(fold);
            
            % Train PLS model
            [~, ~, ~, ~, beta_cv] = plsregress(X_centered(train, :), Y_centered(train, :), ncomp);
            
            % Predict
            Y_pred_cv = [ones(sum(test),1), X_centered(test, :)] * beta_cv;
            
            % Calculate MSE
            cv_mse = cv_mse + sum((Y_pred_cv - Y_centered(test, :)).^2, 'all') / (k_folds * sum(test));
        end
    else
        % Manual k-fold splitting (now uses fixed seed from rng(0))
        n_samples = size(X_centered, 1);
        indices = randperm(n_samples); % Reproducible due to fixed seed
        fold_size = floor(n_samples / k_folds);
        
        for fold = 1:k_folds
            test_start = (fold-1)*fold_size + 1;
            test_end = min(fold*fold_size, n_samples);
            test_indices = indices(test_start:test_end);
            
            test = false(n_samples, 1);
            test(test_indices) = true;
            train = ~test;
            
            % Train PLS model
            [~, ~, ~, ~, beta_cv] = plsregress(X_centered(train, :), Y_centered(train, :), ncomp);
            
            % Predict
            Y_pred_cv = [ones(sum(test),1), X_centered(test, :)] * beta_cv;
            
            % Calculate MSE
            cv_mse = cv_mse + sum((Y_pred_cv - Y_centered(test, :)).^2, 'all') / (k_folds * sum(test));
        end
    end
    
    mse(ncomp) = cv_mse;
end

[~, optimal_ncomp] = min(mse);
fprintf('Optimal PLS components: %d\n', optimal_ncomp);

%% Train Final PLS Model (No changes needed here)
[~, ~, ~, ~, beta] = plsregress(X_centered, Y_centered, optimal_ncomp);
Y_pred_centered = [ones(size(X_centered,1),1), X_centered] * beta;
Y_pred = Y_pred_centered + Y_mean;

%% Model Evaluation (No changes needed here)
% Calculate R² and RMSE
r2_k2co3 = 1 - sum((Y(:,1) - Y_pred(:,1)).^2) / sum((Y(:,1) - mean(Y(:,1))).^2);
rmse_k2co3 = sqrt(mean((Y(:,1) - Y_pred(:,1)).^2));

r2_khco3 = 1 - sum((Y(:,2) - Y_pred(:,2)).^2) / sum((Y(:,2) - mean(Y(:,2))).^2);
rmse_khco3 = sqrt(mean((Y(:,2) - Y_pred(:,2)).^2));

fprintf('K₂CO₃ Performance:\nR² = %.4f\nRMSE = %.4f\n', r2_k2co3, rmse_k2co3);
fprintf('KHCO₃ Performance:\nR² = %.4f\nRMSE = %.4f\n', r2_khco3, rmse_khco3);

%% Plot Results with Enhanced Tooltips (No changes needed here)
% ... (remaining plotting code remains unchanged)

%% Plot Results with Enhanced Tooltips
figure;
hold on;
colors = lines(numel(sampleNames));

for i = 1:numel(sampleNames)
    absorption = spectralData{:, sanitizedSampleNames{i}};
    h = plot(wavenumbers, absorption, ...
             'Color', colors(i,:), ...
             'DisplayName', sampleNames{i});
    
    % Store actual and predicted values
    set(h, 'UserData', {[concentration_K2CO3(i), concentration_KHCO3(i)], Y_pred(i,:)});
end

xlabel('Wavenumber (cm^{-1})');
ylabel('Absorption');
title('Spectra with Concentration Predictions');
legend('Location', 'bestoutside');
grid on;
%% Plot Actual vs. Predicted with Parity Line
figure;

% K₂CO₃ subplot
subplot(1,2,1);
scatter(Y(:,1), Y_pred(:,1), 40, 'filled', 'MarkerFaceColor', [0 0.447 0.741]);
hold on;
minVal = min([Y(:,1); Y_pred(:,1)]);
maxVal = max([Y(:,1); Y_pred(:,1)]);
plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 1.5);
xlabel('Actual K₂CO₃ Concentration (mol/L)');
ylabel('Predicted K₂CO₃ Concentration (mol/L)');
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
xlabel('Actual KHCO₃ Concentration (mol/L)');
ylabel('Predicted KHCO₃ Concentration (mol/L)');
title(sprintf('KHCO₃\nR² = %.4f, RMSE = %.4f', r2_khco3, rmse_khco3));
grid on;
legend('Predictions', '1:1 Line', 'Location', 'southeast');
axis equal tight

%% Plot Residuals
%figure;

% K₂CO₃ residuals
% subplot(1,2,1);
% residuals_k2co3 = Y_pred(:,1) - Y(:,1);
% scatter(Y(:,1), residuals_k2co3, 40, 'filled', 'MarkerFaceColor', [0 0.447 0.741]);
% hold on;
% plot([min(Y(:,1)) max(Y(:,1))], [0 0], 'k--', 'LineWidth', 1.5);
% xlabel('Actual K₂CO₃ Concentration (mol/L)');
% ylabel('Residual (Predicted - Actual)');
% title('K₂CO₃ Residual Plot');
% grid on;

% KHCO₃ residuals
% subplot(1,2,2);
% residuals_khco3 = Y_pred(:,2) - Y(:,2);
% scatter(Y(:,2), residuals_khco3, 40, 'filled', 'MarkerFaceColor', [0.85 0.325 0.098]);
% hold on;
% plot([min(Y(:,2)) max(Y(:,2))], [0 0], 'k--', 'LineWidth', 1.5);
% xlabel('Actual KHCO₃ Concentration (mol/L)');
% ylabel('Residual (Predicted - Actual)');
% title('KHCO₃ Residual Plot');
% grid on;

%% Configure Data Tips
dcm = datacursormode(gcf);
set(dcm, 'UpdateFcn', @enhancedDatatip);

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