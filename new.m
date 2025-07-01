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

%% Classical Least Squares to Estimate Pure Component Spectra
% Estimate pure spectra (S) using concentration data (Y) and spectral data (X)
S = Y \ X;  % Pure component spectra (rows: K2CO3, KHCO3)
pure_K2CO3 = S(1, :);
pure_KHCO3 = S(2, :);

% Reconstruct spectra from pure components and concentrations
X_reconstructed = Y * S;

%% Preprocess Data (Mean Centering)
X_mean = mean(X, 1);
Y_mean = mean(Y, 1);
X_centered = X - X_mean;
Y_centered = Y - Y_mean;

%% Cross-Validation to Determine Optimal PLS Components
max_components = 10;
k_folds = 5;
mse = zeros(max_components, 1);

rng(0); % Set random seed for reproducibility

for ncomp = 1:max_components
    cv_mse = 0;
    if license('test', 'statistics_toolbox')
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
title('Estimated Pure Component Spectra');
legend('K₂CO₃', 'KHCO₃', 'Location', 'best');
grid on;
hold off;

%% Plot Measured vs. Reconstructed Spectra
figure;
hold on;
colors = lines(numel(sampleNames));

for i = 1:numel(sampleNames)
    % Plot measured spectrum
    absorption = spectralData{:, sanitizedSampleNames{i}};
    h = plot(wavenumbers, absorption, 'Color', colors(i,:), 'DisplayName', sampleNames{i});
    
    % Plot reconstructed spectrum
      %reconstructed = Y(i,:) * S;
     % plot(wavenumbers, reconstructed, '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Store data for tooltips
    set(h, 'UserData', {[concentration_K2CO3(i), concentration_KHCO3(i)], Y_pred(i,:)});
end

xlabel('Wavenumber (cm^{-1})');
ylabel('Absorption');
title('Spectra of K₂CO₃ and KHCO₃ Concentration ');
legend('Location', 'bestoutside');
grid on;

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