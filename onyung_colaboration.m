%% Load Data and Preprocess (same as before)
[filename, pathname] = uigetfile('*.xlsx', 'Select the Excel file');
filepath = fullfile(pathname, filename);
data = readtable(filepath, 'Sheet', 'Sheet1');

% Fix column names
originalColNames = data.Properties.VariableNames;
sanitizedColNames = matlab.lang.makeValidName(originalColNames);
data.Properties.VariableNames = sanitizedColNames;

% Extract metadata and spectral data
validRows = ~cellfun(@isempty, data.name);
metadata = data(validRows, :);
spectralData = data(~validRows, :);

% Get concentrations
k2co3_col = find(contains(originalColNames, 'c_K2CO3'));
khco3_col = find(contains(originalColNames, 'c_KHCO3'));
concentration_K2CO3 = metadata.(sanitizedColNames{k2co3_col});
concentration_KHCO3 = metadata.(sanitizedColNames{khco3_col});

wavenumbers = spectralData.Wavenumber;
sampleNames = metadata.name;

%% Univariate Calibration at Specific Wavenumbers
% Example: Predict concentration at a single wavenumber (e.g., ~1000 cm⁻¹)
targetWavenumber = 1000; % Adjust to your peak of interest
[~, idx] = min(abs(wavenumbers - targetWavenumber));
absorptionAtPeak = spectralData{idx, sampleNames}';

% Train a simple linear model (e.g., for K₂CO₃)
model_K2CO3 = fitlm(absorptionAtPeak, concentration_K2CO3);
predicted_K2CO3 = predict(model_K2CO3, absorptionAtPeak);

%% Plot Absorption vs. Wavenumber with Dynamic Predictions
figure;
hold on;
colors = lines(numel(sampleNames));

for i = 1:numel(sampleNames)
    absorption = spectralData{:, sampleNames{i}};
    h = plot(wavenumbers, absorption, 'Color', colors(i,:), ...
             'DisplayName', sampleNames{i});
    
    % Store predictions for ALL wavenumbers (computationally intensive)
    % This is a simplified example; adjust for your use case
    set(h, 'UserData', struct(...
        'Actual_K2CO3', concentration_K2CO3(i), ...
        'Actual_KHCO3', concentration_KHCO3(i), ...
        'Predicted_K2CO3', predicted_K2CO3(i) ...
    ));
end

xlabel('Wavenumber (cm^{-1})');
ylabel('Absorption');
title('Absorption Spectra with Dynamic Predictions');
legend('Location', 'bestoutside');
grid on;

%% Configure Data Tips to Show Local Predictions
dcm = datacursormode(gcf);
set(dcm, 'UpdateFcn', @dynamicDatatip);

%% Custom Data Tip Function
function txt = dynamicDatatip(~, event_obj)
    % Get spectral coordinates
    pos = get(event_obj, 'Position');
    wavenumber = pos(1);
    absorption = pos(2);
    
    % Get the line object and its UserData
    targetLine = get(event_obj, 'Target');
    userdata = get(targetLine, 'UserData');
    
    % Example: Predict concentration at this wavenumber (custom logic)
    % Replace this with your own model/equation
    predicted_K2CO3_local = userdata.Predicted_K2CO3 * (absorption / max(get(targetLine, 'YData')));
    
    % Create text
    txt = { 
        ['Wavenumber: ', num2str(wavenumber)], ...
        ['Absorption: ', num2str(absorption)], ...
        '--- Actual ---', ...
        ['K₂CO₃: ', num2str(userdata.Actual_K2CO3), ' mol/L'], ...
        ['KHCO₃: ', num2str(userdata.Actual_KHCO3), ' mol/L'], ...
        '--- Predicted (Local) ---', ...
        ['K₂CO₃: ', num2str(predicted_K2CO3_local), ' mol/L'] ...
    };
end