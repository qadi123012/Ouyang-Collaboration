function spectral_preprocessor_final()
    %%% Robust Spectral Preprocessing Pipeline %%%
    
    %% Load Data
    [filename, pathname] = uigetfile('*.xlsx', 'Select Raw Spectral Excel File');
    if isequal(filename, 0), error('No file selected'); end
    raw_file = fullfile(pathname, filename);
    fprintf('Loading raw data from: %s\n', raw_file);
    raw_data = readtable(raw_file, 'Sheet', 'Sheet1');
    
    %% Process Column Names
    original_colnames = raw_data.Properties.VariableNames;
    sanitized_colnames = matlab.lang.makeValidName(original_colnames);
    raw_data.Properties.VariableNames = sanitized_colnames;
    
    %% Separate Data
    valid_rows = ~cellfun(@isempty, raw_data.name);
    metadata = raw_data(valid_rows, :);
    spectral_raw = raw_data(~valid_rows, :);
    
    %% Extract Components with Validation
    try
        wavenumbers = spectral_raw.Wavenumber;
        sample_names = metadata.name;
        n_wavenumbers = length(wavenumbers);
        assert(n_wavenumbers >= 10, 'Minimum 10 spectral points required');
    catch ME
        error('Data validation failed: %s', ME.message);
    end
    
    %% Adaptive Preprocessing Parameters
    params = struct(...
        'smoothing',        'wavelet', ... % Wavelet transform for denoising
        'baselineCorrection', true, ...
        'normalization',    'snv', ...
        'derivative',       0 ...
    );
    
    fprintf(['Adjusted parameters:\n'...
             '  Smoothing: %s\n'],...
             params.smoothing);
    
    %% Core Processing
    X = spectral_raw{:, matlab.lang.makeValidName(sample_names)}';
    X_clean = preprocessSpectra(X, wavenumbers, params);
    
    %% Save Results
    [~, name] = fileparts(filename);
    output_file = fullfile(pathname, [name '_PREPROCESSED.xlsx']);
    writetable(metadata, output_file, 'Sheet', 'Metadata');
    writetable(array2table([wavenumbers, X_clean'],...
        'VariableNames', ['Wavenumber', matlab.lang.makeValidName(sample_names)']),...
        output_file, 'Sheet', 'Spectra');
    
    fprintf('Successfully processed. Output: %s\n', output_file);
    
    %% Preprocessing Functions with Enhanced Safeguards
    function X_processed = preprocessSpectra(X, wavenumbers, params)
        X_processed = X;
        [n_samples, n_pts] = size(X_processed);
        
        % 1. Wavelet Transform for Denoising
        if strcmpi(params.smoothing, 'wavelet')
            for i = 1:n_samples
                X_processed(i,:) = wdenoise(X_processed(i,:), 'DenoisingMethod', 'SURE', 'Wavelet', 'sym4');
            end
        end
        
        % 2. Baseline Correction
        if params.baselineCorrection
            for i = 1:n_samples
                baseline = ipf_baseline(wavenumbers, X_processed(i,:)');
                X_processed(i,:) = X_processed(i,:) - baseline';
            end
        end
        
        % 3. Normalization (use SNV normalization)
        switch lower(params.normalization)
            case 'snv'
                X_processed = (X_processed - mean(X_processed,2)) ./ std(X_processed,0,2);
            case 'area'
                area = trapz(wavenumbers, X_processed, 2);
                X_processed = X_processed ./ area;
        end
        
        % 4. Derivatives (only if enough points)
        if params.derivative > 0 && n_pts > 10
            [~, g] = sgolay(2, min(11, n_pts));
            for i = 1:n_samples
                X_processed(i,:) = conv(X_processed(i,:), factorial(params.derivative) * g(:,params.derivative+1), 'same');
            end
        end
    end

    function baseline = ipf_baseline(x, y)
        % Iterative Polynomial Fitting (IPF) for baseline correction
        max_iter = 100;
        tol = 1e-6;
        p = polyfit(x, y, 2);
        baseline = polyval(p, x);
        
        for iter = 1:max_iter
            residual = y - baseline;
            std_res = std(residual);
            idx = abs(residual) < 2 * std_res;
            p = polyfit(x(idx), y(idx), 2);
            new_baseline = polyval(p, x);
            
            if norm(new_baseline - baseline) < tol
                break;
            end
            
            baseline = new_baseline;
        end
    end
end