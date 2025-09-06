%% GLM-based Snow Impact Analysis
% This script analyzes traffic response to snow using a Generalized Linear Model
% with lagged precipitation regressors, constant, and linear trend
% Run after runTelraamAnalysis.m

clc ; close all

%% Check prerequisites
requiredVars = {'locationData', 'weatherData', 'analysis', 'style'};
missingVars = {};

for i = 1:length(requiredVars)
    if ~exist(requiredVars{i}, 'var')
        missingVars{end+1} = requiredVars{i};
    end
end

if ~isempty(missingVars)
    error('Missing required variables: %s\nPlease run runTelraamAnalysis.m first.', ...
        strjoin(missingVars, ', '));
end

%% Configuration
config = struct();
config.winterStart = datetime(2024, 11, 16);
config.winterEnd = datetime(2025, 3, 31);
config.snowTempThreshold = 0; % Temperature threshold for snow (°C)
config.nLags = 10; % Number of lag days for precipitation
config.crossoverDay = 68; % Day where linear trend crosses zero

%% Main Analysis
fprintf('\n========== GLM-BASED SNOW IMPACT ANALYSIS ==========\n\n');
fprintf('Analysis Period: %s to %s\n', ...
    datestr(config.winterStart, 'mmm dd, yyyy'), ...
    datestr(config.winterEnd, 'mmm dd, yyyy'));

% Calculate total days in period
winterDays = days(config.winterEnd - config.winterStart) + 1;
fprintf('Total days in analysis period: %d\n\n', winterDays);

% Perform GLM analysis for all locations and modes
glmResults = performGLMAnalysis(locationData, weatherData, config);

% Generate visualizations
generateGLMVisualizations(glmResults, config, style);

% Generate summary report
generateGLMSummaryReport(glmResults, config);

%% ======================== ANALYSIS FUNCTIONS ========================

function glmResults = performGLMAnalysis(locationData, weatherData, config)
    % Perform GLM analysis for all locations and modes
    
    glmResults = struct();
    locationNames = fieldnames(locationData);
    modes = {'Bike Total', 'Car Total', 'Pedestrian Total'};
    modeNames = {'Bikes', 'Cars', 'Pedestrians'};
    
    % First, prepare winter weather data
    winterWeather = prepareWinterWeatherData(weatherData, config);
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        data = locationData.(locationName);
        
        fprintf('Analyzing Location: %s\n', extractLocationShortName(locationName));
        fprintf('%s\n', repmat('-', 1, 60));
        
        glmResults.(locationName) = struct();
        
        for modeIdx = 1:length(modes)
            currentMode = modes{modeIdx};
            modeName = modeNames{modeIdx};
            
            fprintf('  %s: ', modeName);
            
            % Get daily counts for winter period
            winterCounts = getWinterDailyCounts(data, winterWeather, currentMode);
            
            if length(winterCounts) < length(winterWeather.dates)
                fprintf('Insufficient data\n');
                continue;
            end
            
            % Build design matrix
            X = buildDesignMatrix(winterWeather, config);
            
            % Perform GLM fit
            [beta, stats] = performGLMFit(X, winterCounts);
            
            % Store results
            results = struct();
            results.beta = beta;
            results.stats = stats;
            results.snowResponse = beta(1:config.nLags); % First N betas are snow response
            results.constant = beta(config.nLags + 1);
            results.trend = beta(config.nLags + 2);
            results.dates = winterWeather.dates;
            results.counts = winterCounts;
            results.fitted = X * beta;
            results.residuals = winterCounts - results.fitted;
            
            glmResults.(locationName).(modeName) = results;
            
            % Report key findings
            fprintf('R² = %.3f, ', stats.rsquared);
            
            % Find peak response
            [minResponse, minLag] = min(results.snowResponse);
            fprintf('Peak: %.1f counts/mm at lag %d\n', minResponse, minLag-1);
        end
        fprintf('\n');
    end
end

function winterWeather = prepareWinterWeatherData(weatherData, config)
    % Extract and prepare weather data for winter period
    
    % Find winter period indices
    winterMask = (weatherData.dates >= config.winterStart) & ...
                 (weatherData.dates <= config.winterEnd);
    
    winterWeather = struct();
    winterWeather.dates = weatherData.dates(winterMask);
    winterWeather.temperature = weatherData.temperature(winterMask);
    winterWeather.precipitation = weatherData.precipitation(winterMask);
    
    % Create snow indicator (precipitation when cold)
    winterWeather.snow = winterWeather.precipitation .* ...
                        (winterWeather.temperature < config.snowTempThreshold);
    
    % Report snow statistics
    snowDays = winterWeather.snow > 0;
    fprintf('Snow Statistics:\n');
    fprintf('  Days with snow: %d (%.1f%%)\n', sum(snowDays), 100*sum(snowDays)/length(snowDays));
    fprintf('  Total snowfall: %.1f mm\n', sum(winterWeather.snow));
    fprintf('  Max daily snow: %.1f mm\n\n', max(winterWeather.snow));
end

function winterCounts = getWinterDailyCounts(locationData, winterWeather, modeString)
    % Get daily traffic counts aligned with winter weather data
    
    % Calculate daily totals
    data = locationData.data;
    data.DayOnly = dateshift(data.('Date and Time (Local)'), 'start', 'day');
    
    % Check if the mode exists
    if ~ismember(modeString, data.Properties.VariableNames)
        winterCounts = [];
        return;
    end
    
    groupedData = groupsummary(data, 'DayOnly', 'sum', modeString);
    
    sumColumnName = ['sum_' modeString];
    dailyData = table(groupedData.DayOnly, groupedData.(sumColumnName), ...
        'VariableNames', {'Date', 'Count'});
    
    % Align with winter dates
    winterCounts = zeros(size(winterWeather.dates));
    for i = 1:length(winterWeather.dates)
        idx = find(dailyData.Date == winterWeather.dates(i), 1);
        if ~isempty(idx)
            winterCounts(i) = dailyData.Count(idx);
        else
            winterCounts(i) = NaN;
        end
    end
    
    % Interpolate missing values
    nanIdx = isnan(winterCounts);
    if sum(~nanIdx) > 10 && sum(nanIdx) > 0
        winterCounts(nanIdx) = interp1(find(~nanIdx), winterCounts(~nanIdx), find(nanIdx), 'linear', 'extrap');
    end
end

function X = buildDesignMatrix(winterWeather, config)
    % Build GLM design matrix with lagged precipitation, constant, and trend
    
    n = length(winterWeather.dates);
    nRegressors = config.nLags + 2; % Lags + constant + trend
    X = zeros(n, nRegressors);
    
    % 1. Lagged precipitation regressors (using snow indicator)
    snow = winterWeather.snow;
    for lag = 0:(config.nLags-1)
        if lag == 0
            X(:, lag+1) = snow;
        else
            % Shift by lag days
            X((lag+1):end, lag+1) = snow(1:(end-lag));
            % Pad beginning with zeros
            X(1:lag, lag+1) = 0;
        end
    end
    
    % 2. Constant regressor
    X(:, config.nLags + 1) = 1;
    
    % 3. Linear trend regressor
    % Linear function from -1 to +1, crossing zero at sample config.crossoverDay
    t = (1:n)';
    % Adjust crossover day if it exceeds the number of samples
    actualCrossover = min(config.crossoverDay, n);
    X(:, config.nLags + 2) = 2 * (t - actualCrossover) / n;
end

function [beta, stats] = performGLMFit(X, y)
    % Perform GLM fit using ordinary least squares
    
    % Remove any rows with NaN in y
    validRows = ~isnan(y);
    X = X(validRows, :);
    y = y(validRows);
    
    % Check for sufficient data
    if length(y) < size(X, 2) * 2
        beta = nan(size(X, 2), 1);
        stats = struct('rsquared', NaN, 'pvalues', nan(size(X, 2), 1), ...
                      'se', nan(size(X, 2), 1), 'fstat', NaN);
        return;
    end
    
    % Perform regression
    [beta, ~, residuals, ~, stats_raw] = regress(y, X);
    
    % Calculate additional statistics
    yMean = mean(y);
    SST = sum((y - yMean).^2);
    SSE = sum(residuals.^2);
    SSR = SST - SSE;
    
    stats = struct();
    stats.rsquared = 1 - SSE/SST;
    stats.adjrsquared = 1 - (SSE/(length(y)-length(beta)))/(SST/(length(y)-1));
    
    % Standard errors and p-values
    MSE = SSE / (length(y) - length(beta));
    C = inv(X' * X);
    stats.se = sqrt(MSE * diag(C));
    stats.tstat = beta ./ stats.se;
    stats.pvalues = 2 * (1 - tcdf(abs(stats.tstat), length(y) - length(beta)));
    
    % F-statistic
    stats.fstat = (SSR / (length(beta) - 1)) / MSE;
    stats.fpvalue = 1 - fcdf(stats.fstat, length(beta) - 1, length(y) - length(beta));
end

%% ======================== VISUALIZATION FUNCTIONS ========================

function generateGLMVisualizations(glmResults, config, style)
    % Generate all GLM visualization plots
    
    % 1. Snow response functions
    plotSnowResponseFunctions(glmResults, config, style);
    
    % 2. Model fit comparison
    plotModelFitComparison(glmResults, style);
    
    % 3. Example time series with fit
    plotExampleTimeSeries(glmResults, config, style);
    
    % 4. Coefficient significance
    plotCoefficientSignificance(glmResults, config, style);
end

function plotSnowResponseFunctions(glmResults, config, style)
    % Plot the estimated snow response functions (beta coefficients)
    % Normalized as percentage of baseline (constant term)
    
    figure('Position', [50 50 1200 800]);
    
    locationNames = fieldnames(glmResults);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    colors = {[0 0 1], [1 0 0], [0 0.7 0]};
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        locData = glmResults.(locationName);
        
        subplot(2, length(locationNames), locIdx);
        hold on;
        
        for modeIdx = 1:length(modes)
            if isfield(locData, modes{modeIdx}) && ~isempty(locData.(modes{modeIdx}))
                results = locData.(modes{modeIdx});
                
                if ~isnan(results.snowResponse(1)) && results.constant > 0
                    lags = 0:(config.nLags-1);
                    
                    % Normalize snow response by constant term (baseline)
                    normalizedResponse = (results.snowResponse / results.constant) * 100;
                    
                    h = plot(lags, normalizedResponse, '-o', ...
                        'Color', colors{modeIdx}, ...
                        'LineWidth', 2, ...
                        'MarkerSize', 6, ...
                        'DisplayName', modes{modeIdx});
                    
                    % Add error bars if available (also normalized)
                    if isfield(results.stats, 'se') && length(results.stats.se) >= config.nLags
                        normalizedError = (1.96 * results.stats.se(1:config.nLags) / results.constant) * 100;
                        errorbar(lags, normalizedResponse, normalizedError, ...
                            'Color', colors{modeIdx}, ...
                            'LineStyle', 'none', ...
                            'HandleVisibility', 'off');
                    end
                end
            end
        end
        
        xlabel('Lag (days)');
        ylabel('Response (% of baseline per mm snow)');
        title(sprintf('Snow Response - %s', extractLocationShortName(locationName)));
        grid on;
        legend('Location', 'best');
        xlim([-0.5 config.nLags-0.5]);
        
        % Add zero line
        line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
        
        % Bottom panel: Cumulative response (also normalized)
        subplot(2, length(locationNames), locIdx + length(locationNames));
        hold on;
        
        for modeIdx = 1:length(modes)
            if isfield(locData, modes{modeIdx}) && ~isempty(locData.(modes{modeIdx}))
                results = locData.(modes{modeIdx});
                
                if ~isnan(results.snowResponse(1)) && results.constant > 0
                    lags = 0:(config.nLags-1);
                    
                    % Normalize cumulative response by constant term
                    cumResponse = cumsum(results.snowResponse / results.constant) * 100;
                    
                    plot(lags, cumResponse, '-', ...
                        'Color', colors{modeIdx}, ...
                        'LineWidth', 2, ...
                        'DisplayName', modes{modeIdx});
                end
            end
        end
        
        xlabel('Lag (days)');
        ylabel('Cumulative Response (% of baseline)');
        title('Cumulative Impact');
        grid on;
        xlim([-0.5 config.nLags-0.5]);
        line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
    end
    
    sgtitle('GLM Snow Response Functions (Normalized)', 'FontSize', style.titleFontSize);
end

function plotModelFitComparison(glmResults, style)
    % Compare model fit across locations and modes
    
    figure('Position', [100 100 1000 600]);
    
    locationNames = fieldnames(glmResults);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    
    % Collect R-squared values
    rsquaredMatrix = [];
    labels = {};
    
    for modeIdx = 1:length(modes)
        for locIdx = 1:length(locationNames)
            if isfield(glmResults.(locationNames{locIdx}), modes{modeIdx})
                results = glmResults.(locationNames{locIdx}).(modes{modeIdx});
                if ~isempty(results) && isfield(results.stats, 'rsquared')
                    rsquaredMatrix(modeIdx, locIdx) = results.stats.rsquared;
                else
                    rsquaredMatrix(modeIdx, locIdx) = NaN;
                end
            else
                rsquaredMatrix(modeIdx, locIdx) = NaN;
            end
        end
    end
    
    % Create grouped bar chart
    subplot(1,2,1);
    b = bar(rsquaredMatrix);
    xlabel('Mode');
    ylabel('R²');
    title('Model Fit Comparison');
    set(gca, 'XTickLabel', modes);
    legend(cellfun(@(x) extractLocationShortName(x), locationNames, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;
    ylim([0 max(rsquaredMatrix(:))*1.1]);
    
    % Plot peak responses (as percentage of baseline)
    subplot(1,2,2);
    peakMatrix = [];
    
    for modeIdx = 1:length(modes)
        for locIdx = 1:length(locationNames)
            if isfield(glmResults.(locationNames{locIdx}), modes{modeIdx})
                results = glmResults.(locationNames{locIdx}).(modes{modeIdx});
                if ~isempty(results) && ~isnan(results.snowResponse(1)) && results.constant > 0
                    % Normalize peak response by baseline
                    peakMatrix(modeIdx, locIdx) = (min(results.snowResponse) / results.constant) * 100;
                else
                    peakMatrix(modeIdx, locIdx) = 0;
                end
            else
                peakMatrix(modeIdx, locIdx) = 0;
            end
        end
    end
    
    b = bar(peakMatrix);
    xlabel('Mode');
    ylabel('Peak Response (% per mm)');
    title('Peak Snow Impact');
    set(gca, 'XTickLabel', modes);
    grid on;
    
    sgtitle('GLM Model Performance', 'FontSize', style.titleFontSize);
end

function plotExampleTimeSeries(glmResults, config, style)
    % Plot example time series with actual vs fitted values
    
    figure('Position', [150 150 1200 900]);
    
    locationNames = fieldnames(glmResults);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    
    plotIdx = 1;
    
    % Focus on a 30-day window with significant snow
    windowStart = datetime(2025, 2, 1);
    windowEnd = datetime(2025, 3, 1);
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        
        for modeIdx = 1:length(modes)
            if isfield(glmResults.(locationName), modes{modeIdx})
                results = glmResults.(locationName).(modes{modeIdx});
                
                if ~isempty(results) && ~all(isnan(results.fitted))
                    subplot(length(modes), length(locationNames), plotIdx);
                    
                    % Find window indices
                    windowMask = (results.dates >= windowStart) & (results.dates <= windowEnd);
                    
                    if sum(windowMask) > 10
                        dates = results.dates(windowMask);
                        actual = results.counts(windowMask);
                        fitted = results.fitted(windowMask);
                        
                        % Plot actual vs fitted
                        plot(dates, actual, 'k-', 'LineWidth', 2, 'DisplayName', 'Actual');
                        hold on;
                        plot(dates, fitted, 'r--', 'LineWidth', 2, 'DisplayName', 'GLM Fitted');
                        
                        xlabel('Date');
                        ylabel('Counts');
                        title(sprintf('%s - %s', extractLocationShortName(locationName), modes{modeIdx}));
                        legend('Location', 'best');
                        grid on;
                        xtickformat('MMM dd');
                    end
                end
            end
            
            plotIdx = plotIdx + 1;
        end
    end
    
    sgtitle('GLM Fit Example: February 2025', 'FontSize', style.titleFontSize);
end

function plotCoefficientSignificance(glmResults, config, style)
    % Plot coefficient significance for snow response terms
    
    figure('Position', [200 200 1200 600]);
    
    locationNames = fieldnames(glmResults);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    colors = {[0 0 1], [1 0 0], [0 0.7 0]};
    
    for locIdx = 1:length(locationNames)
        subplot(1, length(locationNames), locIdx);
        hold on;
        
        locationName = locationNames{locIdx};
        
        for modeIdx = 1:length(modes)
            if isfield(glmResults.(locationName), modes{modeIdx})
                results = glmResults.(locationName).(modes{modeIdx});
                
                if ~isempty(results) && isfield(results.stats, 'pvalues')
                    lags = 0:(config.nLags-1);
                    pvals = results.stats.pvalues(1:config.nLags);
                    
                    % Plot -log10(p-value)
                    negLogP = -log10(pvals);
                    plot(lags, negLogP, '-o', ...
                        'Color', colors{modeIdx}, ...
                        'LineWidth', 1.5, ...
                        'MarkerSize', 6, ...
                        'DisplayName', modes{modeIdx});
                end
            end
        end
        
        % Add significance thresholds
        line(xlim, -log10([0.05 0.05]), 'Color', 'k', 'LineStyle', '--', ...
            'DisplayName', 'p=0.05');
        line(xlim, -log10([0.01 0.01]), 'Color', 'k', 'LineStyle', ':', ...
            'DisplayName', 'p=0.01');
        
        xlabel('Lag (days)');
        ylabel('-log10(p-value)');
        title(sprintf('Coefficient Significance - %s', extractLocationShortName(locationName)));
        grid on;
        xlim([-0.5 config.nLags-0.5]);
        
        if locIdx == 1
            legend('Location', 'best');
        end
    end
    
    sgtitle('Statistical Significance of Snow Response Coefficients', 'FontSize', style.titleFontSize);
end

%% ======================== UTILITY FUNCTIONS ========================

function shortName = extractLocationShortName(fullName)
    if contains(fullName, 'Draper')
        shortName = 'Draper';
    elseif contains(fullName, 'King')
        shortName = 'King Edward';
    else
        shortName = fullName;
    end
end

function generateGLMSummaryReport(glmResults, config)
    % Generate summary report of GLM findings
    
    fprintf('\n========== GLM ANALYSIS SUMMARY ==========\n\n');
    
    locationNames = fieldnames(glmResults);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    
    % 1. Model fit summary
    fprintf('1. MODEL FIT (R²):\n');
    for locIdx = 1:length(locationNames)
        fprintf('   %s:\n', extractLocationShortName(locationNames{locIdx}));
        for modeIdx = 1:length(modes)
            if isfield(glmResults.(locationNames{locIdx}), modes{modeIdx})
                results = glmResults.(locationNames{locIdx}).(modes{modeIdx});
                if ~isempty(results) && isfield(results.stats, 'rsquared')
                    fprintf('     %s: %.3f\n', modes{modeIdx}, results.stats.rsquared);
                end
            end
        end
    end
    
    % 2. Peak impact summary
    fprintf('\n2. PEAK SNOW IMPACT (counts per mm):\n');
    for modeIdx = 1:length(modes)
        fprintf('   %s:\n', modes{modeIdx});
        for locIdx = 1:length(locationNames)
            if isfield(glmResults.(locationNames{locIdx}), modes{modeIdx})
                results = glmResults.(locationNames{locIdx}).(modes{modeIdx});
                if ~isempty(results) && ~isnan(results.snowResponse(1))
                    [minVal, minLag] = min(results.snowResponse);
                    fprintf('     %s: %.1f at lag %d days\n', ...
                        extractLocationShortName(locationNames{locIdx}), minVal, minLag-1);
                end
            end
        end
    end
    
    % 3. Recovery time analysis
    fprintf('\n3. RECOVERY PATTERNS:\n');
    for modeIdx = 1:length(modes)
        fprintf('   %s:\n', modes{modeIdx});
        
        % Average cumulative impact across locations
        cumImpacts = [];
        for locIdx = 1:length(locationNames)
            if isfield(glmResults.(locationNames{locIdx}), modes{modeIdx})
                results = glmResults.(locationNames{locIdx}).(modes{modeIdx});
                if ~isempty(results) && ~isnan(results.snowResponse(1))
                    cumImpacts(locIdx, :) = cumsum(results.snowResponse);
                end
            end
        end
        
        if ~isempty(cumImpacts)
            avgCumImpact = mean(cumImpacts, 1);
            totalImpact = avgCumImpact(end);
            
            % Find 90% recovery point
            recoveryThreshold = 0.9 * totalImpact;
            recoveryIdx = find(avgCumImpact <= recoveryThreshold, 1, 'last');
            
            fprintf('     Total impact: %.0f count-days per mm\n', totalImpact);
            if ~isempty(recoveryIdx)
                fprintf('     90%% recovery: %d days\n', recoveryIdx);
            else
                fprintf('     90%% recovery: <1 day\n');
            end
        end
    end
    
    % 4. Trend analysis
    fprintf('\n4. SEASONAL TRENDS:\n');
    for locIdx = 1:length(locationNames)
        fprintf('   %s:\n', extractLocationShortName(locationNames{locIdx}));
        for modeIdx = 1:length(modes)
            if isfield(glmResults.(locationNames{locIdx}), modes{modeIdx})
                results = glmResults.(locationNames{locIdx}).(modes{modeIdx});
                if ~isempty(results) && isfield(results, 'trend')
                    % Calculate total change over period
                    totalChange = results.trend * 2 * mean(results.counts);
                    fprintf('     %s: %.0f counts (%.1f%% of average)\n', ...
                        modes{modeIdx}, totalChange, ...
                        100 * results.trend * 2);
                end
            end
        end
    end
    
    fprintf('\n==========================================\n');
end