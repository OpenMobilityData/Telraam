%% Decomposed Weather Effects Analysis for Telraam Data
% This script performs advanced weather analysis that separates seasonal
% effects from day-to-day weather anomalies.
%
% PREREQUISITE: Run runTelraamAnalysis.m first to generate required data
%
% Required variables from workspace:
%   - locationData: processed traffic count data for all locations
%   - weatherData: weather data including temperature, precipitation, etc.
%   - analysis: analysis parameters including time range and mode
%   - style: plotting style parameters

%% Check for required variables
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

fprintf('\n========== DECOMPOSED WEATHER EFFECTS ANALYSIS ==========\n');
fprintf('Separating seasonal patterns from weather anomalies\n\n');

%% Perform the decomposed analysis
performDecomposedWeatherAnalysis(locationData, weatherData, analysis, style);

%% ======================== MAIN ANALYSIS FUNCTION ========================

function performDecomposedWeatherAnalysis(locationData, weatherData, analysis, style)
    % Perform analysis that separates seasonal from day-to-day weather effects
    
    fprintf('\n=== Decomposed Weather Effects Analysis for %s ===\n', analysis.modeDisplayString);
    
    locationNames = fieldnames(locationData);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        fprintf('\nLocation: %s\n', locationInfo.name);
        fprintf('%s\n', repmat('=', 1, 70));
        
        % Prepare weekly data
        [weeklyData, isValid] = prepareMultivariateData(data, weatherData, analysis);
        
        if sum(isValid) < 20
            fprintf('Insufficient data for analysis (n=%d)\n', sum(isValid));
            continue;
        end
        
        % 1. Calculate seasonal expectations
        seasonalData = calculateSeasonalBaseline(weeklyData, isValid);
        
        % 2. Calculate weather anomalies (deviations from seasonal norms)
        anomalyData = calculateWeatherAnomalies(weeklyData, seasonalData, isValid);
        
        % 3. Perform three separate analyses
        fprintf('\n--- MODEL 1: Seasonal Effects Only ---\n');
        seasonalResults = analyzeSeasonalEffects(weeklyData, seasonalData, isValid);
        
        fprintf('\n--- MODEL 2: Weather Anomalies Only ---\n');
        anomalyResults = analyzeWeatherAnomalies(weeklyData, anomalyData, isValid);
        
        fprintf('\n--- MODEL 3: Combined Model (Seasonal + Anomalies) ---\n');
        combinedResults = analyzeCombinedEffects(weeklyData, seasonalData, anomalyData, isValid);
        
        % 4. Variance decomposition
        performVarianceDecomposition(seasonalResults, anomalyResults, combinedResults, locationInfo);
        
        % 5. Create visualization
        plotDecomposedEffects(weeklyData, seasonalData, anomalyData, isValid, locationInfo, analysis, style);
    end
end

%% ======================== DATA PREPARATION FUNCTIONS ========================

function [weeklyData, isValid] = prepareMultivariateData(locationDataStruct, weatherData, analysis)
    % Prepare weekly aggregated data with all weather variables
    
    % Get weekly traffic counts
    weeklyTraffic = calculateWeeklyTotals(locationDataStruct, analysis);
    
    % Aggregate weather data to weekly
    weeklyWeather = aggregateAllWeatherVariables(weatherData);
    
    % Match weekly traffic with weather data
    [commonWeeks, ia, ib] = intersect(weeklyTraffic.weekStarts, weeklyWeather.weekStarts);
    
    % Initialize output structure
    weeklyData = struct();
    weeklyData.dates = commonWeeks;
    weeklyData.counts = weeklyTraffic.rawCounts(ia);
    weeklyData.temperature = weeklyWeather.avgTemperature(ib);
    weeklyData.precipitation = weeklyWeather.totalPrecipitation(ib);
    weeklyData.windspeed = weeklyWeather.avgWindspeed(ib);
    
    % Add interaction term: precipitation × temperature
    weeklyData.precipTempInteraction = weeklyData.precipitation .* weeklyData.temperature;
    
    % Add derived variables
    weeklyData.weekOfYear = week(weeklyData.dates);
    weeklyData.monthOfYear = month(weeklyData.dates);
    
    % Calculate daylight hours (approximate based on sunrise/sunset)
    if isfield(weeklyWeather, 'avgDaylightHours')
        weeklyData.daylightHours = weeklyWeather.avgDaylightHours(ib);
    else
        % Approximate daylight hours based on month (for Montreal)
        weeklyData.daylightHours = 9 + 6 * sin((weeklyData.monthOfYear - 3) * pi / 6);
    end
    
    % Remove rows with any NaN values
    isValid = ~any(isnan([weeklyData.counts, weeklyData.temperature, ...
                         weeklyData.precipitation, weeklyData.windspeed]), 2);
end

function weeklyData = calculateWeeklyTotals(locationDataStruct, analysis)
    % Calculate weekly totals - copied from main script for independence
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Calculate weekly totals using the existing yearWeekKey
    groupedData = groupsummary(data, 'yearWeekKey', 'sum', analysis.modeString);
    
    % Get week start dates for each yearWeekKey
    weekStartGrouped = groupsummary(data, 'yearWeekKey', 'min', 'weekStartDateTimes');
    
    % Combine results
    weeklyData = struct();
    weeklyData.yearWeekKeys = groupedData.yearWeekKey;
    
    % Use the original column name with groupsummary's 'sum_' prefix
    sumColumnName = ['sum_' analysis.modeString];
    weeklyData.rawCounts = groupedData.(sumColumnName);
    
    % Get week start dates
    [~, ia, ib] = intersect(groupedData.yearWeekKey, weekStartGrouped.yearWeekKey);
    weekStarts = NaT(size(weeklyData.rawCounts));
    weekStarts(ia) = weekStartGrouped.min_weekStartDateTimes(ib);
    weeklyData.weekStarts = weekStarts;
end

function weeklyWeather = aggregateAllWeatherVariables(weatherData)
    % Aggregate all weather variables to weekly
    
    tempTable = table(weatherData.dates, weatherData.temperature, ...
        weatherData.precipitation, weatherData.windspeed, ...
        'VariableNames', {'dates', 'temperature', 'precipitation', ...
        'windspeed'});
    
    % Add week start dates
    tempTable.weekStarts = dateshift(dateshift(tempTable.dates,'dayofweek','Monday','previous'),'start','day');
    
    % Group by week - using appropriate aggregation for each variable
    weeklyGrouped = groupsummary(tempTable, 'weekStarts', ...
        {'mean', 'sum', 'mean'}, ...
        {'temperature', 'precipitation', 'windspeed'});
    
    weeklyWeather = struct();
    weeklyWeather.weekStarts = weeklyGrouped.weekStarts;
    weeklyWeather.avgTemperature = weeklyGrouped.mean_temperature;
    weeklyWeather.totalPrecipitation = weeklyGrouped.sum_precipitation;
    weeklyWeather.avgWindspeed = weeklyGrouped.mean_windspeed;
end

%% ======================== SEASONAL BASELINE FUNCTIONS ========================

function seasonalData = calculateSeasonalBaseline(weeklyData, isValid)
    % Calculate expected seasonal patterns using centered moving average
    % with circular padding to handle edge effects
    
    seasonalData = struct();
    
    % Option 1: Use daylight hours as seasonal proxy (keep existing)
    seasonalData.daylightHours = weeklyData.daylightHours(isValid);
    
    % Option 2: Use sinusoidal seasonal pattern (keep existing)
    weekOfYear = week(weeklyData.dates(isValid));
    seasonalData.seasonalIndex = sin((weekOfYear - 13) * 2 * pi / 52); % Peak at week 13 (late March)
    
    % Option 3: Calculate seasonal baseline using centered moving average
    window = 4; % 4-week window (about a month)
    
    % Sort data by date to ensure proper ordering
    [sortedDates, sortIdx] = sort(weeklyData.dates);
    sortedCounts = weeklyData.counts(sortIdx);
    sortedValid = isValid(sortIdx);
    
    % Find valid data indices
    validIndices = find(sortedValid);
    if isempty(validIndices)
        error('No valid data points for seasonal baseline calculation');
    end
    
    % Extract only valid data for processing
    validCounts = sortedCounts(validIndices);
    validDates = sortedDates(validIndices);
    
    % Apply circular padding for seasonal continuity
    extend_weeks = ceil(window/2);
    
    % Check if we have enough data points for circular padding
    if length(validCounts) > extend_weeks
        % Circular padding: connect end of year to beginning
        extended_counts = [validCounts(end-extend_weeks+1:end); validCounts; validCounts(1:extend_weeks)];
        
        % Apply moving average to extended data
        extended_avg = movmean(extended_counts, window, 'omitnan');
        
        % Extract the middle portion (original data length)
        seasonal_avg_valid = extended_avg(extend_weeks+1:end-extend_weeks);
    else
        % If not enough data for circular padding, use shrinking window
        seasonal_avg_valid = movmean(validCounts, window, 'omitnan', 'Endpoints', 'shrink');
    end
    
    % Map back to full data array
    allSeasonalAvg = nan(size(sortedCounts));
    allSeasonalAvg(validIndices) = seasonal_avg_valid;
    
    % Map back to original order
    invSortIdx(sortIdx) = 1:length(sortIdx);
    seasonalAverage = allSeasonalAvg(invSortIdx);
    
    % Extract only the valid subset for output
    seasonalData.monthlyAverage = seasonalAverage(isValid);
    
    % Keep month information for compatibility with other functions
    months = month(weeklyData.dates(isValid));
    seasonalData.monthOfYear = months;
end

function anomalyData = calculateWeatherAnomalies(weeklyData, seasonalData, isValid)
    % Calculate deviations from seasonal norms using centered moving average
    % with circular padding to handle edge effects
    
    anomalyData = struct();
    window = 4; % 4-week window (about a month)
    
    % Get dates for sorting
    allDates = weeklyData.dates;
    validDates = allDates(isValid);
    
    %% Temperature anomaly calculation
    
    % Sort data by date
    [sortedDates, sortIdx] = sort(allDates);
    sortedTemp = weeklyData.temperature(sortIdx);
    sortedValid = isValid(sortIdx);
    
    % Find valid data indices
    validIndices = find(sortedValid);
    if ~isempty(validIndices)
        % Extract valid temperatures
        validTemp = sortedTemp(validIndices);
        
        % Apply circular padding
        extend_weeks = ceil(window/2);
        if length(validTemp) > extend_weeks
            extended_temp = [validTemp(end-extend_weeks+1:end); validTemp; validTemp(1:extend_weeks)];
            extended_avg = movmean(extended_temp, window, 'omitnan');
            seasonal_temp_valid = extended_avg(extend_weeks+1:end-extend_weeks);
        else
            seasonal_temp_valid = movmean(validTemp, window, 'omitnan', 'Endpoints', 'shrink');
        end
        
        % Map back to full array
        seasonalTemp = nan(size(sortedTemp));
        seasonalTemp(validIndices) = seasonal_temp_valid;
        
        % Map back to original order
        invSortIdx(sortIdx) = 1:length(sortIdx);
        expectedTemp = seasonalTemp(invSortIdx);
        
        % Calculate anomaly for valid data
        anomalyData.temperatureAnomaly = weeklyData.temperature(isValid) - expectedTemp(isValid);
    else
        anomalyData.temperatureAnomaly = [];
    end
    
    %% Precipitation anomaly calculation (log-transformed)
    
    % Reset sort indices for precipitation
    clear invSortIdx;
    sortedPrecip = weeklyData.precipitation(sortIdx);
    
    if ~isempty(validIndices)
        % Extract valid precipitation (add 1 to avoid log(0))
        validPrecip = sortedPrecip(validIndices) + 1;
        
        % Apply circular padding
        if length(validPrecip) > extend_weeks
            extended_precip = [validPrecip(end-extend_weeks+1:end); validPrecip; validPrecip(1:extend_weeks)];
            extended_avg = movmean(extended_precip, window, 'omitnan');
            seasonal_precip_valid = extended_avg(extend_weeks+1:end-extend_weeks);
        else
            seasonal_precip_valid = movmean(validPrecip, window, 'omitnan', 'Endpoints', 'shrink');
        end
        
        % Map back to full array
        seasonalPrecip = nan(size(sortedPrecip));
        seasonalPrecip(validIndices) = seasonal_precip_valid;
        
        % Map back to original order
        invSortIdx(sortIdx) = 1:length(sortIdx);
        expectedPrecip = seasonalPrecip(invSortIdx);
        
        % Calculate log anomaly for valid data
        anomalyData.precipitationAnomaly = log(weeklyData.precipitation(isValid) + 1) - log(expectedPrecip(isValid));
    else
        anomalyData.precipitationAnomaly = [];
    end
    
    %% Wind speed anomaly calculation
    
    % Reset sort indices for wind
    clear invSortIdx;
    sortedWind = weeklyData.windspeed(sortIdx);
    
    if ~isempty(validIndices)
        % Extract valid wind speeds
        validWind = sortedWind(validIndices);
        
        % Apply circular padding
        if length(validWind) > extend_weeks
            extended_wind = [validWind(end-extend_weeks+1:end); validWind; validWind(1:extend_weeks)];
            extended_avg = movmean(extended_wind, window, 'omitnan');
            seasonal_wind_valid = extended_avg(extend_weeks+1:end-extend_weeks);
        else
            seasonal_wind_valid = movmean(validWind, window, 'omitnan', 'Endpoints', 'shrink');
        end
        
        % Map back to full array
        seasonalWind = nan(size(sortedWind));
        seasonalWind(validIndices) = seasonal_wind_valid;
        
        % Map back to original order
        invSortIdx(sortIdx) = 1:length(sortIdx);
        expectedWind = seasonalWind(invSortIdx);
        
        % Calculate anomaly for valid data
        anomalyData.windspeedAnomaly = weeklyData.windspeed(isValid) - expectedWind(isValid);
    else
        anomalyData.windspeedAnomaly = [];
    end
    
    %% Handle empty anomaly data
    
    if isempty(anomalyData.temperatureAnomaly)
        warning('No valid temperature anomaly data calculated');
        anomalyData.temperatureAnomaly = zeros(sum(isValid), 1);
    end
    
    if isempty(anomalyData.precipitationAnomaly)
        warning('No valid precipitation anomaly data calculated');
        anomalyData.precipitationAnomaly = zeros(sum(isValid), 1);
    end
    
    if isempty(anomalyData.windspeedAnomaly)
        warning('No valid wind speed anomaly data calculated');
        anomalyData.windspeedAnomaly = zeros(sum(isValid), 1);
    end
end

%% ======================== REGRESSION ANALYSIS FUNCTIONS ========================

function results = analyzeSeasonalEffects(weeklyData, seasonalData, isValid)
    % Analyze only seasonal effects
    
    y = weeklyData.counts(isValid);
    
    % Use daylight hours and seasonal index as predictors
    X = [seasonalData.daylightHours, seasonalData.seasonalIndex];
    
    % Standardize and add intercept
    [X_std, mu, sigma] = zscore(X);
    X_std = [ones(size(X_std,1),1), X_std];
    
    % Fit model
    [b, ~, ~, ~, stats] = regress(y, X_std);
    
    results = struct();
    results.coefficients = b;
    results.rSquared = stats(1);
    results.predictedValues = X_std * b;
    results.residuals = y - results.predictedValues;
    
    fprintf('Seasonal R²: %.3f\n', results.rSquared);
    fprintf('Coefficients: Daylight=%.3f, SeasonalIndex=%.3f\n', b(2), b(3));
end

function results = analyzeWeatherAnomalies(weeklyData, anomalyData, isValid)
    % Analyze only weather anomaly effects
    
    y = weeklyData.counts(isValid);
    
    % Use weather anomalies as predictors
    X = [anomalyData.temperatureAnomaly, ...
         anomalyData.precipitationAnomaly, ...
         anomalyData.windspeedAnomaly];
    
    % Add interaction term for temp-precip anomaly
    X = [X, anomalyData.temperatureAnomaly .* anomalyData.precipitationAnomaly];
    
    % Standardize and add intercept
    [X_std, mu, sigma] = zscore(X);
    X_std = [ones(size(X_std,1),1), X_std];
    
    % Fit model
    [b, ~, ~, ~, stats] = regress(y, X_std);
    
    results = struct();
    results.coefficients = b;
    results.rSquared = stats(1);
    results.predictedValues = X_std * b;
    results.residuals = y - results.predictedValues;
    
    fprintf('Weather Anomaly R²: %.3f\n', results.rSquared);
    fprintf('Coefficients: TempAnom=%.3f, PrecipAnom=%.3f, WindAnom=%.3f, Interaction=%.3f\n', ...
        b(2), b(3), b(4), b(5));
end

function results = analyzeCombinedEffects(weeklyData, seasonalData, anomalyData, isValid)
    % Analyze combined seasonal + weather anomaly effects
    
    y = weeklyData.counts(isValid);
    
    % Combine seasonal and anomaly predictors
    X = [seasonalData.daylightHours, ...
         seasonalData.seasonalIndex, ...
         anomalyData.temperatureAnomaly, ...
         anomalyData.precipitationAnomaly, ...
         anomalyData.windspeedAnomaly, ...
         anomalyData.temperatureAnomaly .* anomalyData.precipitationAnomaly];
    
    % Standardize and add intercept
    [X_std, mu, sigma] = zscore(X);
    X_std = [ones(size(X_std,1),1), X_std];
    
    % Fit model
    [b, ~, ~, ~, stats] = regress(y, X_std);
    
    results = struct();
    results.coefficients = b;
    results.rSquared = stats(1);
    results.predictedValues = X_std * b;
    results.residuals = y - results.predictedValues;
    
    fprintf('Combined Model R²: %.3f\n', results.rSquared);
    
    % Calculate importance of each component
    results.importance = abs(b(2:end));
    results.predictorNames = {'Daylight', 'SeasonalIndex', 'TempAnomaly', ...
        'PrecipAnomaly', 'WindAnomaly', 'Temp×Precip'};
end

%% ======================== VARIANCE DECOMPOSITION FUNCTIONS ========================

function performVarianceDecomposition(seasonalResults, anomalyResults, combinedResults, locationInfo)
    % Decompose variance into seasonal vs weather components
    
    fprintf('\n--- Variance Decomposition ---\n');
    
    % Calculate unique contributions
    seasonalOnly = seasonalResults.rSquared;
    weatherOnly = anomalyResults.rSquared;
    combined = combinedResults.rSquared;
    
    % Estimate unique and shared variance
    uniqueSeasonal = max(0, combined - weatherOnly);
    uniqueWeather = max(0, combined - seasonalOnly);
    sharedVariance = max(0, seasonalOnly + weatherOnly - combined);
    unexplained = 1 - combined;
    
    fprintf('Variance explained by:\n');
    fprintf('  Seasonal factors only: %.1f%%\n', uniqueSeasonal * 100);
    fprintf('  Weather anomalies only: %.1f%%\n', uniqueWeather * 100);
    fprintf('  Shared (confounded): %.1f%%\n', sharedVariance * 100);
    fprintf('  Unexplained: %.1f%%\n', unexplained * 100);
    fprintf('  Total: %.1f%%\n', combined * 100);
    
    % Create pie chart of variance decomposition
    figure('Position', [408 126 600 600]);
    
    values = [uniqueSeasonal, uniqueWeather, sharedVariance, unexplained];
    labels = {'Seasonal Only', 'Weather Only', 'Shared', 'Unexplained'};
    colors = [0.2 0.7 0.2;  % Green for seasonal
              0.2 0.4 0.8;  % Blue for weather
              0.8 0.8 0.2;  % Yellow for shared
              0.8 0.8 0.8]; % Gray for unexplained
    
    p = pie(values, labels);
    
    % Customize colors
    for i = 1:4
        set(p(2*i-1), 'FaceColor', colors(i,:));
    end
    
    title(sprintf('Variance Decomposition for %s\n%s', ...
        locationInfo.name, 'Weekly Bike Counts'), 'FontSize', 14);
end

%% ======================== VISUALIZATION FUNCTIONS ========================

function plotDecomposedEffects(weeklyData, seasonalData, anomalyData, isValid, locationInfo, analysis, style)
    % Create visualization of decomposed effects
    % Now uses the smooth seasonal baselines already calculated
    
    figure('Position', [408 126 1200 800]);
    
    dates = weeklyData.dates(isValid);
    counts = weeklyData.counts(isValid);
    
    % Subplot 1: Raw counts with seasonal baseline
    subplot(3,1,1);
    plot(dates, counts, 'k-', 'LineWidth', 1.5);
    hold on;
    plot(dates, seasonalData.monthlyAverage, 'g--', 'LineWidth', 2);
    ylabel('Weekly Counts');
    legend('Actual', 'Seasonal Average', 'Location', 'northwest');
    title(sprintf('%s - %s', analysis.modeDisplayString, locationInfo.name));
    grid on;
    
    % Subplot 2: Temperature actual vs expected
    subplot(3,1,2);
    plot(dates, weeklyData.temperature(isValid), 'r-', 'LineWidth', 1.5);
    hold on;
    
    % Calculate expected temperature from the anomaly data
    % Since anomaly = actual - expected, then expected = actual - anomaly
    expectedTemp = weeklyData.temperature(isValid) - anomalyData.temperatureAnomaly;
    plot(dates, expectedTemp, 'r--', 'LineWidth', 2);
    
    ylabel('Temperature (°C)');
    legend('Actual', 'Seasonal Normal', 'Location', 'northwest');
    grid on;
    
    % Subplot 3: Temperature anomaly effect on counts
    subplot(3,1,3);
    scatter(anomalyData.temperatureAnomaly, counts - seasonalData.monthlyAverage, ...
        50, 'filled', 'MarkerFaceAlpha', 0.6);
    xlabel('Temperature Anomaly (°C)');
    ylabel('Count Deviation from Seasonal');
    
    % Add trend line
    p = polyfit(anomalyData.temperatureAnomaly, counts - seasonalData.monthlyAverage, 1);
    xRange = linspace(min(anomalyData.temperatureAnomaly), max(anomalyData.temperatureAnomaly), 100);
    hold on;
    plot(xRange, polyval(p, xRange), 'r-', 'LineWidth', 2);
    
    % Calculate correlation
    r = corr(anomalyData.temperatureAnomaly, counts - seasonalData.monthlyAverage);
    text(0.05, 0.95, sprintf('r = %.3f', r), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 12, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    grid on;
end

%% ======================== UTILITY FUNCTIONS ========================

function out = num2sepstr(numin, format, sep)
    % NUM2SEPSTR Convert to string with separation at thousands (copied from main script)
    if nargin < 2
        format = ''; % we choose a format below when we know numin is scalar and real
    end
    if nargin < 3
        sep = ',';
    end
    if numel(numin)>1
        out = cell(size(numin));
        for ii = 1:numel(numin)
            out{ii} = num2sepstr(numin(ii), format, sep);
        end
        return
    end
    if ~isreal(numin)
        out = sprintf('%s+%si', ...
            num2sepstr(real(numin), format, sep), ...
            num2sepstr(imag(numin), format, sep));
        return
    end
    autoformat = isempty(format);
    if autoformat
        if isinteger(numin) || mod(round(numin, 4), 1) == 0
            format = '%.0f';
        else
            format = '%.4f'; % 4 digits is the num2str default
        end
    end
    str = sprintf(format, numin);
    if isempty(str)
        error('num2sepstr:invalidFormat', ...
            'Invalid format (sprintf could not use "%s").', format)
    end
    out = regexpi(str, '^(\D*\d{0,3})(\d{3})*(\D\d*)?$', 'tokens', 'once');
    if numel(out)
        out = [out{1}, regexprep(out{2}, '(\d{3})', [sep,'$1']), out{3}];
    else
        out = str;
    end
    if autoformat
        % Trim trailing zeros after the decimal. (By checking above for numbers
        % that look like integers using autoformat, we avoid ever having ONLY
        % zeros after the decimal. There will always be at least one nonzero
        % digit following the decimal.)
        out = regexprep(out, '(\.\d*[1-9])(0*)', '$1');
    end
end