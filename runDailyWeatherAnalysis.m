%% Daily Weather Impact Analysis for Telraam Data
% This script analyzes the impact of weather on cycling during the non-winter
% season (April 1 - November 15) using daily data.
%
% PREREQUISITE: Run runTelraamAnalysis.m first to generate required data
%
% Required variables from workspace:
%   - locationData: processed traffic count data for all locations
%   - weatherData: weather data including temperature, precipitation, windspeed
%   - analysis: analysis parameters including time range and mode
%   - style: plotting style parameters
%
% ADVANTAGES OF DAILY ANALYSIS:
% - More observations (~230 days vs ~33 weeks)
% - Better temporal matching (weather is already at daily resolution)
% - Minimizes seasonal confounding by restricting to April-November
% - Captures day-to-day weather variability

close all; clc

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

fprintf('\n========== DAILY WEATHER IMPACT ANALYSIS ==========\n');
fprintf('Analyzing non-winter weather effects (April 1 - November 15)\n');
fprintf('Data filtered to before 3pm to avoid sunset effects\n');
fprintf('Mode: %s\n\n', analysis.modeDisplayString);

%% Configuration for Daily Analysis
dailyConfig = struct();
dailyConfig.seasonStart = [4, 1];   % April 1
dailyConfig.seasonEnd = [11, 15];   % November 15
dailyConfig.rainyThreshold = 1.0;   % mm/day - above this is "rainy"
dailyConfig.heavyRainThreshold = 10.0;  % mm/day - above this is "heavy rain"
dailyConfig.windyThreshold = 20.0;  % km/h - above this is "windy"
dailyConfig.coldThreshold = 10.0;   % °C - below this is "cold"
dailyConfig.hotThreshold = 25.0;    % °C - above this is "hot"

%% Perform the daily analysis for each location
locationNames = fieldnames(locationData);

for i = 1:length(locationNames)
    locationName = locationNames{i};
    data = locationData.(locationName);
    locationInfo = data.locationInfo;
    
    fprintf('\n\n========================================\n');
    fprintf('Location: %s\n', locationInfo.name);
    fprintf('========================================\n\n');
    
    % Perform comprehensive daily analysis
    performDailyWeatherAnalysis(data, weatherData, analysis, style, dailyConfig);
end

fprintf('\n========== ANALYSIS COMPLETE ==========\n\n');

%% ======================== MAIN ANALYSIS FUNCTION ========================

function performDailyWeatherAnalysis(locationDataStruct, weatherData, analysis, style, config)
    % Perform daily weather impact analysis for a single location
    
    % Step 1: Prepare daily data
    [dailyData, isValid] = prepareDailyData(locationDataStruct, weatherData, analysis, config);
    
    if sum(isValid) < 30
        fprintf('WARNING: Insufficient data for analysis (n=%d days)\n', sum(isValid));
        fprintf('Need at least 30 days of valid data.\n');
        return;
    end
    
    fprintf('Valid data points: %d days\n', sum(isValid));
    fprintf('Date range: %s to %s\n\n', ...
        datestr(min(dailyData.dates(isValid))), ...
        datestr(max(dailyData.dates(isValid))));
    
    % Step 2: Descriptive statistics
    reportDescriptiveStatistics(dailyData, isValid, config);
    
    % Step 3: Categorical comparisons
    reportCategoricalComparisons(dailyData, isValid, config);
    
    % Step 4: Regression analysis
    regressionResults = performRegressionAnalysis(dailyData, isValid);
    
    % Step 4b: Test interaction effects
    interactionResults = testInteractionEffects(dailyData, isValid);
    
    % Step 5: Create visualizations
    createDailyWeatherVisualizations(dailyData, isValid, regressionResults, ...
        locationDataStruct.locationInfo, analysis, style, config);
end

%% ======================== DATA PREPARATION ========================

function [dailyData, isValid] = prepareDailyData(locationDataStruct, weatherData, analysis, config)
    % Prepare daily aggregated data with weather variables and controls
    
    % Calculate daily totals from hourly data
    dailyTraffic = calculateDailyTotals(locationDataStruct, analysis);
    
    % Match with weather data
    [commonDates, ia, ib] = intersect(dailyTraffic.dates, weatherData.dates);
    
    if isempty(commonDates)
        error('No overlapping dates between traffic and weather data');
    end
    
    % Initialize daily data structure
    dailyData = struct();
    dailyData.dates = commonDates;
    dailyData.counts = dailyTraffic.rawCounts(ia);
    dailyData.temperature = weatherData.temperature(ib);
    dailyData.precipitation = weatherData.precipitation(ib);
    dailyData.windspeed = weatherData.windspeed(ib);
    
    % Add temporal controls
    dailyData.isWeekend = isweekend(dailyData.dates);
    dailyData.dayOfWeek = weekday(dailyData.dates);
    dailyData.monthOfYear = month(dailyData.dates);
    dailyData.year = year(dailyData.dates);
    
    % Calculate day within season (for trend control)
    % We'll calculate this after filtering to season, so initialize to zero for now
    dailyData.daysIntoSeason = zeros(size(dailyData.dates));
    
    % Create categorical weather variables
    dailyData.isDry = dailyData.precipitation < config.rainyThreshold;
    dailyData.isRainy = (dailyData.precipitation >= config.rainyThreshold) & ...
                         (dailyData.precipitation < config.heavyRainThreshold);
    dailyData.isHeavyRain = dailyData.precipitation >= config.heavyRainThreshold;
    
    dailyData.isCalm = dailyData.windspeed < config.windyThreshold;
    dailyData.isWindy = dailyData.windspeed >= config.windyThreshold;
    
    dailyData.isCold = dailyData.temperature < config.coldThreshold;
    dailyData.isMild = (dailyData.temperature >= config.coldThreshold) & ...
                        (dailyData.temperature < config.hotThreshold);
    dailyData.isHot = dailyData.temperature >= config.hotThreshold;
    
    % Define "good weather" days
    dailyData.isGoodWeather = dailyData.isDry & dailyData.isCalm & dailyData.isMild;
    
    % Filter to non-winter season (April 1 - November 15)
    uniqueYears = unique(dailyData.year);
    inSeason = false(size(dailyData.dates));
    
    for y = uniqueYears'
        seasonStart = datetime(y, config.seasonStart(1), config.seasonStart(2));
        seasonEnd = datetime(y, config.seasonEnd(1), config.seasonEnd(2), 23, 59, 59);
        yearMask = (dailyData.dates >= seasonStart) & (dailyData.dates <= seasonEnd);
        inSeason = inSeason | yearMask;
        
        % Calculate daysIntoSeason for this year
        if any(yearMask)
            dailyData.daysIntoSeason(yearMask) = days(dailyData.dates(yearMask) - seasonStart);
        end
    end
    
    % Debug output
    fprintf('DEBUG: Total days in dataset: %d\n', length(dailyData.dates));
    fprintf('DEBUG: Days in season (Apr 1 - Nov 15): %d\n', sum(inSeason));
    fprintf('DEBUG: Days with NaN counts: %d\n', sum(isnan(dailyData.counts)));
    fprintf('DEBUG: Days with NaN temperature: %d\n', sum(isnan(dailyData.temperature)));
    fprintf('DEBUG: Days with NaN precipitation: %d\n', sum(isnan(dailyData.precipitation)));
    fprintf('DEBUG: Days with NaN windspeed: %d\n', sum(isnan(dailyData.windspeed)));
    fprintf('DEBUG: Days with negative counts: %d\n', sum(dailyData.counts < 0));
    
    % Identify valid data (in season, no NaNs)
    isValid = inSeason & ...
              ~isnan(dailyData.counts) & ...
              ~isnan(dailyData.temperature) & ...
              ~isnan(dailyData.precipitation) & ...
              ~isnan(dailyData.windspeed) & ...
              (dailyData.counts >= 0);  % Ensure non-negative counts
    
    fprintf('DEBUG: Valid days after all filters: %d\n\n', sum(isValid));
end

function dailyData = calculateDailyTotals(locationDataStruct, analysis)
    % Calculate daily totals from hourly data
    % FILTERS OUT DATA AFTER 3PM to avoid sunset effects
    
    data = locationDataStruct.data;
    
    % Debug: check what we have
    fprintf('DEBUG: Data class: %s\n', class(data));
    if istimetable(data)
        fprintf('DEBUG: Timetable size: %d rows\n', height(data));
        fprintf('DEBUG: Variables: %s\n', strjoin(data.Properties.VariableNames, ', '));
    elseif istable(data)
        fprintf('DEBUG: Table size: %d rows\n', height(data));
        fprintf('DEBUG: Variables: %s\n', strjoin(data.Properties.VariableNames, ', '));
    else
        fprintf('DEBUG: Data is not a table or timetable\n');
    end
    
    % Check if the date/time column exists
    if istimetable(data)
        % For timetable, use the row times
        dateTimes = data.Properties.RowTimes;
    elseif istable(data) && any(strcmp(data.Properties.VariableNames, 'Date and Time (Local)'))
        dateTimes = data.('Date and Time (Local)');
    else
        error('Cannot find date/time information in data');
    end
    
    % FILTER: Keep only data before 3pm (15:00)
    hourOfDay = hour(dateTimes);
    before3pm = hourOfDay < 15;
    
    fprintf('DEBUG: Filtering to before 3pm: keeping %d of %d rows (%.1f%%)\n', ...
        sum(before3pm), length(before3pm), 100*sum(before3pm)/length(before3pm));
    
    % Apply filter
    data = data(before3pm, :);
    dateTimes = dateTimes(before3pm);
    
    % Shift to local date (date at midnight local time)
    localDates = dateshift(dateTimes, 'start', 'day');
    
    % Check if mode column exists
    if ~any(strcmp(data.Properties.VariableNames, analysis.modeString))
        error('Mode "%s" not found in data. Available columns: %s', ...
            analysis.modeString, strjoin(data.Properties.VariableNames, ', '));
    end
    
    % Create temporary table with date and counts
    tempTable = table(localDates, data.(analysis.modeString), ...
        'VariableNames', {'localDate', 'counts'});
    
    % Group by date and sum the mode counts
    groupedData = groupsummary(tempTable, 'localDate', 'sum', 'counts');
    
    dailyData = struct();
    dailyData.dates = groupedData.localDate;
    dailyData.rawCounts = groupedData.sum_counts;
    
    fprintf('DEBUG: Daily data created with %d days (before 3pm only)\n', length(dailyData.dates));
    fprintf('DEBUG: Date range: %s to %s\n', ...
        datestr(min(dailyData.dates)), datestr(max(dailyData.dates)));
end

%% ======================== DESCRIPTIVE STATISTICS ========================

function reportDescriptiveStatistics(dailyData, isValid, config)
    % Report descriptive statistics for the daily data
    
    fprintf('--- DESCRIPTIVE STATISTICS ---\n\n');
    
    % Overall statistics
    fprintf('Daily Counts:\n');
    fprintf('  Mean:   %s\n', num2sepstr(mean(dailyData.counts(isValid))));
    fprintf('  Median: %s\n', num2sepstr(median(dailyData.counts(isValid))));
    fprintf('  Std:    %s\n', num2sepstr(std(dailyData.counts(isValid))));
    fprintf('  Min:    %s\n', num2sepstr(min(dailyData.counts(isValid))));
    fprintf('  Max:    %s\n\n', num2sepstr(max(dailyData.counts(isValid))));
    
    % Weather statistics
    fprintf('Temperature (°C):\n');
    fprintf('  Mean:   %.1f\n', mean(dailyData.temperature(isValid)));
    fprintf('  Median: %.1f\n', median(dailyData.temperature(isValid)));
    fprintf('  Range:  %.1f to %.1f\n\n', ...
        min(dailyData.temperature(isValid)), max(dailyData.temperature(isValid)));
    
    fprintf('Precipitation (mm/day):\n');
    fprintf('  Mean:   %.1f\n', mean(dailyData.precipitation(isValid)));
    fprintf('  Median: %.1f\n', median(dailyData.precipitation(isValid)));
    fprintf('  Max:    %.1f\n\n', max(dailyData.precipitation(isValid)));
    
    fprintf('Wind Speed (km/h):\n');
    fprintf('  Mean:   %.1f\n', mean(dailyData.windspeed(isValid)));
    fprintf('  Median: %.1f\n', median(dailyData.windspeed(isValid)));
    fprintf('  Max:    %.1f\n\n', max(dailyData.windspeed(isValid)));
    
    % Day distribution
    nWeekdays = sum(~dailyData.isWeekend(isValid));
    nWeekends = sum(dailyData.isWeekend(isValid));
    fprintf('Days:\n');
    fprintf('  Weekdays: %d (%.1f%%)\n', nWeekdays, 100*nWeekdays/sum(isValid));
    fprintf('  Weekends: %d (%.1f%%)\n\n', nWeekends, 100*nWeekends/sum(isValid));
    
    % Weather conditions distribution
    nDry = sum(dailyData.isDry(isValid));
    nRainy = sum(dailyData.isRainy(isValid));
    nHeavy = sum(dailyData.isHeavyRain(isValid));
    
    fprintf('Precipitation Categories:\n');
    fprintf('  Dry days (<%s mm):              %d (%.1f%%)\n', ...
        num2sepstr(config.rainyThreshold), nDry, 100*nDry/sum(isValid));
    fprintf('  Rainy days (%s-%s mm):    %d (%.1f%%)\n', ...
        num2sepstr(config.rainyThreshold), num2sepstr(config.heavyRainThreshold), ...
        nRainy, 100*nRainy/sum(isValid));
    fprintf('  Heavy rain days (>%s mm):        %d (%.1f%%)\n\n', ...
        num2sepstr(config.heavyRainThreshold), nHeavy, 100*nHeavy/sum(isValid));
end

%% ======================== CATEGORICAL COMPARISONS ========================

function reportCategoricalComparisons(dailyData, isValid, config)
    % Compare average counts across weather categories
    
    fprintf('--- WEATHER CATEGORY COMPARISONS ---\n\n');
    
    % Separate weekdays and weekends for fair comparison
    weekdayMask = isValid & ~dailyData.isWeekend;
    weekendMask = isValid & dailyData.isWeekend;
    
    if sum(weekdayMask) < 10 || sum(weekendMask) < 10
        fprintf('WARNING: Insufficient data for weekday/weekend comparison\n\n');
        return;
    end
    
    %% Precipitation effect
    fprintf('PRECIPITATION IMPACT (Weekdays only):\n');
    
    dryDays = weekdayMask & dailyData.isDry;
    rainyDays = weekdayMask & dailyData.isRainy;
    heavyRainDays = weekdayMask & dailyData.isHeavyRain;
    
    if sum(dryDays) >= 3
        meanDry = mean(dailyData.counts(dryDays));
        fprintf('  Dry days:         %s (n=%d)\n', num2sepstr(meanDry, '%.0f'), sum(dryDays));
    end
    
    if sum(rainyDays) >= 3
        meanRainy = mean(dailyData.counts(rainyDays));
        fprintf('  Rainy days:       %s (n=%d)\n', num2sepstr(meanRainy, '%.0f'), sum(rainyDays));
        
        if sum(dryDays) >= 3
            pctChange = 100 * (meanRainy - meanDry) / meanDry;
            fprintf('    → %.1f%% vs dry days\n', pctChange);
            
            % Statistical test
            [~, p] = ttest2(dailyData.counts(dryDays), dailyData.counts(rainyDays));
            fprintf('    → t-test p-value: %.4f', p);
            if p < 0.001
                fprintf(' ***\n');
            elseif p < 0.01
                fprintf(' **\n');
            elseif p < 0.05
                fprintf(' *\n');
            else
                fprintf(' (not significant)\n');
            end
        end
    end
    
    if sum(heavyRainDays) >= 3
        meanHeavy = mean(dailyData.counts(heavyRainDays));
        fprintf('  Heavy rain days:  %s (n=%d)\n', num2sepstr(meanHeavy, '%.0f'), sum(heavyRainDays));
        
        if sum(dryDays) >= 3
            pctChange = 100 * (meanHeavy - meanDry) / meanDry;
            fprintf('    → %.1f%% vs dry days\n', pctChange);
        end
    end
    fprintf('\n');
    
    %% Wind effect
    fprintf('WIND IMPACT (Weekdays only):\n');
    
    calmDays = weekdayMask & dailyData.isCalm;
    windyDays = weekdayMask & dailyData.isWindy;
    
    if sum(calmDays) >= 3
        meanCalm = mean(dailyData.counts(calmDays));
        fprintf('  Calm days (<%.0f km/h):  %s (n=%d)\n', ...
            config.windyThreshold, num2sepstr(meanCalm, '%.0f'), sum(calmDays));
    end
    
    if sum(windyDays) >= 3
        meanWindy = mean(dailyData.counts(windyDays));
        fprintf('  Windy days (≥%.0f km/h): %s (n=%d)\n', ...
            config.windyThreshold, num2sepstr(meanWindy, '%.0f'), sum(windyDays));
        
        if sum(calmDays) >= 3
            pctChange = 100 * (meanWindy - meanCalm) / meanCalm;
            fprintf('    → %.1f%% vs calm days\n', pctChange);
            
            [~, p] = ttest2(dailyData.counts(calmDays), dailyData.counts(windyDays));
            fprintf('    → t-test p-value: %.4f', p);
            if p < 0.001
                fprintf(' ***\n');
            elseif p < 0.01
                fprintf(' **\n');
            elseif p < 0.05
                fprintf(' *\n');
            else
                fprintf(' (not significant)\n');
            end
        end
    end
    fprintf('\n');
    
    %% Temperature effect
    fprintf('TEMPERATURE IMPACT (Weekdays only):\n');
    
    coldDays = weekdayMask & dailyData.isCold;
    mildDays = weekdayMask & dailyData.isMild;
    hotDays = weekdayMask & dailyData.isHot;
    
    if sum(coldDays) >= 3
        meanCold = mean(dailyData.counts(coldDays));
        fprintf('  Cold days (<%.0f°C):      %s (n=%d)\n', ...
            config.coldThreshold, num2sepstr(meanCold, '%.0f'), sum(coldDays));
    end
    
    if sum(mildDays) >= 3
        meanMild = mean(dailyData.counts(mildDays));
        fprintf('  Mild days (%.0f-%.0f°C): %s (n=%d)\n', ...
            config.coldThreshold, config.hotThreshold, num2sepstr(meanMild, '%.0f'), sum(mildDays));
    end
    
    if sum(hotDays) >= 3
        meanHot = mean(dailyData.counts(hotDays));
        fprintf('  Hot days (>%.0f°C):       %s (n=%d)\n', ...
            config.hotThreshold, num2sepstr(meanHot, '%.0f'), sum(hotDays));
        
        if sum(mildDays) >= 3
            pctChange = 100 * (meanHot - meanMild) / meanMild;
            fprintf('    → %.1f%% vs mild days\n', pctChange);
        end
    end
    fprintf('\n');
    
    %% Combined "good weather" effect
    fprintf('COMBINED WEATHER CONDITIONS (Weekdays only):\n');
    
    goodWeatherDays = weekdayMask & dailyData.isGoodWeather;
    poorWeatherDays = weekdayMask & ~dailyData.isGoodWeather;
    
    if sum(goodWeatherDays) >= 3 && sum(poorWeatherDays) >= 3
        meanGood = mean(dailyData.counts(goodWeatherDays));
        meanPoor = mean(dailyData.counts(poorWeatherDays));
        
        fprintf('  Good weather days:  %s (n=%d)\n', num2sepstr(meanGood, '%.0f'), sum(goodWeatherDays));
        fprintf('  Poor weather days:  %s (n=%d)\n', num2sepstr(meanPoor, '%.0f'), sum(poorWeatherDays));
        
        pctChange = 100 * (meanGood - meanPoor) / meanPoor;
        fprintf('    → Good weather yields %.1f%% more cycling\n', pctChange);
        
        [~, p] = ttest2(dailyData.counts(goodWeatherDays), dailyData.counts(poorWeatherDays));
        fprintf('    → t-test p-value: %.4f', p);
        if p < 0.001
            fprintf(' ***\n');
        elseif p < 0.01
            fprintf(' **\n');
        elseif p < 0.05
            fprintf(' *\n');
        else
            fprintf(' (not significant)\n');
        end
    end
    fprintf('\n');
    
    %% Weekday vs Weekend comparison
    fprintf('WEEKDAY VS WEEKEND:\n');
    meanWeekday = mean(dailyData.counts(weekdayMask));
    meanWeekend = mean(dailyData.counts(weekendMask));
    
    fprintf('  Weekdays: %s (n=%d)\n', num2sepstr(meanWeekday, '%.0f'), sum(weekdayMask));
    fprintf('  Weekends: %s (n=%d)\n', num2sepstr(meanWeekend, '%.0f'), sum(weekendMask));
    
    pctChange = 100 * (meanWeekend - meanWeekday) / meanWeekday;
    fprintf('    → %.1f%% difference\n\n', pctChange);
end

%% ======================== REGRESSION ANALYSIS ========================

function results = performRegressionAnalysis(dailyData, isValid)
    % Perform multiple regression analysis
    
    fprintf('--- REGRESSION ANALYSIS ---\n\n');
    
    % Extract valid data
    y = dailyData.counts(isValid);
    
    % Build design matrix
    % Predictors: Temperature, Precipitation, Windspeed, Weekend, Time trend
    % Note: No daylight hours needed since we filter to before 3pm
    X = [dailyData.temperature(isValid), ...
         dailyData.precipitation(isValid), ...
         dailyData.windspeed(isValid), ...
         double(dailyData.isWeekend(isValid)), ...
         dailyData.daysIntoSeason(isValid)];
    
    % Check for any remaining NaNs
    validRows = ~any(isnan(X), 2) & ~isnan(y);
    X = X(validRows, :);
    y = y(validRows);
    
    if sum(validRows) < 20
        fprintf('WARNING: Insufficient valid data for regression (n=%d)\n\n', sum(validRows));
        results = struct();
        return;
    end
    
    fprintf('Regression using %d days with complete data.\n', sum(validRows));
    
    % Standardize predictors (for interpretable coefficients)
    [X_std, mu, sigma] = zscore(X);
    X_std = [ones(size(X_std,1),1), X_std];  % Add intercept
    
    % Fit model
    [b, bint, r, rint, stats] = regress(y, X_std);
    
    % Store results
    results = struct();
    results.coefficients = b;
    results.coefficientIntervals = bint;
    results.residuals = r;
    results.rSquared = stats(1);
    results.fStat = stats(2);
    results.pValue = stats(3);
    results.errorVariance = stats(4);
    results.predictorNames = {'Intercept', 'Temperature', 'Precipitation', ...
                              'Windspeed', 'Weekend', 'Time Trend'};
    results.predictedValues = X_std * b;
    results.mu = mu;
    results.sigma = sigma;
    
    % Report results
    fprintf('Multiple Regression Results:\n');
    fprintf('  R² = %.3f\n', results.rSquared);
    fprintf('  F-statistic = %.2f\n', results.fStat);
    fprintf('  p-value = %.6f', results.pValue);
    if results.pValue < 0.001
        fprintf(' ***\n\n');
    elseif results.pValue < 0.01
        fprintf(' **\n\n');
    elseif results.pValue < 0.05
        fprintf(' *\n\n');
    else
        fprintf('\n\n');
    end
    
    fprintf('Standardized Coefficients (effect of 1 SD change):\n');
    for i = 2:length(results.predictorNames)  % Skip intercept
        fprintf('  %-15s: %7.1f  [%.1f, %.1f]', ...
            results.predictorNames{i}, ...
            b(i), bint(i,1), bint(i,2));
        
        % Check if 95% CI excludes zero
        if bint(i,1) > 0 || bint(i,2) < 0
            fprintf(' *');
        end
        fprintf('\n');
    end
    fprintf('\n');
    
    % Interpret in practical terms
    fprintf('Practical Interpretation (holding other factors constant):\n');
    
    % Temperature: effect of 5°C increase
    tempEffect = b(2) * (5 / sigma(1));
    fprintf('  +5°C temperature → %+.0f bikes/day (before 3pm)\n', tempEffect);
    
    % Precipitation: effect of 5mm rain
    precipEffect = b(3) * (5 / sigma(2));
    fprintf('  +5mm precipitation → %+.0f bikes/day (before 3pm)\n', precipEffect);
    
    % Windspeed: effect of 10 km/h increase
    windEffect = b(4) * (10 / sigma(3));
    fprintf('  +10 km/h wind → %+.0f bikes/day (before 3pm)\n', windEffect);
    
    % Weekend effect
    weekendEffect = b(5) * (1 / sigma(4));
    fprintf('  Weekend vs weekday → %+.0f bikes/day (before 3pm)\n', weekendEffect);
    
    fprintf('\n');
    
    % Model diagnostics
    fprintf('Model Diagnostics:\n');
    fprintf('  Mean residual: %.1f (should be ≈0)\n', mean(r));
    fprintf('  Residual std: %.1f\n', std(r));
    
    % Check for outliers (residuals > 3 SD)
    outlierThreshold = 3 * std(r);
    nOutliers = sum(abs(r) > outlierThreshold);
    fprintf('  Outliers (|residual| > 3 SD): %d (%.1f%%)\n', ...
        nOutliers, 100*nOutliers/length(r));
    
    fprintf('\n');
end

%% ======================== INTERACTION EFFECTS TESTING ========================

function results = testInteractionEffects(dailyData, isValid)
    % Test all pairwise interaction effects among weather variables
    
    fprintf('\n--- INTERACTION EFFECTS ANALYSIS ---\n\n');
    
    % Extract valid data
    y = dailyData.counts(isValid);
    
    % Base predictors
    temp = dailyData.temperature(isValid);
    precip = dailyData.precipitation(isValid);
    wind = dailyData.windspeed(isValid);
    weekend = double(dailyData.isWeekend(isValid));
    timetrend = dailyData.daysIntoSeason(isValid);
    
    % Check for valid data
    validRows = ~any(isnan([temp, precip, wind, weekend, timetrend]), 2) & ~isnan(y);
    temp = temp(validRows);
    precip = precip(validRows);
    wind = wind(validRows);
    weekend = weekend(validRows);
    timetrend = timetrend(validRows);
    y = y(validRows);
    
    if sum(validRows) < 20
        fprintf('WARNING: Insufficient valid data for interaction testing (n=%d)\n\n', sum(validRows));
        results = struct();
        return;
    end
    
    % Standardize continuous predictors for interpretability
    temp_std = zscore(temp);
    precip_std = zscore(precip);
    wind_std = zscore(wind);
    timetrend_std = zscore(timetrend);
    
    % Base model (no interactions)
    X_base = [ones(length(y),1), temp_std, precip_std, wind_std, weekend, timetrend_std];
    [b_base, ~, ~, ~, stats_base] = regress(y, X_base);
    R2_base = stats_base(1);
    
    fprintf('Base Model (no interactions): R² = %.3f\n\n', R2_base);
    
    % Test each pairwise interaction
    interactions = {
        {'Temp × Precip', temp_std .* precip_std}
        {'Temp × Wind', temp_std .* wind_std}
        {'Temp × Weekend', temp_std .* weekend}
        {'Temp × Time', temp_std .* timetrend_std}
        {'Precip × Wind', precip_std .* wind_std}
        {'Precip × Weekend', precip_std .* weekend}
        {'Precip × Time', precip_std .* timetrend_std}
        {'Wind × Weekend', wind_std .* weekend}
        {'Wind × Time', wind_std .* timetrend_std}
        {'Weekend × Time', weekend .* timetrend_std}
    };
    
    fprintf('Testing individual interaction terms:\n');
    fprintf('%-20s %8s %8s %10s %12s\n', 'Interaction', 'Coef', 'R²', 'ΔR²', 'p-value');
    fprintf('%s\n', repmat('-', 1, 62));
    
    results = struct();
    results.interactions = cell(length(interactions), 1);
    results.coefficients = zeros(length(interactions), 1);
    results.rSquared = zeros(length(interactions), 1);
    results.deltaR2 = zeros(length(interactions), 1);
    results.pValues = zeros(length(interactions), 1);
    
    for i = 1:length(interactions)
        name = interactions{i}{1};
        term = interactions{i}{2};
        
        % Add interaction term to base model
        X_int = [X_base, term];
        [b_int, ~, ~, ~, stats_int] = regress(y, X_int);
        
        R2_int = stats_int(1);
        deltaR2 = R2_int - R2_base;
        
        % F-test for the interaction term
        % F = (SSR_full - SSR_reduced) / (df_full - df_reduced) / MSE_full
        n = length(y);
        k_base = size(X_base, 2) - 1;  % degrees of freedom (excluding intercept)
        k_int = size(X_int, 2) - 1;
        
        SSR_base = stats_base(4) * (n - k_base - 1);  % Residual sum of squares
        SSR_int = stats_int(4) * (n - k_int - 1);
        
        F_stat = ((SSR_base - SSR_int) / (k_int - k_base)) / stats_int(4);
        p_value = 1 - fcdf(F_stat, k_int - k_base, n - k_int - 1);
        
        % Store results
        results.interactions{i} = name;
        results.coefficients(i) = b_int(end);  % Coefficient of interaction term
        results.rSquared(i) = R2_int;
        results.deltaR2(i) = deltaR2;
        results.pValues(i) = p_value;
        
        % Print results
        fprintf('%-20s %8.1f %8.3f %10.4f', name, b_int(end), R2_int, deltaR2);
        if p_value < 0.001
            fprintf(' %12.4f ***\n', p_value);
        elseif p_value < 0.01
            fprintf(' %12.4f **\n', p_value);
        elseif p_value < 0.05
            fprintf(' %12.4f *\n', p_value);
        else
            fprintf(' %12.4f\n', p_value);
        end
    end
    
    fprintf('\n');
    
    % Test full model with all significant interactions
    significantIdx = results.pValues < 0.05;
    if any(significantIdx)
        fprintf('Full model with significant interactions (p < 0.05):\n');
        
        % Build model with all significant interactions
        X_full = X_base;
        sigNames = {};
        for i = 1:length(interactions)
            if significantIdx(i)
                X_full = [X_full, interactions{i}{2}];
                sigNames{end+1} = interactions{i}{1};
            end
        end
        
        [b_full, ~, ~, ~, stats_full] = regress(y, X_full);
        R2_full = stats_full(1);
        F_full = stats_full(2);
        p_full = stats_full(3);
        
        fprintf('  Including: %s\n', strjoin(sigNames, ', '));
        fprintf('  R² = %.3f (vs %.3f base)\n', R2_full, R2_base);
        fprintf('  ΔR² = %.3f\n', R2_full - R2_base);
        fprintf('  F-statistic = %.2f, p = %.6f', F_full, p_full);
        if p_full < 0.001
            fprintf(' ***\n');
        elseif p_full < 0.01
            fprintf(' **\n');
        elseif p_full < 0.05
            fprintf(' *\n');
        else
            fprintf('\n');
        end
        
        results.fullModel = struct();
        results.fullModel.coefficients = b_full;
        results.fullModel.rSquared = R2_full;
        results.fullModel.deltaR2 = R2_full - R2_base;
    else
        fprintf('No significant interactions found (p < 0.05).\n');
        results.fullModel = struct();
    end
    
    fprintf('\n');
end

%% ======================== VISUALIZATION ========================

function createDailyWeatherVisualizations(dailyData, isValid, regressionResults, ...
    locationInfo, analysis, style, config)
    % Create comprehensive visualization of daily weather effects
    
    % Extract valid data for plotting
    dates = dailyData.dates(isValid);
    counts = dailyData.counts(isValid);
    temp = dailyData.temperature(isValid);
    precip = dailyData.precipitation(isValid);
    wind = dailyData.windspeed(isValid);
    isWeekend = dailyData.isWeekend(isValid);
    
    %% Figure 1: Time series overview
    fig1 = figure('Position', [100 100 1400 900]);
    sgtitle(sprintf('%s - %s\nDaily Weather Impact Analysis', ...
        locationInfo.name, analysis.modeDisplayString), ...
        'FontSize', style.titleFontSize, 'FontWeight', 'bold');
    
    % Subplot 1: Daily counts with weather overlay
    subplot(4,1,1);
    plot(dates, counts, 'k-', 'LineWidth', 1);
    hold on;
    
    % Highlight rainy days
    rainyDays = precip >= config.rainyThreshold;
    scatter(dates(rainyDays), counts(rainyDays), 50, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    
    ylabel('Daily Counts', 'FontSize', style.labelFontSize);
    legend('All Days', 'Rainy Days (≥1mm)', 'Location', 'northwest', ...
        'FontSize', style.legendFontSize);
    title('Daily Counts with Rain Indicators', 'FontSize', style.titleFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    % Subplot 2: Temperature
    subplot(4,1,2);
    plot(dates, temp, 'r-', 'LineWidth', 1.5);
    ylabel('Temperature (°C)', 'FontSize', style.labelFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    % Subplot 3: Precipitation
    subplot(4,1,3);
    bar(dates, precip, 'FaceColor', [0.3 0.6 1], 'EdgeColor', 'none');
    ylabel('Precipitation (mm)', 'FontSize', style.labelFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    % Subplot 4: Wind speed
    subplot(4,1,4);
    plot(dates, wind, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    ylabel('Wind Speed (km/h)', 'FontSize', style.labelFontSize);
    xlabel('Date', 'FontSize', style.labelFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    %% Figure 2: Weather relationships
    fig2 = figure('Position', [150 150 1400 900]);
    sgtitle(sprintf('%s - %s\nWeather Relationships', ...
        locationInfo.name, analysis.modeDisplayString), ...
        'FontSize', style.titleFontSize, 'FontWeight', 'bold');
    
    % Subplot 1: Temperature vs Counts
    subplot(2,2,1);
    
    % Separate weekdays and weekends
    weekdayMask = ~isWeekend;
    weekendMask = isWeekend;
    
    scatter(temp(weekdayMask), counts(weekdayMask), 60, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(temp(weekendMask), counts(weekendMask), 60, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    
    % Add regression line (weekdays only)
    if sum(weekdayMask) >= 10
        p = polyfit(temp(weekdayMask), counts(weekdayMask), 1);
        tempRange = linspace(min(temp), max(temp), 100);
        plot(tempRange, polyval(p, tempRange), 'b-', 'LineWidth', 2);
        
        % Calculate correlation
        r = corr(temp(weekdayMask), counts(weekdayMask));
        text(0.05, 0.95, sprintf('Weekday r = %.3f', r), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    xlabel('Temperature (°C)', 'FontSize', style.labelFontSize);
    ylabel('Daily Counts', 'FontSize', style.labelFontSize);
    legend('Weekdays', 'Weekends', 'Weekday Trend', 'Location', 'northwest', ...
        'FontSize', style.legendFontSize);
    title('Temperature Effect', 'FontSize', style.titleFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    % Subplot 2: Precipitation vs Counts
    subplot(2,2,2);
    
    scatter(precip(weekdayMask), counts(weekdayMask), 60, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(precip(weekendMask), counts(weekendMask), 60, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    
    xlabel('Precipitation (mm)', 'FontSize', style.labelFontSize);
    ylabel('Daily Counts', 'FontSize', style.labelFontSize);
    legend('Weekdays', 'Weekends', 'Location', 'northwest', ...
        'FontSize', style.legendFontSize);
    title('Precipitation Effect', 'FontSize', style.titleFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    % Subplot 3: Wind vs Counts
    subplot(2,2,3);
    
    scatter(wind(weekdayMask), counts(weekdayMask), 60, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(wind(weekendMask), counts(weekendMask), 60, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.5);
    
    xlabel('Wind Speed (km/h)', 'FontSize', style.labelFontSize);
    ylabel('Daily Counts', 'FontSize', style.labelFontSize);
    legend('Weekdays', 'Weekends', 'Location', 'northwest', ...
        'FontSize', style.legendFontSize);
    title('Wind Effect', 'FontSize', style.titleFontSize);
    grid on;
    set(gca, 'FontSize', style.axisFontSize);
    
    % Subplot 4: Predicted vs Actual
    if ~isempty(regressionResults) && isfield(regressionResults, 'predictedValues')
        subplot(2,2,4);
        
        % Need to match predicted values to valid data
        validRegression = ~any(isnan([temp, precip, wind, double(isWeekend), ...
                                      dailyData.daysIntoSeason(isValid)]), 2);
        
        actualCounts = counts(validRegression);
        predictedCounts = regressionResults.predictedValues;
        
        scatter(predictedCounts, actualCounts, 60, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.5);
        hold on;
        
        % Add perfect prediction line
        minVal = min([actualCounts; predictedCounts]);
        maxVal = max([actualCounts; predictedCounts]);
        plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 2);
        
        xlabel('Predicted Counts', 'FontSize', style.labelFontSize);
        ylabel('Actual Counts', 'FontSize', style.labelFontSize);
        title(sprintf('Model Fit (R² = %.3f)', regressionResults.rSquared), ...
            'FontSize', style.titleFontSize);
        grid on;
        set(gca, 'FontSize', style.axisFontSize);
        axis equal;
        xlim([minVal maxVal]);
        ylim([minVal maxVal]);
    end
    
    %% Figure 3: Categorical comparisons
    fig3 = figure('Position', [200 200 1400 600]);
    sgtitle(sprintf('%s - %s\nWeather Category Comparisons', ...
        locationInfo.name, analysis.modeDisplayString), ...
        'FontSize', style.titleFontSize, 'FontWeight', 'bold');
    
    % Subplot 1: Precipitation categories
    subplot(1,3,1);
    
    dryDays = precip < config.rainyThreshold;
    rainyDays = (precip >= config.rainyThreshold) & (precip < config.heavyRainThreshold);
    heavyDays = precip >= config.heavyRainThreshold;
    
    % Collect data for boxplot
    groupData = [];
    groupLabels = [];
    
    if sum(dryDays & weekdayMask) >= 3
        dryData = counts(dryDays & weekdayMask);
        groupData = [groupData; dryData(:)];  % Force column vector
        groupLabels = [groupLabels; repmat({'Dry (WD)'}, length(dryData), 1)];
    end
    if sum(rainyDays & weekdayMask) >= 3
        rainyData = counts(rainyDays & weekdayMask);
        groupData = [groupData; rainyData(:)];  % Force column vector
        groupLabels = [groupLabels; repmat({'Rainy (WD)'}, length(rainyData), 1)];
    end
    if sum(heavyDays & weekdayMask) >= 3
        heavyData = counts(heavyDays & weekdayMask);
        groupData = [groupData; heavyData(:)];  % Force column vector
        groupLabels = [groupLabels; repmat({'Heavy (WD)'}, length(heavyData), 1)];
    end
    
    if ~isempty(groupData)
        boxplot(groupData, groupLabels);
        ylabel('Daily Counts', 'FontSize', style.labelFontSize);
        title('Precipitation Impact', 'FontSize', style.titleFontSize);
        set(gca, 'FontSize', style.axisFontSize);
        grid on;
    end
    
    % Subplot 2: Wind categories
    subplot(1,3,2);
    
    calmDays = wind < config.windyThreshold;
    windyDays = wind >= config.windyThreshold;
    
    groupData = [];
    groupLabels = [];
    
    if sum(calmDays & weekdayMask) >= 3
        calmData = counts(calmDays & weekdayMask);
        groupData = [groupData; calmData(:)];
        groupLabels = [groupLabels; repmat({'Calm (WD)'}, length(calmData), 1)];
    end
    if sum(windyDays & weekdayMask) >= 3
        windyData = counts(windyDays & weekdayMask);
        groupData = [groupData; windyData(:)];
        groupLabels = [groupLabels; repmat({'Windy (WD)'}, length(windyData), 1)];
    end
    
    if ~isempty(groupData)
        boxplot(groupData, groupLabels);
        ylabel('Daily Counts', 'FontSize', style.labelFontSize);
        title('Wind Impact', 'FontSize', style.titleFontSize);
        set(gca, 'FontSize', style.axisFontSize);
        grid on;
    end
    
    % Subplot 3: Combined weather quality
    subplot(1,3,3);
    
    goodWeather = dailyData.isGoodWeather(isValid);
    
    groupData = [];
    groupLabels = [];
    
    if sum(goodWeather & weekdayMask) >= 3
        goodWDData = counts(goodWeather & weekdayMask);
        groupData = [groupData; goodWDData(:)];
        groupLabels = [groupLabels; repmat({'Good (WD)'}, length(goodWDData), 1)];
    end
    if sum(~goodWeather & weekdayMask) >= 3
        poorWDData = counts(~goodWeather & weekdayMask);
        groupData = [groupData; poorWDData(:)];
        groupLabels = [groupLabels; repmat({'Poor (WD)'}, length(poorWDData), 1)];
    end
    if sum(goodWeather & weekendMask) >= 3
        goodWEData = counts(goodWeather & weekendMask);
        groupData = [groupData; goodWEData(:)];
        groupLabels = [groupLabels; repmat({'Good (WE)'}, length(goodWEData), 1)];
    end
    if sum(~goodWeather & weekendMask) >= 3
        poorWEData = counts(~goodWeather & weekendMask);
        groupData = [groupData; poorWEData(:)];
        groupLabels = [groupLabels; repmat({'Poor (WE)'}, length(poorWEData), 1)];
    end
    
    if ~isempty(groupData)
        boxplot(groupData, groupLabels);
        ylabel('Daily Counts', 'FontSize', style.labelFontSize);
        title('Overall Weather Quality', 'FontSize', style.titleFontSize);
        set(gca, 'FontSize', style.axisFontSize);
        grid on;
    end
end

%% ======================== UTILITY FUNCTIONS ========================

function out = num2sepstr(numin, format, sep)
    % NUM2SEPSTR Convert to string with separation at thousands
    % Copied from main script for independence
    if nargin < 2
        format = '';
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
            format = '%.4f';
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
        out = regexprep(out, '(\.\d*[1-9])(0*)', '$1');
    end
end