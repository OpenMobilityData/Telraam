%% Modular Telraam Analysis - Single Script Version
% This version can be run as a single script while maintaining modularity

clear all; close all; clc;

%% Configuration Setup

westernSegmentName = 'rue de Terrebonne @ King Edward';
easternSegmentName = 'rue de Terrebonne @ Draper';

% Locations to analyze
locations = {
    struct('name', easternSegmentName, ...
           'fileStem2024', 'telraam-raw-data-2024East60', ...
           'fileStem2025', 'telraam-raw-data-2025East60', ...
           'plotColor', [0 0 1]);  % Blue
    struct('name', westernSegmentName, ...
           'fileStem2024', 'telraam-raw-data-2024West60', ...
           'fileStem2025', 'telraam-raw-data-2025West60', ...
           'plotColor', [0 0 0]);  % Black
};

% Analysis parameters

modeString = 'Bike Total'; modeDisplayString = 'Bike Counts';
%modeString = 'Pedestrian Total'; modeDisplayString = 'Pedestrian Counts';
%modeString = 'Car Total'; modeDisplayString = 'Car Counts';
%modeString = 'Large vehicle Total'; modeDisplayString = 'Heavy Truck Counts';

analysis = struct( ...
    'startTime', datetime(2024,08,01,00,00,01), ...
    'endTime', datetime(2025,06,01,23,59,59), ...
    'modeString', modeString, ...
    'modeDisplayString', modeDisplayString, ...
    'uptimeThreshold', 0.0, ...
    'maxUptimeCorrection', 1.0, ...
    'truncationCutoffTime', timeofday(datetime('today')+hours(15)), ...
    'daylightCorrectionRatioWD', 1.65, ...
    'daylightCorrectionRatioWE', 1.37, ...
    'includePartialMonths', false ...
);

% Plot configuration - easy to turn elements on/off
plots = struct( ...
    'showRawCounts', true, ...
    'showAdjustedCounts', false, ...  % Turn off for first example
    'showTruncatedCounts', false, ...
    'showWeather', true, ...
    'showPrecipitationBubbles', true, ...  % Turn off for first example
    'plotTypes', {{'daily','weekly','monthly'}}, ...  % Only daily for now
    'combinedPlots', true ...
);

% Plotting style parameters
style = struct( ...
    'plotLineWidth', 10.0, ...
    'axisFontSize', 16.0, ...
    'labelFontSize', 20.0, ...
    'titleFontSize', 24.0, ...
    'legendFontSize', 16.0, ...
    'axisBackgroundColor', 0.8.*[1 1 1], ...
    'legendBackgroundAlpha', 0.2 ...
);

% Multi-modal analysis parameters
multiModal = struct( ...
    'enabled', true, ...
    'location', easternSegmentName, ...  % Which location to analyze
    'modes', {{'Bike Total', 'Pedestrian Total', 'Car Total'}}, ...
    'modeDisplayNames', {{'Bike Counts', 'Pedestrian Counts', 'Car Counts'}}, ...
    'modeColors', {{[0 0 1], [0 0.8 0], [1 0 0], [0 0.8 0.8]}}, ...  % Note the double braces
    'plotWeather', true ...
);


%% Load and Process Data for All Locations
locationData = struct();

for i = 1:length(locations)
    location = locations{i};
    fprintf('Loading data for %s...\n', location.name);
    
    % Load raw data
    rawData = loadSingleLocationData(location, analysis);
    
    % Process data
    processedData = processTelraamData(rawData, analysis);
    
    % Store in structure
    fieldName = matlab.lang.makeValidName(location.name);
    locationData.(fieldName).data = processedData;
    locationData.(fieldName).locationInfo = location;
end

%% Get Weather Data (once for all locations) - with caching
% Get unique days across all locations
allDates = [];
locationNames = fieldnames(locationData);
for i = 1:length(locationNames)
    dates = locationData.(locationNames{i}).data.('Date and Time (Local)');
    allDates = [allDates; dates];
end

uniqueDays = unique(dateshift(allDates, 'start', 'day'));
dailyNoonTimes = uniqueDays + hours(12);

% Create cache filename based on date range
startDateStr = datestr(uniqueDays(1), 'yyyy-mm-dd');
endDateStr = datestr(uniqueDays(end), 'yyyy-mm-dd');
cacheFilename = sprintf('weatherCache_%s_to_%s.mat', startDateStr, endDateStr);

% Check if cached data exists
if exist(cacheFilename, 'file')
    fprintf('Loading cached weather data from %s...\n', cacheFilename);
    load(cacheFilename, 'weatherData');
    fprintf('Loaded weather data for %d days from cache.\n', length(weatherData.dates));
else
    % Get weather data from API
    fprintf('Getting weather data for %d days from API...\n', length(uniqueDays));
    [precipitationData, temperatureData, sunriseData, sunsetData, sunhoursData, snowData, windspeedData, feelslikeData] = ...
        getWeatherstackData('Montreal', dailyNoonTimes);
    
    % Store weather data correctly - extract the actual vectors
    weatherData = struct();
    weatherData.dates = uniqueDays;
    weatherData.precipitation = precipitationData;
    weatherData.temperature = temperatureData;
    weatherData.feelslike = feelslikeData;
    weatherData.windspeed = windspeedData;
    weatherData.sunrise = sunriseData;
    weatherData.sunset = sunsetData;
    weatherData.sunhours = sunhoursData;
    weatherData.snow = snowData;
    
    % Save to cache
    fprintf('Saving weather data to cache: %s\n', cacheFilename);
    save(cacheFilename, 'weatherData');
end

%% Generate Combined Daily Plot
plotCombinedDaily(locationData, weatherData, analysis, plots, style);

%% Generate Combined Weekly Plot
plotCombinedWeekly(locationData, weatherData, analysis, plots, style);

%% Generate Combined Monthly Plot
plotCombinedMonthly(locationData, weatherData, analysis, plots, style);

%% Generate Multi-Modal Plots (if enabled)
if multiModal.enabled
    plotMultiModalDaily(locationData, weatherData, analysis, plots, style, multiModal);
    plotMultiModalWeekly(locationData, weatherData, analysis, plots, style, multiModal);
    plotMultiModalMonthly(locationData, weatherData, analysis, plots, style, multiModal);
end

%% ======================== FUNCTIONS ========================

function inputTable = loadSingleLocationData(location, analysis)
    % Load 2024 data
    inputTable2024 = loadYearData(location.fileStem2024, 2024);
    
    % Load 2025 data
    inputTable2025 = loadYearData(location.fileStem2025, 2025);
    
    % Combine and convert to timetable
    columnsToKeep = {'Uptime','Date and Time (Local)','Bike Total','Pedestrian Total',...
                     'Night Total','Speed V85 km/h','Car Total','Large vehicle Total'};
    inputTable2024 = inputTable2024(:,columnsToKeep);
    inputTable2025 = inputTable2025(:,columnsToKeep);
    
    combinedTable = vertcat(inputTable2024, inputTable2025);
    inputTable = table2timetable(combinedTable);
    
    % Apply date range filter
    inputTable = inputTable(inputTable.('Date and Time (Local)') >= analysis.startTime & ...
                           inputTable.('Date and Time (Local)') <= analysis.endTime, :);
end

function yearTable = loadYearData(fileStem, targetYear)
    excelFileName = [fileStem '.xlsx'];
    
    % Read table with preserved variable names to avoid the warning
    opts = detectImportOptions(excelFileName);
    opts.VariableNamingRule = 'preserve';
    yearTable = readtable(excelFileName, opts);
    
    % Uptime is no longer numeric - should remove
    if ismember('Uptime', yearTable.Properties.VariableNames)
        yearTable.Uptime = ones(size(yearTable.Uptime));
    else
        yearTable.Uptime = ones(height(yearTable), 1);
    end

    % Convert dates and filter by year
    yearTable.('Date and Time (Local)') = datetime(char(yearTable.('Date and Time (Local)')),'Format','yy-MM-dd HH:mm');
    yearTable = yearTable(year(yearTable.('Date and Time (Local)')) == targetYear, :);
end

function processedTable = processTelraamData(inputTable, analysis)
    % Add descriptive columns
    processedTable = addTemporalColumns(inputTable);
    
    % Apply uptime corrections
    processedTable = applyUptimeCorrections(processedTable, analysis);
    
    % Apply daylight corrections
    processedTable = applyDaylightCorrections(processedTable, analysis);
end

function inputTable = addTemporalColumns(inputTable)
    weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
    
    inputTable.dayOfWeek = string(day(inputTable.('Date and Time (Local)'),'name'));
    inputTable.dayOfWeekCat = categorical(inputTable.dayOfWeek);
    inputTable.weekOfYear = week(inputTable.('Date and Time (Local)'),'iso-weekofyear');
    inputTable.isWeekday = ismember(inputTable.dayOfWeekCat, weekdays);
    inputTable.weekStartDateTimes = dateshift(dateshift(inputTable.('Date and Time (Local)'),'dayofweek','Monday','previous'),'start','day');
    inputTable.Daylight = ~((inputTable.('Night Total') > 0) | (isnan(inputTable.('Speed V85 km/h'))));
    inputTable.DaylightUptime = inputTable.Daylight .* inputTable.Uptime;
    
    % Fix week numbering
    inputTable = fixWeekNumbering(inputTable);
end

function inputTable = fixWeekNumbering(inputTable)
    inputTable.weekOfYear((month(inputTable.('Date and Time (Local)'))==12) & (inputTable.weekOfYear==1)) = 53;
    inputTable.weekOfYear((month(inputTable.('Date and Time (Local)'))==1) & (inputTable.weekOfYear==1)) = 53;
    inputTable.yearOfMondayInWeek = year(inputTable.weekStartDateTimes);
    januaryIndicesToChange = (month(inputTable.weekStartDateTimes)==1) & (inputTable.weekOfYear==53);
    inputTable.yearOfMondayInWeek(januaryIndicesToChange) = inputTable.yearOfMondayInWeek(januaryIndicesToChange) - 1;
    inputTable.yearWeekKey = inputTable.yearOfMondayInWeek + inputTable.weekOfYear./100;
end

function inputTable = applyUptimeCorrections(inputTable, analysis)
    if ismember(analysis.modeString, inputTable.Properties.VariableNames)
        uptimeCorrection = 1./inputTable.Uptime;
        uptimeCorrection(uptimeCorrection > analysis.maxUptimeCorrection) = analysis.maxUptimeCorrection;
        inputTable.AdjustedCountsUptime = inputTable.(analysis.modeString) .* uptimeCorrection;
    else
        warning('Mode string %s not found in data', analysis.modeString);
        inputTable.AdjustedCountsUptime = zeros(height(inputTable), 1);
    end
end

function inputTable = applyDaylightCorrections(inputTable, analysis)
    % Split into weekday/weekend
    inputTableWD = inputTable(inputTable.isWeekday,:);
    inputTableWE = inputTable(~inputTable.isWeekday,:);
    
    % Apply daylight corrections
    inputTableWD.AdjustedCountsUptimeDaylight = inputTableWD.AdjustedCountsUptime .* analysis.daylightCorrectionRatioWD;
    inputTableWE.AdjustedCountsUptimeDaylight = inputTableWE.AdjustedCountsUptime .* analysis.daylightCorrectionRatioWE;
    
    % Recombine
    inputTable = sortrows([inputTableWD; inputTableWE], 'Date and Time (Local)');
end

function plotCombinedDaily(locationData, weatherData, analysis, plots, style)
    figure('Position', [408 126 1132 921]);
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    if plots.showWeather
        yyaxis right
        weatherHandles = plotWeatherData(weatherData, plots, style);
        %plotHandles = [plotHandles, weatherHandles];
    end
    
    % Plot traffic data on left axis
    yyaxis left
    hold on
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate daily totals
        dailyData = calculateDailyTotals(data, analysis);
        
        % Plot raw counts if enabled
        if plots.showRawCounts
            h1 = plot(dailyData.dates, dailyData.rawCounts, '-', ...
                'LineWidth', style.plotLineWidth * 0.5, ...
                'Color', locationInfo.plotColor, ...
                'DisplayName', sprintf('%s Raw ( Min = %s ; Max = %s ; Total = %s )', ...
                    locationInfo.name, ...
                    num2sepstr(min(dailyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(dailyData.rawCounts), '%.0f'), ...
                    num2sepstr(sum(dailyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
        
        % Plot adjusted counts if enabled
        if plots.showAdjustedCounts
            h2 = plot(dailyData.dates, dailyData.adjustedCounts, '--', ...
                'LineWidth', style.plotLineWidth * 0.3, ...
                'Color', locationInfo.plotColor * 0.7, ...
                'DisplayName', sprintf('%s Adjusted ( Min = %s ; Max = %s ; Total = %s )', ...
                    locationInfo.name, ...
                    num2sepstr(min(dailyData.adjustedCounts), '%.0f'), ...
                    num2sepstr(max(dailyData.adjustedCounts), '%.0f'), ...
                    num2sepstr(sum(dailyData.adjustedCounts), '%.0f')));
            plotHandles = [plotHandles, h2];
        end
    end
    
    % Format plot
    formatCombinedPlot(analysis, plots, style, weatherData);
    
    if plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end

    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'north', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function dailyData = calculateDailyTotals(locationDataStruct, analysis)
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Create daily grouping
    data.DayOnly = dateshift(data.('Date and Time (Local)'), 'start', 'day');
    
    % Calculate raw daily totals
    groupedData = groupsummary(data, 'DayOnly', 'sum', analysis.modeString);
    
    % Calculate adjusted daily totals (for times up to cutoff)
    truncatedData = data(timeofday(data.('Date and Time (Local)')) <= analysis.truncationCutoffTime, :);
    adjustedGrouped = groupsummary(truncatedData, 'DayOnly', 'sum', 'AdjustedCountsUptimeDaylight');
    
    % Also get daylight data counts per day to identify days with no daylight data
    daylightGrouped = groupsummary(data, 'DayOnly', 'sum', 'Daylight');
    
    % Combine results
    dailyData = struct();
    dailyData.dates = groupedData.DayOnly;
    
    % Use the original column name with groupsummary's 'sum_' prefix
    sumColumnName = ['sum_' analysis.modeString];
    dailyData.rawCounts = groupedData.(sumColumnName);
    
    % Match adjusted counts to raw count dates
    [~, ia, ib] = intersect(groupedData.DayOnly, adjustedGrouped.DayOnly);
    adjustedCounts = nan(size(dailyData.rawCounts));
    adjustedCounts(ia) = adjustedGrouped.sum_AdjustedCountsUptimeDaylight(ib);
    
    % Ensure adjusted is at least as large as raw
    dailyData.adjustedCounts = max(dailyData.rawCounts, adjustedCounts);
    
    % Filter out days with no daylight data (sum_Daylight = 0)
    [~, ic, id] = intersect(groupedData.DayOnly, daylightGrouped.DayOnly);
    daylightCounts = zeros(size(dailyData.rawCounts));
    daylightCounts(ic) = daylightGrouped.sum_Daylight(id);
    
    % Keep only days that have some daylight data
    validDays = daylightCounts > 0;
    dailyData.dates = dailyData.dates(validDays);
    dailyData.rawCounts = dailyData.rawCounts(validDays);
    dailyData.adjustedCounts = dailyData.adjustedCounts(validDays);
end

function weatherHandles = plotWeatherData(weatherData, plots, style)
    % Plot temperature line - use the same approach as original script
    h1 = plot(weatherData.dates, weatherData.temperature, '-', ...
        'LineWidth', style.plotLineWidth * 0.5, ...
        'Color', [0 0.4471 0.7412 0.3], ...
        'DisplayName', 'Temperature (°C)');
    
    weatherHandles = h1;
    
    % Add precipitation bubbles if enabled
    if plots.showPrecipitationBubbles
        hold on
        h2 = bubblechart(weatherData.dates, weatherData.temperature, weatherData.precipitation, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.3, ...
            'DisplayName', 'Precipitation (bubble size)');
        weatherHandles = [weatherHandles, h2];
        hold off
    end
    
    ylabel('Temperature (°C)', 'FontSize', style.labelFontSize, ...
        'Color', [1 0 0 0.5] + 0.5.*[0 1 1 0], 'FontWeight', 'bold');
end

function formatCombinedPlot(analysis, plots, style, weatherData)
    % Left axis formatting
    yyaxis left
    ylabel(['Total ' analysis.modeDisplayString ' per Day'], ...
        'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    if plots.showWeather
        ax.YAxis(2).Color = [0 0.4471 0.7412];
    end
    
    % Title and formatting
    title(['Daily ' analysis.modeDisplayString], ...
        'FontSize', style.titleFontSize);
    
    xlabel('Date', 'FontSize', style.labelFontSize);
    set(gca, 'Color', style.axisBackgroundColor);
    grid on;
    xtickangle(45);
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    if ~isempty(weatherData)
        xlim([weatherData.dates(1) weatherData.dates(end)]);
    end
end

% Keep your existing utility functions
function [precipitationValues, temperatureValues, sunriseValues, sunsetValues, sunhoursValues, snowValues, windspeedValues, feelslikeValues] = getWeatherstackData(location, dates)
    % Define the base URL for Weatherstack API
    baseURL = 'http://api.weatherstack.com/historical';
    
    % Replace with your actual Weatherstack access key (if needed)
    accessKey = 'fe2d67122ba14cc9e0b2c931f6105b4b'; % Leave empty if access key is not required
    
    % Iterate over each date in the range
    numDates = length(dates);
    precipitationValues = [];
    windspeedValues = [];
    feelslikeValues = [];
    snowValues = [];
    temperatureValues = [];
    sunhoursValues = [];
    
    for ix=1:numDates
        currentDate = dates(ix);
        
        % Convert the current date to string format
        dateStr = datestr(currentDate, 'yyyy-mm-dd');
        
        % Construct the API request URL
        requestURL = sprintf('%s?access_key=%s&query=%s&historical_date=%s&hourly=1&interval=24&units=m', baseURL, accessKey, location, dateStr);
        
        % Submit the request
        response = webread(requestURL);
        
        % Extract precipitation data
        precipitation = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.precip;
        windspeed = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.wind_speed;
        feelslike = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.feelslike;
        temperature = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.temperature;
        avgtemp = response.historical.(['x' strrep(dateStr,'-','_')]).avgtemp;
        totalsnow = response.historical.(['x' strrep(dateStr,'-','_')]).totalsnow;
        sunhours = response.historical.(['x' strrep(dateStr,'-','_')]).sunhour;
        sunriseTimeStr = response.historical.(['x' strrep(dateStr,'-','_')]).astro.sunrise;
        sunsetTimeStr = response.historical.(['x' strrep(dateStr,'-','_')]).astro.sunset;
        
        precipitationValues(ix) = precipitation;
        windspeedValues(ix) = windspeed;
        feelslikeValues(ix) = feelslike;
        temperatureValues(ix) = avgtemp;
        snowValues(ix) = totalsnow;
        sunhoursValues(ix) = sunhours;
        sunriseTime = datetime(sunriseTimeStr, 'InputFormat', 'hh:mm a');
        sunsetTime = datetime(sunsetTimeStr, 'InputFormat', 'hh:mm a');
        sunriseValues{ix} = dateshift(currentDate, 'start', 'day') + timeofday(sunriseTime);
        sunsetValues{ix} = dateshift(currentDate, 'start', 'day') + timeofday(sunsetTime);
    end
    
    precipitationValues = precipitationValues.';
    windspeedValues = windspeedValues.';
    feelslikeValues = feelslikeValues.';
    snowValues = snowValues.';
    temperatureValues = temperatureValues.';
    sunriseValues = sunriseValues.';
    sunsetValues = sunsetValues.';
    sunhoursValues = sunhoursValues.';
end

function out = num2sepstr(numin, format, sep)
    % NUM2SEPSTR Convert to string with separation at thousands.
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

function plotCombinedWeekly(locationData, weatherData, analysis, plots, style)
    figure('Position', [408 126 1132 921]);
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    if plots.showWeather
        yyaxis right
        weatherHandles = plotWeeklyWeatherData(weatherData, plots, style);
    end
    
    % Plot traffic data on left axis
    yyaxis left
    hold on
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate weekly totals
        weeklyData = calculateWeeklyTotals(data, analysis);
        
        % Plot raw counts if enabled
        if plots.showRawCounts
            h1 = plot(weeklyData.weekStarts, weeklyData.rawCounts, '-', ...
                'LineWidth', style.plotLineWidth * 0.5, ...
                'Color', locationInfo.plotColor, ...
                'DisplayName', sprintf('%s Raw ( Min = %s ; Max = %s ; Total = %s )', ...
                    locationInfo.name, ...
                    num2sepstr(min(weeklyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(weeklyData.rawCounts), '%.0f'), ...
                    num2sepstr(sum(weeklyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
        
        % Plot adjusted counts if enabled
        if plots.showAdjustedCounts
            h2 = plot(weeklyData.weekStarts, weeklyData.adjustedCounts, '--', ...
                'LineWidth', style.plotLineWidth * 0.3, ...
                'Color', locationInfo.plotColor * 0.7, ...
                'DisplayName', sprintf('%s Adjusted ( Min = %s ; Max = %s ; Total = %s )', ...
                    locationInfo.name, ...
                    num2sepstr(min(weeklyData.adjustedCounts), '%.0f'), ...
                    num2sepstr(max(weeklyData.adjustedCounts), '%.0f'), ...
                    num2sepstr(sum(weeklyData.adjustedCounts), '%.0f')));
            plotHandles = [plotHandles, h2];
        end
    end
    
    % Format plot
    formatCombinedWeeklyPlot(analysis, plots, style, weatherData);

    if plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end

    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'north', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function weeklyData = calculateWeeklyTotals(locationDataStruct, analysis)
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Calculate weekly totals using the existing yearWeekKey
    groupedData = groupsummary(data, 'yearWeekKey', 'sum', analysis.modeString);
    
    % Calculate adjusted weekly totals (for times up to cutoff)
    truncatedData = data(timeofday(data.('Date and Time (Local)')) <= analysis.truncationCutoffTime, :);
    adjustedGrouped = groupsummary(truncatedData, 'yearWeekKey', 'sum', 'AdjustedCountsUptimeDaylight');
    
    % Also get daylight data counts per week to identify weeks with no daylight data
    daylightGrouped = groupsummary(data, 'yearWeekKey', 'sum', 'Daylight');
    
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
    
    % Match adjusted counts to raw count weeks
    [~, ic, id] = intersect(groupedData.yearWeekKey, adjustedGrouped.yearWeekKey);
    adjustedCounts = nan(size(weeklyData.rawCounts));
    adjustedCounts(ic) = adjustedGrouped.sum_AdjustedCountsUptimeDaylight(id);
    
    % Ensure adjusted is at least as large as raw
    weeklyData.adjustedCounts = max(weeklyData.rawCounts, adjustedCounts);
    
    % Filter out weeks with no daylight data (sum_Daylight = 0)
    [~, ie, if_] = intersect(groupedData.yearWeekKey, daylightGrouped.yearWeekKey);
    daylightCounts = zeros(size(weeklyData.rawCounts));
    daylightCounts(ie) = daylightGrouped.sum_Daylight(if_);
    
    % Keep only weeks that have some daylight data
    validWeeks = daylightCounts > 0;
    weeklyData.yearWeekKeys = weeklyData.yearWeekKeys(validWeeks);
    weeklyData.weekStarts = weeklyData.weekStarts(validWeeks);
    weeklyData.rawCounts = weeklyData.rawCounts(validWeeks);
    weeklyData.adjustedCounts = weeklyData.adjustedCounts(validWeeks);
end

function weatherHandles = plotWeeklyWeatherData(weatherData, plots, style)
    % Aggregate weather data to weekly averages
    weeklyWeatherData = aggregateWeatherToWeekly(weatherData);
    
    % Plot temperature line
    h1 = plot(weeklyWeatherData.weekStarts, weeklyWeatherData.avgTemperature, '-', ...
        'LineWidth', style.plotLineWidth * 0.5, ...
        'Color', [0 0.4471 0.7412 0.3], ...
        'DisplayName', 'Avg Temperature (°C)');
    
    weatherHandles = h1;
    
    % Add precipitation bubbles if enabled
    if plots.showPrecipitationBubbles
        hold on
        h2 = bubblechart(weeklyWeatherData.weekStarts, weeklyWeatherData.avgTemperature, weeklyWeatherData.totalPrecipitation, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.3, ...
            'DisplayName', 'Total Precipitation (bubble size)');
        weatherHandles = [weatherHandles, h2];
        hold off
    end
    
    ylabel('Avg Temperature (°C)', 'FontSize', style.labelFontSize, ...
        'Color', [1 0 0 0.5] + 0.5.*[0 1 1 0], 'FontWeight', 'bold');
end

function weeklyWeatherData = aggregateWeatherToWeekly(weatherData)
    % Create week grouping for weather data
    tempTable = table(weatherData.dates, weatherData.temperature, weatherData.precipitation, ...
        'VariableNames', {'dates', 'temperature', 'precipitation'});
    
    % Add week start dates
    tempTable.weekStarts = dateshift(dateshift(tempTable.dates,'dayofweek','Monday','previous'),'start','day');
    
    % Group by week
    weeklyGrouped = groupsummary(tempTable, 'weekStarts', {'mean', 'sum'}, {'temperature', 'precipitation'});
    
    weeklyWeatherData = struct();
    weeklyWeatherData.weekStarts = weeklyGrouped.weekStarts;
    weeklyWeatherData.avgTemperature = weeklyGrouped.mean_temperature;
    weeklyWeatherData.totalPrecipitation = weeklyGrouped.sum_precipitation;
end

function formatCombinedWeeklyPlot(analysis, plots, style, weatherData)
    % Left axis formatting
    yyaxis left
    ylabel(['Total ' analysis.modeDisplayString ' per Week'], ...
        'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    if plots.showWeather
        ax.YAxis(2).Color = [0 0.4471 0.7412];
    end
    
    % Title and formatting
    title(['Weekly ' analysis.modeDisplayString], ...
        'FontSize', style.titleFontSize);
    
    xlabel('Week Starting Date', 'FontSize', style.labelFontSize);
    set(gca, 'Color', style.axisBackgroundColor);
    grid on;
    xtickangle(45);
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    if ~isempty(weatherData)
        xlim([weatherData.dates(1) weatherData.dates(end)]);
    end
end

function plotCombinedMonthly(locationData, weatherData, analysis, plots, style)
    figure('Position', [408 126 1132 921]);
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    

    % Plot traffic data on left axis
    yyaxis left
    hold on

    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;

        % Calculate monthly totals
        monthlyData = calculateMonthlyTotals(data, analysis);

        % Plot raw counts if enabled
        if plots.showRawCounts
            h1 = plot(monthlyData.monthStarts, monthlyData.rawCounts, '-', ...
                'LineWidth', style.plotLineWidth * 0.5, ...
                'MarkerSize', style.plotLineWidth * 2, ...
                'Color', locationInfo.plotColor, ...
                'MarkerFaceColor', locationInfo.plotColor, ...
                'DisplayName', sprintf('%s Raw ( Min = %s ; Max = %s ; Total = %s )', ...
                locationInfo.name, ...
                num2sepstr(min(monthlyData.rawCounts), '%.0f'), ...
                num2sepstr(max(monthlyData.rawCounts), '%.0f'), ...
                num2sepstr(sum(monthlyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end

        % Plot adjusted counts if enabled
        if plots.showAdjustedCounts
            h2 = plot(monthlyData.monthStarts, monthlyData.adjustedCounts, '--', ...
                'LineWidth', style.plotLineWidth * 0.3, ...
                'MarkerSize', style.plotLineWidth * 1.5, ...
                'Color', locationInfo.plotColor * 0.7, ...
                'MarkerFaceColor', locationInfo.plotColor * 0.7, ...
                'DisplayName', sprintf('%s Adjusted ( Min = %s ; Max = %s ; Total = %s )', ...
                locationInfo.name, ...
                num2sepstr(min(monthlyData.adjustedCounts), '%.0f'), ...
                num2sepstr(max(monthlyData.adjustedCounts), '%.0f'), ...
                num2sepstr(sum(monthlyData.adjustedCounts), '%.0f')));
            plotHandles = [plotHandles, h2];
        end
    end

    % Plot weather on right axis if enabled
    if plots.showWeather
        yyaxis right
        weatherHandles = plotMonthlyWeatherData(weatherData, plots, style, monthlyData.monthStarts);
    end

    % Format plot
    formatCombinedMonthlyPlot(analysis, plots, style, weatherData);

    if plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end

    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'north', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function monthlyData = calculateMonthlyTotals(locationDataStruct, analysis)
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Create monthly grouping
    data.monthStartDateTimes = dateshift(data.('Date and Time (Local)'), 'start', 'month');
    
    % Determine which months to include based on includePartialMonths flag
    if analysis.includePartialMonths
        % Include all months that have any data
        monthsToInclude = unique(data.monthStartDateTimes);
    else
        % Only include complete months based on actual data coverage
        monthsToInclude = getCompleteMonthsFromData(data);
    end
    
    % Filter data to only include selected months
    filteredData = data(ismember(data.monthStartDateTimes, monthsToInclude), :);
    
    if isempty(filteredData)
        % Return empty structure if no valid months
        monthlyData = struct();
        monthlyData.monthStarts = datetime.empty;
        monthlyData.rawCounts = [];
        monthlyData.adjustedCounts = [];
        return;
    end
    
    % Calculate monthly totals
    groupedData = groupsummary(filteredData, 'monthStartDateTimes', 'sum', analysis.modeString);
    
    % Calculate adjusted monthly totals (for times up to cutoff)
    truncatedData = filteredData(timeofday(filteredData.('Date and Time (Local)')) <= analysis.truncationCutoffTime, :);
    adjustedGrouped = groupsummary(truncatedData, 'monthStartDateTimes', 'sum', 'AdjustedCountsUptimeDaylight');
    
    % Also get daylight data counts per month to identify months with no daylight data
    daylightGrouped = groupsummary(filteredData, 'monthStartDateTimes', 'sum', 'Daylight');
    
    % Combine results
    monthlyData = struct();
    monthlyData.monthStarts = groupedData.monthStartDateTimes;
    
    % Use the original column name with groupsummary's 'sum_' prefix
    sumColumnName = ['sum_' analysis.modeString];
    monthlyData.rawCounts = groupedData.(sumColumnName);
    
    % Match adjusted counts to raw count months
    [~, ia, ib] = intersect(groupedData.monthStartDateTimes, adjustedGrouped.monthStartDateTimes);
    adjustedCounts = nan(size(monthlyData.rawCounts));
    adjustedCounts(ia) = adjustedGrouped.sum_AdjustedCountsUptimeDaylight(ib);
    
    % Ensure adjusted is at least as large as raw
    monthlyData.adjustedCounts = max(monthlyData.rawCounts, adjustedCounts);
    
    % Filter out months with no daylight data (sum_Daylight = 0)
    [~, ic, id] = intersect(groupedData.monthStartDateTimes, daylightGrouped.monthStartDateTimes);
    daylightCounts = zeros(size(monthlyData.rawCounts));
    daylightCounts(ic) = daylightGrouped.sum_Daylight(id);
    
    % Keep only months that have some daylight data
    validMonths = daylightCounts > 0;
    monthlyData.monthStarts = monthlyData.monthStarts(validMonths);
    monthlyData.rawCounts = monthlyData.rawCounts(validMonths);
    monthlyData.adjustedCounts = monthlyData.adjustedCounts(validMonths);
end

function completeMonths = getCompleteMonthsFromData(data)
    % Determine which months have complete coverage based on actual data
    
    % Get the actual date range of the data
    minDate = min(data.('Date and Time (Local)'));
    maxDate = max(data.('Date and Time (Local)'));
    
    % Get all unique months in the data
    allMonths = unique(data.monthStartDateTimes);
    
    % Check each month for completeness
    completeMonths = [];
    
    for i = 1:length(allMonths)
        monthStart = allMonths(i);
        monthEnd = dateshift(monthStart, 'end', 'month');
        
        % A month is complete if:
        % 1. The month start is after or equal to the first full month of data, OR
        % 2. The month end is before or equal to the last full month of data
        
        % First, check if this is the first or last month in the data
        isFirstMonth = (monthStart == dateshift(minDate, 'start', 'month'));
        isLastMonth = (monthStart == dateshift(maxDate, 'start', 'month'));
        
        if isFirstMonth && isLastMonth
            % Special case: only one month of data
            % Include it only if we have data spanning most of the month
            daysCovered = days(maxDate - minDate) + 1;
            daysInMonth = day(monthEnd);
            if daysCovered >= daysInMonth * 0.8  % At least 80% of the month
                completeMonths = [completeMonths; monthStart];
            end
        elseif isFirstMonth
            % First month: include only if data starts near the beginning
            if day(minDate) <= 3  % Data starts within first 3 days of month
                completeMonths = [completeMonths; monthStart];
            end
        elseif isLastMonth
            % Last month: include only if data goes near the end
            if day(maxDate) >= day(monthEnd) - 2  % Data goes to within 2 days of month end
                completeMonths = [completeMonths; monthStart];
            end
        else
            % Middle months: these should be complete
            completeMonths = [completeMonths; monthStart];
        end
    end
end

function weatherHandles = plotMonthlyWeatherData(weatherData, plots, style, monthStarts)
    % Aggregate weather data to monthly averages
    monthlyWeatherData = aggregateWeatherToMonthly(weatherData);
    
    % Prune to match range of counting data bins
    keepIx = find(monthlyWeatherData.monthStarts>=min(monthStarts) & monthlyWeatherData.monthStarts<=max(monthStarts));

    % Plot temperature line
    h1 = plot(monthlyWeatherData.monthStarts(keepIx), monthlyWeatherData.avgTemperature(keepIx), '-', ...
        'LineWidth', style.plotLineWidth * 0.5, ...
        'Color', [0 0.4471 0.7412 0.3], ...
        'DisplayName', 'Avg Temperature (°C)');
    
    weatherHandles = h1;
    
    % Add precipitation bubbles if enabled
    if plots.showPrecipitationBubbles
        hold on
        h2 = bubblechart(monthlyWeatherData.monthStarts(keepIx), monthlyWeatherData.avgTemperature(keepIx), monthlyWeatherData.totalPrecipitation(keepIx), ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.3, ...
            'DisplayName', 'Total Precipitation (bubble size)');
        weatherHandles = [weatherHandles, h2];
        hold off
    end
    
    ylabel('Avg Temperature (°C)', 'FontSize', style.labelFontSize, ...
        'Color', [1 0 0 0.5] + 0.5.*[0 1 1 0], 'FontWeight', 'bold');
end

function monthlyWeatherData = aggregateWeatherToMonthly(weatherData)
    % Create month grouping for weather data
    tempTable = table(weatherData.dates, weatherData.temperature, weatherData.precipitation, ...
        'VariableNames', {'dates', 'temperature', 'precipitation'});
    
    % Add month start dates
    tempTable.monthStarts = dateshift(tempTable.dates, 'start', 'month');
    
    % Group by month
    monthlyGrouped = groupsummary(tempTable, 'monthStarts', {'mean', 'sum'}, {'temperature', 'precipitation'});
    
    monthlyWeatherData = struct();
    monthlyWeatherData.monthStarts = monthlyGrouped.monthStarts;
    monthlyWeatherData.avgTemperature = monthlyGrouped.mean_temperature;
    monthlyWeatherData.totalPrecipitation = monthlyGrouped.sum_precipitation;
end

function formatCombinedMonthlyPlot(analysis, plots, style, weatherData)
    % Left axis formatting
    yyaxis left
    ylabel(['Total ' analysis.modeDisplayString ' per Month'], ...
        'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    if plots.showWeather
        ax.YAxis(2).Color = [0 0.4471 0.7412];
    end
    
    % Title and formatting
    if analysis.includePartialMonths
        monthsTitle = 'Monthly ' + string(analysis.modeDisplayString) + ' (with partial months)';
    else
        monthsTitle = 'Monthly ' + string(analysis.modeDisplayString) + ' on rue de Terrebonne';
    end
    title(monthsTitle, 'FontSize', style.titleFontSize);
    
    xlabel('Month', 'FontSize', style.labelFontSize);
    set(gca, 'Color', style.axisBackgroundColor);
    grid on;
    xtickangle(45);
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    % Set x-axis limits based on traffic data range (not weather data)
    if ~isempty(weatherData)
        xlim([weatherData.dates(1) weatherData.dates(end)]);
    end

end

function plotMultiModalDaily(locationData, weatherData, analysis, plots, style, multiModal)
    figure('Position', [408 126 1132 921]);
    
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    if multiModal.plotWeather && plots.showWeather
        yyaxis right
        weatherHandles = plotWeatherData(weatherData, plots, style);
    end
    
    % Plot traffic data on left axis
    yyaxis left
    hold on
    
    % Get the specified location data
    locationFieldName = matlab.lang.makeValidName(multiModal.location);
    if ~isfield(locationData, locationFieldName)
        error('Location "%s" not found in data', multiModal.location);
    end
    
    data = locationData.(locationFieldName);
    
    % Plot each mode
    for i = 1:length(multiModal.modes)
        currentMode = multiModal.modes{i};
        currentModeDisplay = multiModal.modeDisplayNames{i};
        currentColor = multiModal.modeColors{i};
        
        % Create temporary analysis structure for this mode
        tempAnalysis = analysis;
        tempAnalysis.modeString = currentMode;
        tempAnalysis.modeDisplayString = currentModeDisplay;
        
        % Calculate daily totals for this mode
        dailyData = calculateDailyTotals(data, tempAnalysis);
        
        if ~isempty(dailyData.dates)
            % Plot raw counts
            h1 = plot(dailyData.dates, dailyData.rawCounts, '-', ...
                'LineWidth', style.plotLineWidth * 0.5, ...
                'Color', currentColor, ...
                'DisplayName', sprintf('%s ( Min = %s ; Max = %s ; Total = %s )', ...
                    currentModeDisplay, ...
                    num2sepstr(min(dailyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(dailyData.rawCounts), '%.0f'), ...
                    num2sepstr(sum(dailyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
    end
    
    % Format plot
    formatMultiModalPlot('Daily', multiModal, plots, style, weatherData);
    
    if multiModal.plotWeather && plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'north', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function plotMultiModalWeekly(locationData, weatherData, analysis, plots, style, multiModal)
    figure('Position', [408 126 1132 921]);
    
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    if multiModal.plotWeather && plots.showWeather
        yyaxis right
        weatherHandles = plotWeeklyWeatherData(weatherData, plots, style);
    end
    
    % Plot traffic data on left axis
    yyaxis left
    hold on
    
    % Get the specified location data
    locationFieldName = matlab.lang.makeValidName(multiModal.location);
    data = locationData.(locationFieldName);
    
    % Plot each mode
    for i = 1:length(multiModal.modes)
        currentMode = multiModal.modes{i};
        currentModeDisplay = multiModal.modeDisplayNames{i};
        currentColor = multiModal.modeColors{i};
        
        % Create temporary analysis structure for this mode
        tempAnalysis = analysis;
        tempAnalysis.modeString = currentMode;
        tempAnalysis.modeDisplayString = currentModeDisplay;
        
        % Calculate weekly totals for this mode
        weeklyData = calculateWeeklyTotals(data, tempAnalysis);
        
        if ~isempty(weeklyData.weekStarts)
            % Plot raw counts
            h1 = plot(weeklyData.weekStarts, weeklyData.rawCounts, '-', ...
                'LineWidth', style.plotLineWidth * 0.5, ...
                'Color', currentColor, ...
                'DisplayName', sprintf('%s ( Min = %s ; Max = %s ; Total = %s )', ...
                    currentModeDisplay, ...
                    num2sepstr(min(weeklyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(weeklyData.rawCounts), '%.0f'), ...
                    num2sepstr(sum(weeklyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
    end
    
    % Format plot
    formatMultiModalPlot('Weekly', multiModal, plots, style, weatherData);
    
    if multiModal.plotWeather && plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'north', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function plotMultiModalMonthly(locationData, weatherData, analysis, plots, style, multiModal)
    figure('Position', [408 126 1132 921]);
    
    plotHandles = [];
    
    
    
    % Plot traffic data on left axis
    yyaxis left
    hold on
    
    % Get the specified location data
    locationFieldName = matlab.lang.makeValidName(multiModal.location);
    data = locationData.(locationFieldName);
    
    % Plot each mode
    for i = 1:length(multiModal.modes)
        currentMode = multiModal.modes{i};
        currentModeDisplay = multiModal.modeDisplayNames{i};
        currentColor = multiModal.modeColors{i};
        
        % Create temporary analysis structure for this mode
        tempAnalysis = analysis;
        tempAnalysis.modeString = currentMode;
        tempAnalysis.modeDisplayString = currentModeDisplay;
        
        % Calculate monthly totals for this mode
        monthlyData = calculateMonthlyTotals(data, tempAnalysis);
        
        if ~isempty(monthlyData.monthStarts)
            % Plot raw counts
            h1 = plot(monthlyData.monthStarts, monthlyData.rawCounts, '-', ...
                'LineWidth', style.plotLineWidth * 0.5, ...
                'MarkerSize', style.plotLineWidth * 1.5, ...
                'Color', currentColor, ...
                'MarkerFaceColor', currentColor, ...
                'DisplayName', sprintf('%s ( Min = %s ; Max = %s ; Total = %s )', ...
                    currentModeDisplay, ...
                    num2sepstr(min(monthlyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(monthlyData.rawCounts), '%.0f'), ...
                    num2sepstr(sum(monthlyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
    end
    
    % Plot weather on right axis if enabled
    if multiModal.plotWeather && plots.showWeather
        yyaxis right
        weatherHandles = plotMonthlyWeatherData(weatherData, plots, style, monthlyData.monthStarts);
    end

    % Format plot
    formatMultiModalPlot('Monthly', multiModal, plots, style, weatherData);
    
    if multiModal.plotWeather && plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'north', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function formatMultiModalPlot(timeScale, multiModal, plots, style, weatherData)
    % Left axis formatting
    yyaxis left
    ylabel(['Total Counts per ' timeScale], ...
        'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    if multiModal.plotWeather && plots.showWeather
        ax.YAxis(2).Color = [0 0.4471 0.7412];
    end
    
    % Title and formatting
    title([timeScale ' Traffic Counts by Mode (' multiModal.location ')'], ...
        'FontSize', style.titleFontSize);
    
    if strcmp(timeScale, 'Daily')
        xlabel('Date', 'FontSize', style.labelFontSize);
    elseif strcmp(timeScale, 'Weekly')
        xlabel('Week Starting Date', 'FontSize', style.labelFontSize);
    else
        xlabel('Month', 'FontSize', style.labelFontSize);
    end
    
    set(gca, 'Color', style.axisBackgroundColor);
    grid on;
    xtickangle(45);
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    if ~isempty(weatherData)
        xlim([weatherData.dates(1) weatherData.dates(end)]);
    end
end

