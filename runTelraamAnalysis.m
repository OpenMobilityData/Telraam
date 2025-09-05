%% Modular Telraam Analysis - Single Script Version
% This version can be run as a single script while maintaining modularity

clear all; close all; clc;

%% Configuration Setup

westernSegmentName = 'rue de Terrebonne @ King Edward';
easternSegmentName = 'rue de Terrebonne @ Draper';

% Locations to analyze
locations = {
    struct('name', easternSegmentName, ...
           'fileStem2024', 'raw-data-2024East60', ...
           'fileStem2025', 'raw-data-2025East60', ...
           'plotColor', [0 0 1]);  % Blue
    struct('name', westernSegmentName, ...
           'fileStem2024', 'raw-data-2024West60', ...
           'fileStem2025', 'raw-data-2025West60', ...
           'plotColor', [0 0 0]);  % Black
};

%locations = {locations{1}};

% Analysis parameters

modeString = 'Bike Total'; modeDisplayString = 'Bike Counts';
%modeString = 'Pedestrian Total'; modeDisplayString = 'Pedestrian Counts';
%modeString = 'Car Total'; modeDisplayString = 'Car Counts';
%modeString = 'Large vehicle Total'; modeDisplayString = 'Heavy Truck Counts';

analysis = struct( ...
    'startTime', datetime(2024,08,15,00,00,01), ...
    'endTime', datetime(2025,09,03,23,59,59), ...
    'modeString', modeString, ...
    'modeDisplayString', modeDisplayString, ...
    'uptimeThreshold', 0.0, ...
    'maxUptimeCorrection', 1.0, ...
    'truncationCutoffTime', timeofday(datetime('today')+hours(15)), ...
    'daylightCorrectionRatioWD', 1.65, ...
    'daylightCorrectionRatioWE', 1.37, ...
    'includePartialMonths', false, ...
    'outlierDetection', struct( ...
        'enabled', true, ...
        'method', 'iqr', ...  % 'iqr', 'zscore', or 'manual'
        'threshold', 3.0, ...  % IQR multiplier or Z-score threshold
        'reportOnly', false, ...  % If true, report but don't remove
        'manualExclusions', [] ...  % Manual datetime exclusions
    ) ...
);

%dateSpan = 'winter';
%dateSpan = 'springSummer';
%dateSpan = 'lastWeek';
%dateSpan = 'lastMonth';

if exist('dateSpan', 'var')
    if strcmp(dateSpan,'winter')
        analysis.startTime = datetime(2024,11,16,00,00,01);
        analysis.endTime = datetime(2025,03,31,23,59,59);
    elseif strcmp(dateSpan,'springSummer')
        analysis.startTime = datetime(2025,04,01,23,59,59);
    elseif strcmp(dateSpan,'lastWeek')
        analysis.startTime = analysis.endTime - days(7) + seconds(2);
    elseif strcmp(dateSpan,'lastMonth')
        analysis.startTime = dateshift(analysis.endTime,'start','month') + seconds(2);
    end
end

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

if plots.showWeather
    % Check if cached data exists
    if exist(cacheFilename, 'file')
        fprintf('Loading cached weather data from %s...\n', cacheFilename);
        load(cacheFilename, 'weatherData');
        fprintf('Loaded weather data for %d days from cache.\n', length(weatherData.dates));
    else
        % Get weather data from API
        fprintf('Getting weather data for %d days from API...\n', length(uniqueDays));
        [precipitationData, averageTemperatureData, minTemperatureData, maxTemperatureData, sunriseData, sunsetData, sunhoursData, snowData, windspeedData, feelslikeData, visibilityData] = ...
            getWeatherstackData('Montreal', dailyNoonTimes);

        % Store weather data correctly - extract the actual vectors
        weatherData = struct();
        weatherData.dates = uniqueDays;
        weatherData.precipitation = precipitationData;
        weatherData.temperature = averageTemperatureData;
        %weatherData.temperature = maxTemperatureData;
        weatherData.feelslike = feelslikeData;
        weatherData.visibility = visibilityData;
        weatherData.windspeed = windspeedData;
        weatherData.sunrise = sunriseData;
        weatherData.sunset = sunsetData;
        weatherData.sunhours = sunhoursData;
        weatherData.snow = snowData;

        % Save to cache
        fprintf('Saving weather data to cache: %s\n', cacheFilename);
        save(cacheFilename, 'weatherData');
    end
end

%% Generate Hourly Raw Count Plots
filteredLocationData = filterNightData(locationData, weatherData);
plotCombinedHourlyRaw(filteredLocationData, weatherData, analysis, plots, style);

%% Generate Hourly Count Histogram
plotHourlyCountHistogram(filteredLocationData, analysis, style);

%% Analyze Zero Count Intervals
analyzeZeroCountIntervals(filteredLocationData, analysis);

%% Optional: Visualize Zero Intervals
% Uncomment to create a visual plot of zero intervals
%plotZeroIntervalAnalysis(filteredLocationData, analysis, style);

%% Generate Combined Daily Plot
plotCombinedDaily(locationData, weatherData, analysis, plots, style);

%% Generate Combined Weekly Plot
plotCombinedWeekly(locationData, weatherData, analysis, plots, style);

%% Generate Day-of-Week Pattern Plots
plotDayOfWeekPatterns(locationData, analysis, plots, style);

%% Generate Combined Monthly Plot
plotCombinedMonthly(locationData, weatherData, analysis, plots, style);

%% Generate Hourly Pattern Plots (if enabled)
plotHourlyPatterns(locationData, analysis, plots, style);

%% Generate Temperature Scatter Plot
plotTemperatureScatterWeekly(locationData, weatherData, analysis, plots, style);
%plotTemperatureScatterSeasonal(locationData, weatherData, analysis, plots, style);

%% Generate Visibility Scatter Plot
%plotVisibilityScatterWeekly(locationData, weatherData, analysis, plots, style);

%% Multivariate Weather Analysis (MOVED TO SEPARATE SCRIPT)
% performMultivariateWeatherAnalysis(locationData, weatherData, analysis, style);

%% Generate Scatter Plot comparing counts at the two locations
plotLocationCorrelation(locationData, analysis, plots, style, 'daily');

% Generate Bike vs Other Modalities Correlation Plots
plotBikeModalityCorrelation(locationData, analysis, plots, style);

%% Generate Modality Pie and Bar Charts
plotModalityPieCharts(locationData, analysis, style);
plotModalityBarChart(locationData, analysis, style);

%% Generate Multi-Modal Plots (if enabled)
if multiModal.enabled
    plotMultiModalDaily(locationData, weatherData, analysis, plots, style, multiModal);
    plotMultiModalWeekly(locationData, weatherData, analysis, plots, style, multiModal);
    plotMultiModalMonthly(locationData, weatherData, analysis, plots, style, multiModal);
    %plotMultiModalHourlyPatterns(locationData, analysis, plots, style, multiModal);
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

    % Apply outlier detection and removal
    if isfield(analysis, 'outlierDetection') && analysis.outlierDetection.enabled
        processedTable = detectAndHandleOutliers(processedTable, analysis);
    end

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
function [precipitationValues, averageTemperatureValues, minTemperatureValues, maxTemperatureValues, sunriseValues, sunsetValues, sunhoursValues, snowValues, windspeedValues, feelslikeValues, visibilityValues] = getWeatherstackData(location, dates)
    % Define the base URL for Weatherstack API
    baseURL = 'http://api.weatherstack.com/historical';
    
    % Replace with your actual Weatherstack access key (if needed)
    %accessKey = 'fe2d67122ba14cc9e0b2c931f6105b4b'; % Leave empty if access key is not required
    accessKey = '9cbcc044ac6abbe7f78b2fe3e26ed4b9'; % Leave empty if access key is not required
    
    
    % Iterate over each date in the range
    numDates = length(dates);
    precipitationValues = [];
    windspeedValues = [];
    feelslikeValues = [];
    visibilityValues = [];
    snowValues = [];
    averageTemperatureValues = [];
    minTemperatureValues = [];
    maxTemperatureValues = [];
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
        visibility = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.visibility;
        temperature = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.temperature;
        avgtemp = response.historical.(['x' strrep(dateStr,'-','_')]).avgtemp;
        mintemp = response.historical.(['x' strrep(dateStr,'-','_')]).mintemp;
        maxtemp = response.historical.(['x' strrep(dateStr,'-','_')]).maxtemp;
        totalsnow = response.historical.(['x' strrep(dateStr,'-','_')]).totalsnow;
        sunhours = response.historical.(['x' strrep(dateStr,'-','_')]).sunhour;
        sunriseTimeStr = response.historical.(['x' strrep(dateStr,'-','_')]).astro.sunrise;
        sunsetTimeStr = response.historical.(['x' strrep(dateStr,'-','_')]).astro.sunset;
        
        precipitationValues(ix) = precipitation;
        windspeedValues(ix) = windspeed;
        feelslikeValues(ix) = feelslike;
        visibilityValues(ix) = visibility;
        averageTemperatureValues(ix) = avgtemp;
        minTemperatureValues(ix) = mintemp;
        maxTemperatureValues(ix) = maxtemp;
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
    visibilityValues = visibilityValues.';
    snowValues = snowValues.';
    averageTemperatureValues = averageTemperatureValues.';
    minTemperatureValues = minTemperatureValues.';
    maxTemperatureValues = maxTemperatureValues.';
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

    % Check if any location has monthly data before proceeding
    locationNames = fieldnames(locationData);
    hasMonthlyData = false;
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        monthlyData = calculateMonthlyTotals(data, analysis);
        
        if ~isempty(monthlyData.monthStarts)
            hasMonthlyData = true;
            break;
        end
    end
    
    if ~hasMonthlyData
        fprintf('Skipping monthly plot: No complete months found in date range (%s to %s).\n', ...
            datestr(analysis.startTime, 'dd-mmm-yyyy'), datestr(analysis.endTime, 'dd-mmm-yyyy'));
        if ~analysis.includePartialMonths
            fprintf('Consider setting analysis.includePartialMonths = true to include partial months.\n');
        end
        return;
    end

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
    
    % Check if we have any months to include
    if isempty(monthsToInclude)
        % Return empty structure if no valid months
        monthlyData = struct();
        monthlyData.monthStarts = datetime.empty(0,1);  % Ensure it's a datetime array
        monthlyData.rawCounts = double.empty(0,1);      % Ensure it's a numeric array
        monthlyData.adjustedCounts = double.empty(0,1); % Ensure it's a numeric array
        
        fprintf('Warning: No complete months found in date range. ');
        if ~analysis.includePartialMonths
            fprintf('Consider setting analysis.includePartialMonths = true.\n');
        else
            fprintf('No data available for monthly analysis.\n');
        end
        return;
    end
    
    % Filter data to only include selected months
    filteredData = data(ismember(data.monthStartDateTimes, monthsToInclude), :);
    
    if isempty(filteredData)
        % Return empty structure if no data after filtering
        monthlyData = struct();
        monthlyData.monthStarts = datetime.empty(0,1);
        monthlyData.rawCounts = double.empty(0,1);
        monthlyData.adjustedCounts = double.empty(0,1);
        
        fprintf('Warning: No data remaining after month filtering.\n');
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

    % Check if the specified location has monthly data
    locationFieldName = matlab.lang.makeValidName(multiModal.location);
    if ~isfield(locationData, locationFieldName)
        fprintf('Skipping multi-modal monthly plot: Location "%s" not found.\n', multiModal.location);
        return;
    end
    
    data = locationData.(locationFieldName);
    
    % Check if any mode has monthly data
    hasMonthlyData = false;
    for i = 1:length(multiModal.modes)
        tempAnalysis = analysis;
        tempAnalysis.modeString = multiModal.modes{i};
        monthlyData = calculateMonthlyTotals(data, tempAnalysis);
        
        if ~isempty(monthlyData.monthStarts)
            hasMonthlyData = true;
            break;
        end
    end
    
    if ~hasMonthlyData
        fprintf('Skipping multi-modal monthly plot: No complete months found for location "%s".\n', multiModal.location);
        if ~analysis.includePartialMonths
            fprintf('Consider setting analysis.includePartialMonths = true to include partial months.\n');
        end
        return;
    end

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

function plotHourlyPatterns(locationData, analysis, plots, style)
    % Plot average hourly traffic patterns for weekdays and weekends
    % Includes both grand averages and monthly segregation
    
    locationNames = fieldnames(locationData);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate hourly patterns
        hourlyData = calculateHourlyPatterns(data, analysis);
        monthlyHourlyData = calculateMonthlyHourlyPatterns(data, analysis);
        
        % Create grand average plot for each location
        plotGrandAverageHourly(hourlyData, analysis, style, locationInfo.name);
        
        % Create monthly segregated plots for weekdays and weekends
        plotMonthlyHourlyPatterns(monthlyHourlyData, analysis, style, locationInfo.name, 'Weekdays');
        plotMonthlyHourlyPatterns(monthlyHourlyData, analysis, style, locationInfo.name, 'Weekends');
    end
end

function hourlyData = calculateHourlyPatterns(locationDataStruct, analysis)
    % Calculate average hourly traffic patterns for weekdays and weekends
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Define weekdays
    weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
    weekends = {'Saturday', 'Sunday'};
    
    % Initialize arrays for weekday and weekend traffic
    weekdayTraffic = [];
    weekdayTimes = [];
    weekendTraffic = [];
    weekendTimes = [];
    
    % Process weekdays
    for ix = 1:length(weekdays)
        targetDay = weekdays{ix};
        
        % Extract rows that match the target day
        dayRows = strcmp(string(data.dayOfWeek), targetDay);
        
        % Filter the timetable for the selected day
        filteredTable = data(dayRows, :);
        
        if ~isempty(filteredTable)
            % Extract the time component and traffic counts
            timeOfDay = timeofday(filteredTable.('Date and Time (Local)'));
            trafficCounts = filteredTable.(analysis.modeString);
            
            % Add to the weekday arrays
            weekdayTraffic = [weekdayTraffic; trafficCounts];
            weekdayTimes = [weekdayTimes; timeOfDay];
        end
    end
    
    % Process weekends
    for ix = 1:length(weekends)
        targetDay = weekends{ix};
        
        % Extract rows that match the target day
        dayRows = strcmp(string(data.dayOfWeek), targetDay);
        
        % Filter the timetable for the selected day
        filteredTable = data(dayRows, :);
        
        if ~isempty(filteredTable)
            % Extract the time component and traffic counts
            timeOfDay = timeofday(filteredTable.('Date and Time (Local)'));
            trafficCounts = filteredTable.(analysis.modeString);
            
            % Add to the weekend arrays
            weekendTraffic = [weekendTraffic; trafficCounts];
            weekendTimes = [weekendTimes; timeOfDay];
        end
    end
    
    % Calculate averages for each unique time of day
    hourlyData = struct();
    
    if ~isempty(weekdayTimes)
        [uniqueWeekdayTimes, ~, idx] = unique(weekdayTimes);
        avgWeekdayTraffic = accumarray(idx, weekdayTraffic, [], @mean);
        hourlyData.weekdayTimes = uniqueWeekdayTimes;
        hourlyData.weekdayAverage = avgWeekdayTraffic;
    else
        hourlyData.weekdayTimes = [];
        hourlyData.weekdayAverage = [];
    end
    
    if ~isempty(weekendTimes)
        [uniqueWeekendTimes, ~, idx] = unique(weekendTimes);
        avgWeekendTraffic = accumarray(idx, weekendTraffic, [], @mean);
        hourlyData.weekendTimes = uniqueWeekendTimes;
        hourlyData.weekendAverage = avgWeekendTraffic;
    else
        hourlyData.weekendTimes = [];
        hourlyData.weekendAverage = [];
    end
end

% function formatHourlyPlot(analysis, style, locationName)
%     % Format the hourly pattern plot
% 
%     ylabel('Hourly Count', 'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
%     xlabel('Time of Day', 'FontSize', style.labelFontSize);
% 
%     title(['Average Hourly ' analysis.modeDisplayString ' Patterns (' locationName ')'], ...
%         'FontSize', style.titleFontSize);
% 
%     set(gca, 'Color', style.axisBackgroundColor);
%     set(gca, 'FontSize', style.axisFontSize);
%     grid on;
% 
%     % Format y-axis with separators
%     ytick_positions = yticks;
%     ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
%     yticklabels(ytick_labels);
% 
%     ylim([0 max(ylim) * 1.1]);
% 
%     % Set reasonable x-axis limits (24 hours)
%     xlim([hours(0) hours(24)]);
% 
%     % Set x-axis ticks every 4 hours
%     xticks(hours(0:4:24));
%     xticklabels({'00:00', '04:00', '08:00', '12:00', '16:00', '20:00', '24:00'});
% end

function plotGrandAverageHourly(hourlyData, analysis, style, locationName)
    % Plot the grand average hourly patterns (original functionality)
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    plotHandles = [];
    
    % Plot weekday pattern
    if ~isempty(hourlyData.weekdayTimes)
        h1 = plot(hourlyData.weekdayTimes, hourlyData.weekdayAverage, '-', ...
            'LineWidth', style.plotLineWidth, ...
            'Color', [0 0 1], ...
            'DisplayName', sprintf('Weekdays (total = %s per day)', ...
                num2sepstr(sum(hourlyData.weekdayAverage), '%.0f')));
        plotHandles = [plotHandles, h1];
    end
    
    % Plot weekend pattern
    if ~isempty(hourlyData.weekendTimes)
        h2 = plot(hourlyData.weekendTimes, hourlyData.weekendAverage, '-', ...
            'LineWidth', style.plotLineWidth, ...
            'Color', [1 0 0], ...
            'DisplayName', sprintf('Weekends (total = %s per day)', ...
                num2sepstr(sum(hourlyData.weekendAverage), '%.0f')));
        plotHandles = [plotHandles, h2];
    end
    
    % Format plot
    formatHourlyPlot(analysis, style, locationName, 'Average');
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function monthlyHourlyData = calculateMonthlyHourlyPatterns(locationDataStruct, analysis)
    % Calculate hourly patterns segregated by month
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Define weekdays and add month information
    weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
    weekends = {'Saturday', 'Sunday'};
    
    % Add month information to the data
    data.monthStartDateTimes = dateshift(data.('Date and Time (Local)'), 'start', 'month');
    uniqueMonths = unique(data.monthStartDateTimes);
    
    % Initialize storage for monthly patterns
    monthlyHourlyData = struct();
    monthlyHourlyData.months = uniqueMonths;
    monthlyHourlyData.weekdayPatterns = {};
    monthlyHourlyData.weekendPatterns = {};
    monthlyHourlyData.weekdayTimes = {};
    monthlyHourlyData.weekendTimes = {};
    
    % Process each month
    for monthIdx = 1:length(uniqueMonths)
        currentMonth = uniqueMonths(monthIdx);
        monthData = data(data.monthStartDateTimes == currentMonth, :);
        
        if isempty(monthData)
            continue;
        end
        
        % Process weekdays for this month
        weekdayTraffic = [];
        weekdayTimes = [];
        
        for dayIdx = 1:length(weekdays)
            targetDay = weekdays{dayIdx};
            dayRows = strcmp(string(monthData.dayOfWeek), targetDay);
            filteredTable = monthData(dayRows, :);
            
            if ~isempty(filteredTable)
                timeOfDay = timeofday(filteredTable.('Date and Time (Local)'));
                trafficCounts = filteredTable.(analysis.modeString);
                weekdayTraffic = [weekdayTraffic; trafficCounts];
                weekdayTimes = [weekdayTimes; timeOfDay];
            end
        end
        
        % Calculate weekday averages for this month
        if ~isempty(weekdayTimes)
            [uniqueWeekdayTimes, ~, idx] = unique(weekdayTimes);
            avgWeekdayTraffic = accumarray(idx, weekdayTraffic, [], @mean);
            monthlyHourlyData.weekdayTimes{monthIdx} = uniqueWeekdayTimes;
            monthlyHourlyData.weekdayPatterns{monthIdx} = avgWeekdayTraffic;
        else
            monthlyHourlyData.weekdayTimes{monthIdx} = [];
            monthlyHourlyData.weekdayPatterns{monthIdx} = [];
        end
        
        % Process weekends for this month
        weekendTraffic = [];
        weekendTimes = [];
        
        for dayIdx = 1:length(weekends)
            targetDay = weekends{dayIdx};
            dayRows = strcmp(string(monthData.dayOfWeek), targetDay);
            filteredTable = monthData(dayRows, :);
            
            if ~isempty(filteredTable)
                timeOfDay = timeofday(filteredTable.('Date and Time (Local)'));
                trafficCounts = filteredTable.(analysis.modeString);
                weekendTraffic = [weekendTraffic; trafficCounts];
                weekendTimes = [weekendTimes; timeOfDay];
            end
        end
        
        % Calculate weekend averages for this month
        if ~isempty(weekendTimes)
            [uniqueWeekendTimes, ~, idx] = unique(weekendTimes);
            avgWeekendTraffic = accumarray(idx, weekendTraffic, [], @mean);
            monthlyHourlyData.weekendTimes{monthIdx} = uniqueWeekendTimes;
            monthlyHourlyData.weekendPatterns{monthIdx} = avgWeekendTraffic;
        else
            monthlyHourlyData.weekendTimes{monthIdx} = [];
            monthlyHourlyData.weekendPatterns{monthIdx} = [];
        end
    end
end

function plotMonthlyHourlyPatterns(monthlyHourlyData, analysis, style, locationName, dayType)
    % Plot hourly patterns segregated by month
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    plotHandles = [];
    colorMap = lines(length(monthlyHourlyData.months));
    
    % Determine which patterns to plot
    if strcmp(dayType, 'Weekdays')
        patterns = monthlyHourlyData.weekdayPatterns;
        times = monthlyHourlyData.weekdayTimes;
    else
        patterns = monthlyHourlyData.weekendPatterns;
        times = monthlyHourlyData.weekendTimes;
    end
    
    % Plot each month's pattern
    for monthIdx = 1:length(monthlyHourlyData.months)
        if ~isempty(patterns{monthIdx}) && ~isempty(times{monthIdx})
            monthTotal = sum(patterns{monthIdx});
            
            % Use subtle colors for individual months, emphasize most recent
            if monthIdx == length(monthlyHourlyData.months)
                % Emphasize most recent month
                alpha = 1.0;
                lineWidth = style.plotLineWidth * 0.8;
            else
                % Subtle for other months
                alpha = 0.3;
                lineWidth = style.plotLineWidth * 0.5;
            end
            
            h = plot(times{monthIdx}, patterns{monthIdx}, '-', ...
                'LineWidth', lineWidth, ...
                'Color', [colorMap(monthIdx, :) alpha], ...
                'DisplayName', sprintf('%s (%s/day)', ...
                    datestr(monthlyHourlyData.months(monthIdx), 'mmm yyyy'), ...
                    num2sepstr(monthTotal, '%.0f')));
            plotHandles = [plotHandles, h];
        end
    end
    
    % Format plot
    formatHourlyPlot(analysis, style, locationName, [dayType ' by Month']);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function formatHourlyPlot(analysis, style, locationName, plotType)
    % Format the hourly pattern plot
    
    ylabel('Hourly Count', 'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    xlabel('Time of Day', 'FontSize', style.labelFontSize);
    
    title([plotType ' Hourly ' analysis.modeDisplayString ' (' locationName ')'], ...
        'FontSize', style.titleFontSize);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    % Set reasonable x-axis limits (24 hours)
    xlim([hours(0) hours(24)]);
    
    % Set x-axis ticks every 4 hours
    xticks(hours(0:4:24));
    xticklabels({'00:00', '04:00', '08:00', '12:00', '16:00', '20:00', '24:00'});
end

function inputTable = detectAndHandleOutliers(inputTable, analysis)
    % Detect and optionally remove outliers from traffic data using context-aware comparison
    % with minimum absolute difference thresholds to avoid false positives from low-count periods
    
    fprintf('\n=== Context-Aware Outlier Detection for %s ===\n', analysis.modeDisplayString);
    
    % Get the traffic data column
    trafficData = inputTable.(analysis.modeString);
    originalCount = height(inputTable);
    
    % Add temporal context columns if they don't exist
    if ~ismember('hourOfDay', inputTable.Properties.VariableNames)
        inputTable.hourOfDay = hour(inputTable.('Date and Time (Local)'));
    end
    if ~ismember('monthStartDateTimes', inputTable.Properties.VariableNames)
        inputTable.monthStartDateTimes = dateshift(inputTable.('Date and Time (Local)'), 'start', 'month');
    end
    
    % Calculate expected values for each hour/day-type/month combination
    expectedValues = calculateExpectedHourlyCounts(inputTable, analysis);
    
    % Initialize outlier flags
    outlierFlags = false(size(trafficData));
    
    % Apply manual exclusions first
    if ~isempty(analysis.outlierDetection.manualExclusions)
        for i = 1:length(analysis.outlierDetection.manualExclusions)
            excludeTime = analysis.outlierDetection.manualExclusions(i);
            manualFlags = abs(inputTable.('Date and Time (Local)') - excludeTime) < hours(1);
            outlierFlags = outlierFlags | manualFlags;
            if any(manualFlags)
                fprintf('Manual exclusion: %s\n', datestr(excludeTime));
            end
        end
    end
    
    % Apply context-aware statistical outlier detection with absolute difference filtering
    switch lower(analysis.outlierDetection.method)
        case 'iqr'
            statisticalOutliers = detectContextualIQROutliers(inputTable, expectedValues, analysis);
        case 'zscore'
            statisticalOutliers = detectContextualZScoreOutliers(inputTable, expectedValues, analysis);
        case 'manual'
            statisticalOutliers = false(size(trafficData));
    end
    
    % Combine manual and statistical outliers
    outlierFlags = outlierFlags | statisticalOutliers;
    
    % Report outliers with context
    reportContextualOutliers(inputTable, expectedValues, outlierFlags, analysis);
    
    % Remove outliers if not in report-only mode
    if ~analysis.outlierDetection.reportOnly && any(outlierFlags)
        inputTable = inputTable(~outlierFlags, :);
        removedCount = sum(outlierFlags);
        fprintf('Removed %d outlier row(s). Kept %d of %d rows (%.1f%%).\n', ...
            removedCount, height(inputTable), originalCount, ...
            100 * height(inputTable) / originalCount);
    elseif analysis.outlierDetection.reportOnly && any(outlierFlags)
        fprintf('Report-only mode: outliers flagged but not removed.\n');
    end
    
    fprintf('\n');
end

function expectedValues = calculateExpectedHourlyCounts(inputTable, analysis)
    % Calculate expected hourly counts for each hour/day-type/month combination
    
    % Define weekdays
    weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
    
    % Get unique combinations of month, day type, and hour
    uniqueMonths = unique(inputTable.monthStartDateTimes);
    uniqueHours = unique(inputTable.hourOfDay);
    
    % Initialize expected values array
    expectedValues = nan(height(inputTable), 1);
    
    fprintf('Calculating expected hourly patterns...\n');
    
    for monthIdx = 1:length(uniqueMonths)
        currentMonth = uniqueMonths(monthIdx);
        monthData = inputTable(inputTable.monthStartDateTimes == currentMonth, :);
        
        if isempty(monthData)
            continue;
        end
        
        % Process weekdays and weekends separately
        for isWeekdayGroup = [true, false]
            if isWeekdayGroup
                dayTypeData = monthData(ismember(monthData.dayOfWeek, weekdays), :);
                dayTypeName = 'weekday';
            else
                dayTypeData = monthData(~ismember(monthData.dayOfWeek, weekdays), :);
                dayTypeName = 'weekend';
            end
            
            if isempty(dayTypeData)
                continue;
            end
            
            % Calculate robust expected value for each hour in this month/day-type combination
            for hourIdx = 1:length(uniqueHours)
                currentHour = uniqueHours(hourIdx);
                hourData = dayTypeData(dayTypeData.hourOfDay == currentHour, :);
                
                if height(hourData) >= 2  % Need at least 2 observations for meaningful average
                    hourCounts = hourData.(analysis.modeString);
                    hourCounts = hourCounts(~isnan(hourCounts));  % Remove NaN values
                    
                    if length(hourCounts) >= 2
                        % Use robust statistics to calculate expected value
                        robustExpected = calculateRobustExpected(hourCounts);
                        
                        % Find all matching rows in the original table and set expected value
                        matchingRows = (inputTable.monthStartDateTimes == currentMonth) & ...
                                       (inputTable.hourOfDay == currentHour) & ...
                                       (ismember(inputTable.dayOfWeek, weekdays) == isWeekdayGroup);
                        
                        expectedValues(matchingRows) = robustExpected;
                    end
                end
            end
        end
    end
    
    % Report summary of expected value calculation with diagnostics
    validExpected = ~isnan(expectedValues);
    fprintf('Calculated expected values for %d of %d observations (%.1f%%)\n', ...
        sum(validExpected), length(expectedValues), 100 * sum(validExpected) / length(expectedValues));
    
    if sum(validExpected) > 0
        fprintf('Expected value range: %.1f to %.1f (median: %.1f)\n', ...
            min(expectedValues(validExpected)), max(expectedValues(validExpected)), ...
            median(expectedValues(validExpected)));
    end
    
    % Diagnostic: Show details for May weekend 18:00 if it exists in the data
    mayMask = month(inputTable.('Date and Time (Local)')) == 5 & ...
              year(inputTable.('Date and Time (Local)')) == 2025;
    if any(mayMask)
        diagnosticMayWeekend18(inputTable, analysis, mayMask);
    end
end

function outliers = detectContextualIQROutliers(inputTable, expectedValues, analysis)
    % Detect outliers using absolute differences with minimum thresholds
    
    trafficData = inputTable.(analysis.modeString);
    
    % Calculate absolute differences where we have expected values
    validComparisons = ~isnan(expectedValues) & ~isnan(trafficData);
    
    if sum(validComparisons) < 4
        outliers = false(size(trafficData));
        return;
    end
    
    % Use absolute difference: actual - expected
    absoluteDiffs = zeros(size(trafficData));
    absoluteDiffs(validComparisons) = trafficData(validComparisons) - expectedValues(validComparisons);
    
    % Define minimum thresholds based on data characteristics
    overallMedian = median(trafficData(~isnan(trafficData)));
    minAbsoluteDiff = max(5, overallMedian * 0.1);  % At least 5 or 10% of median, whichever is larger
    minExpectedValue = max(2, overallMedian * 0.05); % Only consider periods with expected > 2 or 5% of median
    
    fprintf('Outlier detection thresholds:\n');
    fprintf('  Minimum absolute difference: %.1f\n', minAbsoluteDiff);
    fprintf('  Minimum expected value: %.1f\n', minExpectedValue);
    
    % Only consider observations where:
    % 1. Expected value is above minimum threshold (to avoid low-count periods)
    % 2. Absolute difference is above minimum threshold (to avoid trivial differences)
    candidateOutliers = validComparisons & ...
                       expectedValues >= minExpectedValue & ...
                       abs(absoluteDiffs) >= minAbsoluteDiff;
    
    if sum(candidateOutliers) < 4
        fprintf('Insufficient candidate outliers for IQR analysis (need >= 4, found %d)\n', sum(candidateOutliers));
        outliers = false(size(trafficData));
        return;
    end
    
    % Apply IQR method to absolute differences of candidate outliers
    candidateDiffs = absoluteDiffs(candidateOutliers);
    Q1 = prctile(candidateDiffs, 25);
    Q3 = prctile(candidateDiffs, 75);
    IQR = Q3 - Q1;
    
    if IQR <= 0
        fprintf('Insufficient variation in candidate outliers for IQR analysis\n');
        outliers = false(size(trafficData));
        return;
    end
    
    lowerBound = Q1 - analysis.outlierDetection.threshold * IQR;
    upperBound = Q3 + analysis.outlierDetection.threshold * IQR;
    
    fprintf('IQR analysis on %d candidates: Q1=%.1f, Q3=%.1f, bounds=[%.1f, %.1f]\n', ...
        sum(candidateOutliers), Q1, Q3, lowerBound, upperBound);
    
    % Mark outliers (only among candidate outliers)
    outliers = false(size(trafficData));
    outliers(candidateOutliers) = (candidateDiffs < lowerBound) | (candidateDiffs > upperBound);
end

function outliers = detectContextualZScoreOutliers(inputTable, expectedValues, analysis)
    % Detect outliers using Z-score on absolute differences with minimum thresholds
    
    trafficData = inputTable.(analysis.modeString);
    
    % Calculate absolute differences where we have expected values
    validComparisons = ~isnan(expectedValues) & ~isnan(trafficData);
    
    if sum(validComparisons) < 3
        outliers = false(size(trafficData));
        return;
    end
    
    % Use absolute difference: actual - expected
    absoluteDiffs = zeros(size(trafficData));
    absoluteDiffs(validComparisons) = trafficData(validComparisons) - expectedValues(validComparisons);
    
    % Define minimum thresholds
    overallMedian = median(trafficData(~isnan(trafficData)));
    minAbsoluteDiff = max(5, overallMedian * 0.1);
    minExpectedValue = max(2, overallMedian * 0.05);
    
    fprintf('Outlier detection thresholds:\n');
    fprintf('  Minimum absolute difference: %.1f\n', minAbsoluteDiff);
    fprintf('  Minimum expected value: %.1f\n', minExpectedValue);
    
    % Only consider significant deviations from meaningful expected values
    candidateOutliers = validComparisons & ...
                       expectedValues >= minExpectedValue & ...
                       abs(absoluteDiffs) >= minAbsoluteDiff;
    
    if sum(candidateOutliers) < 3
        fprintf('Insufficient candidate outliers for Z-score analysis (need >= 3, found %d)\n', sum(candidateOutliers));
        outliers = false(size(trafficData));
        return;
    end
    
    % Apply Z-score method to absolute differences of candidate outliers
    candidateDiffs = absoluteDiffs(candidateOutliers);
    meanDiff = mean(candidateDiffs);
    stdDiff = std(candidateDiffs);
    
    if stdDiff == 0
        fprintf('No variation in candidate outlier differences\n');
        outliers = false(size(trafficData));
        return;
    end
    
    fprintf('Z-score analysis on %d candidates: mean=%.1f, std=%.1f\n', ...
        sum(candidateOutliers), meanDiff, stdDiff);
    
    % Mark outliers (only among candidate outliers)
    outliers = false(size(trafficData));
    zScores = abs(candidateDiffs - meanDiff) / stdDiff;
    outliers(candidateOutliers) = zScores > analysis.outlierDetection.threshold;
end

function reportContextualOutliers(inputTable, expectedValues, outlierFlags, analysis)
    % Report outliers with contextual information, focusing on meaningful deviations
    
    trafficData = inputTable.(analysis.modeString);
    outlierIndices = find(outlierFlags);
    
    if ~isempty(outlierIndices)
        fprintf('Found %d significant outlier(s):\n', length(outlierIndices));
        
        % Calculate absolute differences for reporting
        absoluteDiffs = nan(size(trafficData));
        validComparisons = ~isnan(expectedValues);
        absoluteDiffs(validComparisons) = trafficData(validComparisons) - expectedValues(validComparisons);
        
        % Sort outliers by magnitude of absolute difference
        outlierAbsDiffs = abs(absoluteDiffs(outlierIndices));
        [~, sortIdx] = sort(outlierAbsDiffs, 'descend', 'MissingPlacement', 'last');
        sortedIndices = outlierIndices(sortIdx);
        
        % Report all outliers (since we've filtered to meaningful ones)
        for i = 1:length(sortedIndices)
            idx = sortedIndices(i);
            actual = trafficData(idx);
            expected = expectedValues(idx);
            
            if ~isnan(expected)
                absDiff = actual - expected;
                if expected > 0
                    relDiff = absDiff / expected;
                    relDiffStr = sprintf('%.0f%% %s', abs(relDiff)*100, ternary(relDiff > 0, 'above', 'below'));
                else
                    relDiffStr = 'N/A (expected=0)';
                end
                
                dayType = inputTable.isWeekday(idx);
                dayTypeStr = ternary(dayType, 'WD', 'WE');
                
                fprintf('  %s (%s, %02d:00): actual=%d, expected=%.1f, diff=%+.1f (%s)\n', ...
                    datestr(inputTable.('Date and Time (Local)')(idx), 'dd-mmm-yyyy'), ...
                    dayTypeStr, inputTable.hourOfDay(idx), ...
                    actual, expected, absDiff, relDiffStr);
            else
                fprintf('  %s: actual=%d (no expected value available)\n', ...
                    datestr(inputTable.('Date and Time (Local)')(idx)), actual);
            end
        end
        
    else
        fprintf('No significant outliers detected.\n');
    end
end

function result = ternary(condition, trueValue, falseValue)
    % Simple ternary operator implementation
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end

function robustExpected = calculateRobustExpected(counts)
    % Calculate robust expected value using multiple methods and choose the most appropriate
    
    n = length(counts);
    
    if n < 2
        robustExpected = mean(counts);
        return;
    elseif n < 4
        % For small samples, use median
        robustExpected = median(counts);
        return;
    end
    
    % For larger samples, use multiple robust methods and choose based on data characteristics
    medianVal = median(counts);
    trimmedMean = trimmean(counts, 20);  % 20% trimmed mean (removes top/bottom 10%)
    
    % Use Median Absolute Deviation (MAD) to assess variability
    mad = median(abs(counts - medianVal));
    
    % Calculate coefficient of variation using robust statistics
    if medianVal > 0
        robustCV = mad / medianVal;
    else
        robustCV = inf;
    end
    
    % Choose method based on data characteristics
    if robustCV < 0.5
        % Low variability: trimmed mean is appropriate
        robustExpected = trimmedMean;
    else
        % High variability or potential outliers: use median
        robustExpected = medianVal;
    end
    
    % Additional check: if trimmed mean differs dramatically from median,
    % there are likely outliers, so prefer median
    if abs(trimmedMean - medianVal) > 2 * mad && mad > 0
        robustExpected = medianVal;
    end
end

function diagnosticMayWeekend18(inputTable, analysis, mayMask)
    % Diagnostic function to show what's happening with May weekend 18:00 values
    
    weekends = {'Saturday', 'Sunday'};
    may18WeekendMask = mayMask & ...
                       inputTable.hourOfDay == 18 & ...
                       ismember(inputTable.dayOfWeek, weekends);
    
    if any(may18WeekendMask)
        fprintf('\nDiagnostic - May 2025 Weekend 18:00 values:\n');
        may18Data = inputTable(may18WeekendMask, :);
        
        dates = may18Data.('Date and Time (Local)');
        counts = may18Data.(analysis.modeString);
        
        for i = 1:length(dates)
            fprintf('  %s: %d bikes\n', datestr(dates(i), 'dd-mmm-yyyy'), counts(i));
        end
        
        validCounts = counts(~isnan(counts));
        if ~isempty(validCounts)
            fprintf('  Raw statistics: mean=%.1f, median=%.1f, range=[%d, %d]\n', ...
                mean(validCounts), median(validCounts), min(validCounts), max(validCounts));
            
            robustExp = calculateRobustExpected(validCounts);
            fprintf('  Robust expected value: %.1f\n', robustExp);
        end
        fprintf('\n');
    end
end

function plotDayOfWeekPatterns(locationData, analysis, plots, style)
    % Plot day-of-week traffic patterns with monthly segregation
    
    locationNames = fieldnames(locationData);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate day-of-week patterns by month
        dayOfWeekData = calculateDayOfWeekPatterns(data, analysis);
        
        % Create monthly segregated plot
        plotMonthlyDayOfWeekPatterns(dayOfWeekData, analysis, style, locationInfo.name);
        
        % Create grand average plot
        plotGrandAverageDayOfWeek(dayOfWeekData, analysis, style, locationInfo.name);
    end
end

function dayOfWeekData = calculateDayOfWeekPatterns(locationDataStruct, analysis)
    % Calculate day-of-week patterns segregated by month
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Add month information if not present
    if ~ismember('monthStartDateTimes', data.Properties.VariableNames)
        data.monthStartDateTimes = dateshift(data.('Date and Time (Local)'), 'start', 'month');
    end
    
    % Define day order for consistent plotting
    dayOrder = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'};
    
    % Get unique months
    uniqueMonths = unique(data.monthStartDateTimes);
    
    % Initialize storage for monthly patterns
    dayOfWeekData = struct();
    dayOfWeekData.months = uniqueMonths;
    dayOfWeekData.monthlyPatterns = {};
    dayOfWeekData.dayNames = dayOrder;
    
    % Process each month
    for monthIdx = 1:length(uniqueMonths)
        currentMonth = uniqueMonths(monthIdx);
        monthData = data(data.monthStartDateTimes == currentMonth, :);
        
        if isempty(monthData)
            dayOfWeekData.monthlyPatterns{monthIdx} = nan(1, 7);
            continue;
        end
        
        % Calculate daily totals for this month
        monthData.DayOnly = dateshift(monthData.('Date and Time (Local)'), 'start', 'day');
        dailyTotals = groupsummary(monthData, 'DayOnly', 'sum', analysis.modeString);
        
        % Add day of week information to daily totals
        dailyTotals.dayOfWeek = string(day(dailyTotals.DayOnly, 'name'));
        
        % Calculate average for each day of week in this month
        monthPattern = nan(1, 7);
        for dayIdx = 1:7
            targetDay = dayOrder{dayIdx};
            dayRows = strcmp(dailyTotals.dayOfWeek, targetDay);
            
            if any(dayRows)
                sumColumnName = ['sum_' analysis.modeString];
                dayTotals = dailyTotals.(sumColumnName)(dayRows);
                monthPattern(dayIdx) = mean(dayTotals, 'omitnan');
            end
        end
        
        dayOfWeekData.monthlyPatterns{monthIdx} = monthPattern;
    end
    
    % Calculate grand averages across all months
    allPatterns = cell2mat(dayOfWeekData.monthlyPatterns');
    dayOfWeekData.grandAverage = mean(allPatterns, 1, 'omitnan');
end

function plotMonthlyDayOfWeekPatterns(dayOfWeekData, analysis, style, locationName)
    % Plot day-of-week patterns segregated by month
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    plotHandles = [];
    colorMap = lines(length(dayOfWeekData.months));
    
    % Plot each month's pattern
    for monthIdx = 1:length(dayOfWeekData.months)
        monthPattern = dayOfWeekData.monthlyPatterns{monthIdx};
        
        if ~all(isnan(monthPattern))
            monthTotal = sum(monthPattern, 'omitnan'); % Weekly total
            
            % Use subtle colors for individual months, emphasize most recent
            if monthIdx == length(dayOfWeekData.months)
                % Emphasize most recent month
                alpha = 1.0;
                lineWidth = style.plotLineWidth * 0.8;
                markerSize = 8;
            else
                % Subtle for other months
                alpha = 0.4;
                lineWidth = style.plotLineWidth * 0.5;
                markerSize = 6;
            end
            
            % h = plot(1:7, monthPattern, '-o', ...
            %     'LineWidth', lineWidth, ...
            %     'MarkerSize', markerSize, ...
            %     'Color', [colorMap(monthIdx, :) alpha], ...
            %     'MarkerFaceColor', [colorMap(monthIdx, :)],... % alpha], ...
            %     'DisplayName', sprintf('%s (weekly total ≈ %s)', ...
            %         datestr(dayOfWeekData.months(monthIdx), 'mmm yyyy'), ...
            %         num2sepstr(monthTotal, '%.0f')));
            h = plot(1:7, monthPattern, '-', ...
                'LineWidth', lineWidth, ...
                'MarkerSize', markerSize, ...
                'Color', [colorMap(monthIdx, :) alpha], ...
                'DisplayName', sprintf('%s (daily average ≈ %s)', ...
                    datestr(dayOfWeekData.months(monthIdx), 'mmm yyyy'), ...
                    num2sepstr(monthTotal./7, '%.0f')));
            plotHandles = [plotHandles, h];
        end
    end
    
    % Format plot
    formatDayOfWeekPlot(analysis, style, locationName, 'by Month', dayOfWeekData.dayNames);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function plotGrandAverageDayOfWeek(dayOfWeekData, analysis, style, locationName)
    % Plot the grand average day-of-week pattern
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    plotHandles = [];
    
    % Plot grand average pattern
    if ~all(isnan(dayOfWeekData.grandAverage))
        weeklyTotal = sum(dayOfWeekData.grandAverage, 'omitnan') * 7;
        
        h = plot(1:7, dayOfWeekData.grandAverage, '-o', ...
            'LineWidth', style.plotLineWidth, ...
            'MarkerSize', 10, ...
            'Color', [0 0 1], ...
            'MarkerFaceColor', [0 0 1], ...
            'DisplayName', sprintf('Average (weekly total ≈ %s)', ...
                num2sepstr(weeklyTotal, '%.0f')));
        plotHandles = [plotHandles, h];
        
        % Add error bars or confidence intervals if desired
        % (would require storing standard deviations in calculateDayOfWeekPatterns)
    end
    
    % Format plot
    formatDayOfWeekPlot(analysis, style, locationName, 'Average', dayOfWeekData.dayNames);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function formatDayOfWeekPlot(analysis, style, locationName, plotType, dayNames)
    % Format the day-of-week pattern plot
    
    ylabel('Average Daily Count', 'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    xlabel('Day of Week', 'FontSize', style.labelFontSize);
    
    title([plotType ' Day-of-Week ' analysis.modeDisplayString ' (' locationName ')'], ...
        'FontSize', style.titleFontSize);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    % Set x-axis to show day names
    xlim([0.5 7.5]);
    xticks(1:7);
    xticklabels(dayNames);
    xtickangle(45);
end

function plotTemperatureScatter(locationData, weatherData, analysis, plots, style)
    % Plot scatter plot of daily counts versus air temperature for all locations
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot each location with different colors
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate daily totals and match with weather data
        [dailyCounts, dailyTemperatures, validDates] = prepareDailyTemperatureData(data, weatherData, analysis);
        
        if ~isempty(dailyCounts)
            % Create scatter plot
            h = scatter(dailyTemperatures, dailyCounts, 50, ...
                'MarkerEdgeColor', locationInfo.plotColor, ...
                'MarkerFaceColor', locationInfo.plotColor, ...
                'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeAlpha', 0.8, ...
                'DisplayName', sprintf('%s (%d days)', locationInfo.name, length(dailyCounts)));
            plotHandles = [plotHandles, h];
            
            % Add piecewise linear trend line with breakpoint at 0°C
            if length(dailyTemperatures) > 5
                [trendParams, trendStats] = fitPiecewiseLinear(dailyTemperatures, dailyCounts);
                
                if ~isempty(trendParams)
                    % Plot the piecewise linear trend
                    tempRange = linspace(min(dailyTemperatures), max(dailyTemperatures), 100);
                    trendLine = evaluatePiecewiseLinear(tempRange, trendParams);
                    
                    plot(tempRange, trendLine, '--', ...
                        'Color', locationInfo.plotColor, ...
                        'LineWidth', 2, ...
                        'HandleVisibility', 'off');  % Don't show in legend
                    
                    % Store trend stats for later display
                    locationInfo.trendStats = trendStats;
                end
            end
        end
    end
    
    % Format plot
    formatTemperatureScatterPlot(analysis, style);
    
    % Add correlation statistics
    addCorrelationStats(locationData, weatherData, analysis, style);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function [dailyCounts, dailyTemperatures, validDates] = prepareDailyTemperatureData(locationDataStruct, weatherData, analysis)
    % Prepare matched daily counts and temperature data
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Calculate daily totals
    data.DayOnly = dateshift(data.('Date and Time (Local)'), 'start', 'day');
    groupedData = groupsummary(data, 'DayOnly', 'sum', analysis.modeString);
    
    % Get the sum column name
    sumColumnName = ['sum_' analysis.modeString];
    dailyData = table(groupedData.DayOnly, groupedData.(sumColumnName), ...
        'VariableNames', {'Date', 'Count'});
    
    % Match with weather data
    [~, ia, ib] = intersect(dailyData.Date, weatherData.dates);
    
    if isempty(ia)
        % No matching dates
        dailyCounts = [];
        dailyTemperatures = [];
        validDates = [];
        return;
    end
    
    % Extract matched data
    dailyCounts = dailyData.Count(ia);
    dailyTemperatures = weatherData.temperature(ib);
    %dailyTemperatures = weatherData.feelslike(ib);
    validDates = dailyData.Date(ia);
    
    % Remove any NaN values
    validIdx = ~isnan(dailyCounts) & ~isnan(dailyTemperatures);
    dailyCounts = dailyCounts(validIdx);
    dailyTemperatures = dailyTemperatures(validIdx);
    validDates = validDates(validIdx);
end

function formatTemperatureScatterPlot(analysis, style)
    % Format the temperature scatter plot
    
    xlabel('Air Temperature (°C)', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel(['Daily ' analysis.modeDisplayString], 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(['Daily ' analysis.modeDisplayString ' vs Air Temperature'], ...
        'FontSize', style.titleFontSize);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    % Ensure y-axis starts at 0
    ylim([0 max(ylim) * 1.05]);
    
    % Add temperature reference lines (optional)
    xLimits = xlim;
    
    % Add vertical lines for key temperatures
    line([0 0], ylim, 'Color', [0.7 0.7 1], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    line([10 10], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    line([20 20], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Add temperature labels
    text(0, max(ylim) * 0.95, '0°C', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
    text(10, max(ylim) * 0.95, '10°C', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
    text(20, max(ylim) * 0.95, '20°C', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
end

function addCorrelationStats(locationData, weatherData, analysis, style)
    % Add piecewise linear model statistics as text annotation
    
    locationNames = fieldnames(locationData);
    statsText = {};
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Get daily data for analysis
        [dailyCounts, dailyTemperatures, ~] = prepareDailyTemperatureData(data, weatherData, analysis);
        
        if length(dailyCounts) > 5
            % Fit piecewise linear model
            [trendParams, trendStats] = fitPiecewiseLinear(dailyTemperatures, dailyCounts);
            
            if ~isempty(trendParams)
                % Format statistics for display
                belowFreezing = sprintf('T<0°C: slope=%.2f', trendStats.slopeBelow);
                aboveFreezing = sprintf('T≥0°C: slope=%.2f', trendStats.slopeAbove);
                rSquared = sprintf('R²=%.3f', trendStats.rSquared);
                
                statsText{end+1} = sprintf('%s: %s, %s, %s', ...
                    locationInfo.name, belowFreezing, aboveFreezing, rSquared);
            else
                % Fallback to simple correlation if piecewise fit fails
                corrCoeff = corr(dailyTemperatures, dailyCounts, 'type', 'Pearson');
                statsText{end+1} = sprintf('%s: r = %.3f (simple)', locationInfo.name, corrCoeff);
            end
        end
    end
    
    % Display statistics
    if ~isempty(statsText)
        % Position text box in upper left
        xLimits = xlim;
        yLimits = ylim;
        
        textX = xLimits(1) + 0.05 * (xLimits(2) - xLimits(1));
        textY = yLimits(2) - 0.1 * (yLimits(2) - yLimits(1));
        
        % Create text box
        textStr = strjoin(statsText, '\n');
        text(textX, textY, textStr, ...
            'FontSize', style.axisFontSize * 0.9, ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.9], ...
            'EdgeColor', [0.7 0.7 0.7], ...
            'Margin', 5);
    end
end

function plotTemperatureScatterSeasonal(locationData, weatherData, analysis, plots, style)
    % Alternative version: Plot with seasonal color coding
    % Call this instead of plotTemperatureScatter if you want seasonal analysis
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Define seasons
    seasonColors = [
        0.2, 0.7, 0.2;    % Spring (green)
        1.0, 0.6, 0.0;    % Summer (orange)
        0.8, 0.4, 0.0;    % Fall (brown)
        0.2, 0.4, 0.8     % Winter (blue)
    ];
    seasonNames = {'Spring', 'Summer', 'Fall', 'Winter'};
    
    % Plot each location
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Get daily data
        [dailyCounts, dailyTemperatures, validDates] = prepareDailyTemperatureData(data, weatherData, analysis);
        
        if ~isempty(dailyCounts)
            % Assign seasons
            seasons = getSeason(validDates);
            
            % Plot each season separately
            for seasonIdx = 1:4
                seasonMask = seasons == seasonIdx;
                if any(seasonMask)
                    % Adjust marker style for location
                    if i == 1
                        markerStyle = 'o';  % Circle for first location
                    else
                        markerStyle = 's';  % Square for second location
                    end
                    
                    h = scatter(dailyTemperatures(seasonMask), dailyCounts(seasonMask), 50, ...
                        'Marker', markerStyle, ...
                        'MarkerEdgeColor', seasonColors(seasonIdx, :), ...
                        'MarkerFaceColor', seasonColors(seasonIdx, :), ...
                        'MarkerFaceAlpha', 0.6, ...
                        'MarkerEdgeAlpha', 0.8, ...
                        'DisplayName', sprintf('%s %s', locationInfo.name, seasonNames{seasonIdx}));
                    
                    if seasonIdx == 1  % Only add to legend for first season
                        plotHandles = [plotHandles, h];
                    end
                end
            end
        end
    end
    
    % Format plot
    formatTemperatureScatterPlot(analysis, style);
    title(['Daily ' analysis.modeDisplayString ' vs Air Temperature (by Season)'], ...
        'FontSize', style.titleFontSize);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function seasons = getSeason(dates)
    % Assign seasons based on month
    % 1 = Spring (Mar-May), 2 = Summer (Jun-Aug), 3 = Fall (Sep-Nov), 4 = Winter (Dec-Feb)
    
    months = month(dates);
    seasons = zeros(size(months));
    
    seasons(ismember(months, [3, 4, 5])) = 1;  % Spring
    seasons(ismember(months, [6, 7, 8])) = 2;  % Summer
    seasons(ismember(months, [9, 10, 11])) = 3; % Fall
    seasons(ismember(months, [12, 1, 2])) = 4;  % Winter
end

function [params, stats] = fitPiecewiseLinear(temperatures, counts)
    % Fit piecewise linear model with breakpoint at 0°C
    % Model: y = a1*x + b1 for x < 0, y = a2*x + b2 for x >= 0
    % with continuity constraint at x = 0: b1 = b2
    
    try
        % Separate data into below and above freezing
        belowFreezing = temperatures < 0;
        aboveFreezing = temperatures >= 0;
        
        % Need at least 2 points in each segment for fitting
        if sum(belowFreezing) < 2 || sum(aboveFreezing) < 2
            params = [];
            stats = [];
            return;
        end
        
        % Fit two separate linear models
        tempBelow = temperatures(belowFreezing);
        countsBelow = counts(belowFreezing);
        tempAbove = temperatures(aboveFreezing);
        countsAbove = counts(aboveFreezing);
        
        % Linear fit for below freezing: y = a1*x + b1
        pBelow = polyfit(tempBelow, countsBelow, 1);
        a1 = pBelow(1);  % slope
        b1 = pBelow(2);  % intercept
        
        % Linear fit for above freezing: y = a2*x + b2
        pAbove = polyfit(tempAbove, countsAbove, 1);
        a2 = pAbove(1);  % slope
        b2 = pAbove(2);  % intercept
        
        % Store parameters
        params = struct();
        params.slopeBelow = a1;
        params.interceptBelow = b1;
        params.slopeAbove = a2;
        params.interceptAbove = b2;
        params.breakpoint = 0;
        
        % Calculate goodness of fit
        predictedCounts = evaluatePiecewiseLinear(temperatures, params);
        
        % Calculate R-squared
        SSres = sum((counts - predictedCounts).^2);
        SStot = sum((counts - mean(counts)).^2);
        rSquared = 1 - SSres/SStot;
        
        % Calculate individual segment statistics
        predictedBelow = a1 * tempBelow + b1;
        predictedAbove = a2 * tempAbove + b2;
        
        SSresBelow = sum((countsBelow - predictedBelow).^2);
        SStotBelow = sum((countsBelow - mean(countsBelow)).^2);
        rSquaredBelow = 1 - SSresBelow/SStotBelow;
        
        SSresAbove = sum((countsAbove - predictedAbove).^2);
        SStotAbove = sum((countsAbove - mean(countsAbove)).^2);
        rSquaredAbove = 1 - SSresAbove/SStotAbove;
        
        % Store statistics
        stats = struct();
        stats.slopeBelow = a1;
        stats.slopeAbove = a2;
        stats.interceptBelow = b1;
        stats.interceptAbove = b2;
        stats.rSquared = rSquared;
        stats.rSquaredBelow = rSquaredBelow;
        stats.rSquaredAbove = rSquaredAbove;
        stats.nPointsBelow = sum(belowFreezing);
        stats.nPointsAbove = sum(aboveFreezing);
        
        % Check if the breakpoint makes sense (significant difference in slopes)
        if abs(a1 - a2) < 0.1  % Slopes are very similar
            % Maybe a simple linear model would be better
            fprintf('Warning: Piecewise slopes are very similar (%.2f vs %.2f)\n', a1, a2);
        end
        
    catch ME
        fprintf('Error fitting piecewise linear model: %s\n', ME.message);
        params = [];
        stats = [];
    end
end

function predictedValues = evaluatePiecewiseLinear(temperatures, params)
    % Evaluate piecewise linear model at given temperatures
    
    predictedValues = zeros(size(temperatures));
    
    belowFreezing = temperatures < params.breakpoint;
    aboveFreezing = temperatures >= params.breakpoint;
    
    % Below freezing: y = a1*x + b1
    predictedValues(belowFreezing) = params.slopeBelow * temperatures(belowFreezing) + params.interceptBelow;
    
    % Above freezing: y = a2*x + b2
    predictedValues(aboveFreezing) = params.slopeAbove * temperatures(aboveFreezing) + params.interceptAbove;
end

function plotModalityPieCharts(locationData, analysis, style)
    % Generate pie charts showing breakdown of counts by modality for each location
    
    locationNames = fieldnames(locationData);
    
    % Define the modalities to include in the pie charts
    modalityInfo = getModalityInfo();
    
    % Create subplot layout - one pie chart per location
    numLocations = length(locationNames);
    if numLocations == 1
        figure('Position', [408 126 600 600]);
    elseif numLocations == 2
        figure('Position', [408 126 1200 600]);
    else
        figure('Position', [408 126 1200 800]);
    end
    
    for i = 1:numLocations
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate totals for each modality
        [modalityCounts, modalityLabels, modalityColors] = calculateModalityTotals(data, modalityInfo, analysis);
        
        % Create subplot
        if numLocations <= 2
            subplot(1, numLocations, i);
        else
            subplot(2, ceil(numLocations/2), i);
        end
        
        % Generate pie chart
        if ~isempty(modalityCounts) && sum(modalityCounts) > 0
            createModalityPieChart(modalityCounts, modalityLabels, modalityColors, locationInfo.name, analysis, style);
        else
            % Handle case with no data
            text(0.5, 0.5, 'No Data Available', 'HorizontalAlignment', 'center', ...
                'FontSize', style.labelFontSize, 'FontWeight', 'bold');
            title(locationInfo.name, 'FontSize', style.titleFontSize);
            axis off;
        end
    end
    
    % Add overall figure title
    sgtitle(sprintf('Traffic Modality Breakdown (%s to %s)', ...
        datestr(analysis.startTime, 'mmm dd yyyy'), datestr(analysis.endTime, 'mmm dd yyyy')), ...
        'FontSize', style.titleFontSize + 2, 'FontWeight', 'bold');
end

function modalityInfo = getModalityInfo()
    % Define modalities and their display properties
    
    modalityInfo = struct();
    
    % Define each modality with column name, display name, and color
    modalityInfo.modalities = {
        struct('columnName', 'Bike Total', 'displayName', 'Bikes', 'color', [0.2, 0.6, 0.2]);
        struct('columnName', 'Pedestrian Total', 'displayName', 'Pedestrians', 'color', [0.8, 0.4, 0.2]);
        struct('columnName', 'Car Total', 'displayName', 'Cars', 'color', [0.6, 0.6, 0.6]);
        struct('columnName', 'Large vehicle Total', 'displayName', 'Trucks', 'color', [0.4, 0.4, 0.4]);
    };
end

function [counts, labels, colors] = calculateModalityTotals(locationDataStruct, modalityInfo, analysis)
    % Calculate total counts for each modality over the analysis period
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Initialize outputs
    counts = [];
    labels = {};
    colors = [];
    
    % Process each modality
    for i = 1:length(modalityInfo.modalities)
        modality = modalityInfo.modalities{i};
        
        % Check if this modality column exists in the data
        if ismember(modality.columnName, data.Properties.VariableNames)
            % Calculate total counts for this modality
            modalityData = data.(modality.columnName);
            totalCount = sum(modalityData, 'omitnan');
            
            % Only include modalities with non-zero counts
            if totalCount > 0
                counts = [counts; totalCount];
                labels{end+1} = modality.displayName;
                colors = [colors; modality.color];
            end
        end
    end
end

function createModalityPieChart(counts, labels, colors, locationName, analysis, style)
    % Create a pie chart with custom formatting
    
    % Calculate percentages
    totalCount = sum(counts);
    percentages = (counts / totalCount) * 100;
    
    % Create labels with both count and percentage
    pieLabels = cell(size(labels));
    for i = 1:length(labels)
        pieLabels{i} = sprintf('%s\n%s (%.1f%%)', ...
            labels{i}, ...
            num2sepstr(counts(i), '%.0f'), ...
            percentages(i));
    end
    
    % Create pie chart
    p = pie(counts);
    
    % Customize pie chart appearance
    for i = 1:2:length(p)  % Every other element is a patch (the pie slice)
        patchIndex = (i+1)/2;
        if patchIndex <= size(colors, 1)
            set(p(i), 'FaceColor', colors(patchIndex, :));
            set(p(i), 'EdgeColor', 'white');
            set(p(i), 'LineWidth', 2);
        end
    end
    
    % Customize text labels
    for i = 2:2:length(p)  % Every other element starting from 2 is text
        set(p(i), 'FontSize', style.axisFontSize * 0.9);
        set(p(i), 'FontWeight', 'bold');
        set(p(i), 'Color', [0.2 0.2 0.2]);
    end
    
    % Set custom labels
    for i = 1:length(pieLabels)
        if 2*i <= length(p)
            set(p(2*i), 'String', pieLabels{i});
        end
    end
    
    % Add title with total count
    title(sprintf('%s\nTotal: %s', locationName, num2sepstr(totalCount, '%.0f')), ...
        'FontSize', style.titleFontSize * 0.9, 'FontWeight', 'bold');
    
    % Ensure equal aspect ratio for circular pie
    axis equal;
end

function plotModalityPieChartsComparative(locationData, analysis, style)
    % Alternative version: Create a single comparative pie chart showing all locations
    
    figure('Position', [408 126 800 600]);
    
    locationNames = fieldnames(locationData);
    modalityInfo = getModalityInfo();
    
    % Combine data from all locations
    combinedCounts = [];
    combinedLabels = {};
    combinedColors = [];
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate totals for each modality at this location
        [modalityCounts, modalityLabels, modalityColors] = calculateModalityTotals(data, modalityInfo, analysis);
        
        % Add location prefix to labels
        for j = 1:length(modalityLabels)
            combinedLabels{end+1} = sprintf('%s - %s', extractLocationShortName(locationInfo.name), modalityLabels{j});
            combinedCounts = [combinedCounts; modalityCounts(j)];
            combinedColors = [combinedColors; modalityColors(j, :)];
        end
    end
    
    % Create combined pie chart
    if ~isempty(combinedCounts) && sum(combinedCounts) > 0
        createModalityPieChart(combinedCounts, combinedLabels, combinedColors, 'All Locations', analysis, style);
    else
        text(0.5, 0.5, 'No Data Available', 'HorizontalAlignment', 'center', ...
            'FontSize', style.labelFontSize, 'FontWeight', 'bold');
        axis off;
    end
    
    % Add overall title
    sgtitle(sprintf('Combined Traffic Modality Breakdown (%s to %s)', ...
        datestr(analysis.startTime, 'mmm yyyy'), datestr(analysis.endTime, 'mmm yyyy')), ...
        'FontSize', style.titleFontSize + 2, 'FontWeight', 'bold');
end

function shortName = extractLocationShortName(fullName)
    % Extract a short name from the full location name for labeling
    
    % Look for common patterns to shorten
    if contains(fullName, '@')
        parts = split(fullName, '@');
        if length(parts) >= 2
            shortName = strtrim(parts{2});
        else
            shortName = strtrim(parts{1});
        end
    elseif contains(fullName, 'rue de Terrebonne')
        % Extract the cross street
        if contains(fullName, 'King Edward')
            shortName = 'King Edward';
        elseif contains(fullName, 'Draper')
            shortName = 'Draper';
        else
            shortName = fullName;
        end
    else
        % Use first few words
        words = split(fullName);
        if length(words) >= 2
            shortName = sprintf('%s %s', words{1}, words{2});
        else
            shortName = fullName;
        end
    end
    
    % Limit length
    if length(shortName) > 15
        shortName = shortName(1:15);
    end
end

function plotModalityBarChart(locationData, analysis, style)
    % Alternative visualization: Stacked bar chart instead of pie charts
    
    figure('Position', [408 126 1000 600]);
    
    locationNames = fieldnames(locationData);
    modalityInfo = getModalityInfo();
    
    % Prepare data matrix
    dataMatrix = [];
    locationLabels = {};
    modalityLabels = {};
    modalityColors = [];
    
    % Get modality labels from first location that has data
    for i = 1:length(locationNames)
        [~, labels, colors] = calculateModalityTotals(locationData.(locationNames{i}), modalityInfo, analysis);
        if ~isempty(labels)
            modalityLabels = labels;
            modalityColors = colors;
            break;
        end
    end
    
    % Build data matrix
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        [counts, labels, ~] = calculateModalityTotals(data, modalityInfo, analysis);
        
        % Align counts with standard modality order
        alignedCounts = zeros(1, length(modalityLabels));
        for j = 1:length(labels)
            idx = find(strcmp(modalityLabels, labels{j}));
            if ~isempty(idx)
                alignedCounts(idx) = counts(j);
            end
        end
        
        dataMatrix = [dataMatrix; alignedCounts];
        locationLabels{end+1} = extractLocationShortName(locationInfo.name);
    end
    
    % Create stacked bar chart
    if ~isempty(dataMatrix)
        b = bar(dataMatrix, 'stacked');
        
        % Apply colors
        for i = 1:length(b)
            if i <= size(modalityColors, 1)
                set(b(i), 'FaceColor', modalityColors(i, :));
                set(b(i), 'EdgeColor', 'white');
                set(b(i), 'LineWidth', 1);
            end
        end
        
        % Format axes
        set(gca, 'XTickLabel', locationLabels);
        ylabel('Total Count', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
        xlabel('Location', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
        
        % Add legend
        legend(modalityLabels, 'Location', 'northeast', 'FontSize', style.legendFontSize);
        
        % Format y-axis with separators
        ytick_positions = yticks;
        ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
        yticklabels(ytick_labels);
        
        set(gca, 'FontSize', style.axisFontSize);
        grid on;
        set(gca, 'Color', style.axisBackgroundColor);
        
        title(sprintf('Traffic Modality Comparison (%s to %s)', ...
            datestr(analysis.startTime, 'mmm dd yyyy'), datestr(analysis.endTime, 'mmm dd yyyy')), ...
            'FontSize', style.titleFontSize, 'FontWeight', 'bold');
    end
end

function plotCombinedHourlyRaw(locationData, weatherData, analysis, plots, style)
    % Plot raw hourly counts as scatter points for all locations
    
    figure('Position', [408 126 1400 921]);  % Wider figure for dense time data
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    % if plots.showWeather
    %     yyaxis right
    %     weatherHandles = plotWeatherData(weatherData, plots, style);
    % end
    
    % Plot traffic data on left axis
    %yyaxis left
    hold on
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Extract hourly data
        hourlyData = extractHourlyData(data, analysis);
        
        if ~isempty(hourlyData.dateTimes)
            % Plot raw counts as points
            if plots.showRawCounts
                
                % Square root scaling makes differences more visible for count data
                minSize = 10;
                maxSize = 75;
                counts = hourlyData.rawCounts;

                % Apply square root transformation
                sqrtCounts = sqrt(counts);
                minSqrt = min(sqrtCounts);
                maxSqrt = max(sqrtCounts);

                if maxSqrt > minSqrt
                    scaledSizes = minSize + (maxSize - minSize) * (sqrtCounts - minSqrt) / (maxSqrt - minSqrt);
                else
                    scaledSizes = repmat(minSize, size(counts));
                end

                h1 = scatter(hourlyData.dateTimes, hourlyData.rawCounts, scaledSizes, ...
                    'MarkerFaceColor', locationInfo.plotColor, ...
                    'MarkerEdgeColor', 'none', ...
                    'MarkerFaceAlpha', 1.0, ...  % 50% transparency
                    'DisplayName', sprintf('%s (%s hours counted;  total counts:  %s; average count:  %s per hour;  range: %s-%s per hour)', ...
                    locationInfo.name, ...
                    num2sepstr(length(hourlyData.rawCounts),'%.0f'), ...
                    num2sepstr(sum(hourlyData.rawCounts), '%.0f'), ...
                    num2sepstr(mean(hourlyData.rawCounts), '%.1f'), ...
                    num2sepstr(min(hourlyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(hourlyData.rawCounts), '%.0f')));

                    ytvals = [0 10:10:150];
                    yticks(ytvals);
                
                plotHandles = [plotHandles, h1];
            end

            % Plot adjusted counts if enabled (slightly offset for visibility)
            if plots.showAdjustedCounts
                % Add small time offset to avoid overlapping points
                offsetHours = (i-1) * 0.1;  % 6-minute offset per location
                offsetTimes = hourlyData.dateTimes + hours(offsetHours);
                
                h2 = plot(offsetTimes, hourlyData.adjustedCounts, '.', ...
                    'MarkerSize', 2, ...
                    'Color', locationInfo.plotColor * 0.7, ...
                    'DisplayName', sprintf('%s Adjusted (n=%d, range: %s-%s)', ...
                        locationInfo.name, ...
                        length(hourlyData.adjustedCounts), ...
                        num2sepstr(min(hourlyData.adjustedCounts), '%.0f'), ...
                        num2sepstr(max(hourlyData.adjustedCounts), '%.0f')));
                plotHandles = [plotHandles, h2];
            end
        end
    end
    
    % Format plot
    formatCombinedHourlyRawPlot(analysis, plots, style, weatherData);
    
    % if plots.showWeather
    %     plotHandles = [plotHandles, weatherHandles];
    % end

    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function hourlyData = extractHourlyData(locationDataStruct, analysis)
    % Extract hourly raw data for plotting
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Get the relevant columns
    dateTimes = data.('Date and Time (Local)');
    rawCounts = data.(analysis.modeString);
    
    % Get adjusted counts if available
    if ismember('AdjustedCountsUptimeDaylight', data.Properties.VariableNames)
        adjustedCounts = data.AdjustedCountsUptimeDaylight;
    else
        adjustedCounts = rawCounts;  % Fallback to raw if adjusted not available
    end
    
    % Remove NaN values
    validIdx = ~isnan(rawCounts);
    
    hourlyData = struct();
    hourlyData.dateTimes = dateTimes(validIdx);
    hourlyData.rawCounts = rawCounts(validIdx);
    hourlyData.adjustedCounts = adjustedCounts(validIdx);
end

function formatCombinedHourlyRawPlot(analysis, plots, style, weatherData)
    % Format the hourly raw count plot
    
    % Left axis formatting
    %yyaxis left
    ylabel(['Hourly ' analysis.modeDisplayString], ...
        'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    % if plots.showWeather
    %     ax.YAxis(2).Color = [0 0.4471 0.7412];
    % end
    
    % Title and formatting
    title(['Hourly ' analysis.modeDisplayString ' (Raw Telraam Data)'], ...
        'FontSize', style.titleFontSize);
    
    subtitle(sprintf('%s to %s', ...
        datestr(analysis.startTime, 'dd-mmm-yyyy'), ...
        datestr(analysis.endTime, 'dd-mmm-yyyy')), ...
        'FontSize', style.axisFontSize, 'Color', 0.3.*[1 1 1]);

    xlabel('Date and Time', 'FontSize', style.labelFontSize);
    set(gca, 'Color', style.axisBackgroundColor);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([-1 max(ylim) * 1.1]);
    
    % Set appropriate x-axis limits and ticks
    if ~isempty(weatherData)
        xlim([weatherData.dates(1) weatherData.dates(end)]);
    end
    
    % Optimize x-axis tick spacing for dense time data
    xLimits = xlim;
    dateRange = days(xLimits(2) - xLimits(1));
    
    if dateRange <= 7
        % For week or less: show every day
        xtickformat('dd-MMM HH:mm');
        xtickangle(45);
    elseif dateRange <= 31
        % For month or less: show every few days
        ax.XAxis.TickLabelFormat = 'dd-MMM';
        xtickangle(45);
    elseif dateRange <= 90
        % For season or less: show weeks
        ax.XAxis.TickLabelFormat = 'dd-MMM';
        xtickangle(45);
    else
        % For longer periods: show months
        ax.XAxis.TickLabelFormat = 'MMM yyyy';
        xtickangle(45);
    end
end

function plotMultiModalHourlyRaw(locationData, weatherData, analysis, plots, style, multiModal)
    % Plot hourly raw counts for multiple modes at a single location
    
    figure('Position', [408 126 1400 921]);
    
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
        
        % Extract hourly data for this mode
        hourlyData = extractHourlyData(data, tempAnalysis);
        
        if ~isempty(hourlyData.dateTimes)
            % Add small time offset to separate overlapping points
            offsetHours = (i-1) * 0.05;  % 3-minute offset per mode
            offsetTimes = hourlyData.dateTimes + hours(offsetHours);
            
            h1 = plot(offsetTimes, hourlyData.rawCounts, '.', ...
                'MarkerSize', 4, ...
                'Color', currentColor, ...
                'DisplayName', sprintf('%s (n=%d, range: %s-%s)', ...
                    currentModeDisplay, ...
                    length(hourlyData.rawCounts), ...
                    num2sepstr(min(hourlyData.rawCounts), '%.0f'), ...
                    num2sepstr(max(hourlyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
    end
    
    % Format plot
    formatMultiModalHourlyRawPlot(multiModal, plots, style, weatherData);
    
    if multiModal.plotWeather && plots.showWeather
        plotHandles = [plotHandles, weatherHandles];
    end
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function formatMultiModalHourlyRawPlot(multiModal, plots, style, weatherData)
    % Format the multi-modal hourly raw plot
    
    % Left axis formatting
    yyaxis left
    ylabel('Hourly Counts', ...
        'FontSize', style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    if multiModal.plotWeather && plots.showWeather
        ax.YAxis(2).Color = [0 0.4471 0.7412];
    end
    
    % Title and formatting
    title(['Hourly Traffic by Mode (' multiModal.location ') - Raw Time Series'], ...
        'FontSize', style.titleFontSize);
    
    xlabel('Date and Time', 'FontSize', style.labelFontSize);
    set(gca, 'Color', style.axisBackgroundColor);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    ylim([0 max(ylim) * 1.1]);
    
    % Set x-axis limits and formatting
    if ~isempty(weatherData)
        xlim([weatherData.dates(1) weatherData.dates(end)]);
    end
    
    % Optimize x-axis tick spacing
    xLimits = xlim;
    dateRange = days(xLimits(2) - xLimits(1));
    
    if dateRange <= 7
        xtickformat('dd-MMM HH:mm');
        xtickangle(45);
    elseif dateRange <= 31
        ax.XAxis.TickLabelFormat = 'dd-MMM';
        xtickangle(45);
    else
        ax.XAxis.TickLabelFormat = 'MMM yyyy';
        xtickangle(45);
    end
end

function filteredLocationData = filterNightData(locationData, weatherData)
    % Filter out nighttime data points where bike/pedestrian detection is not possible
    % 
    % This function removes observations where:
    % - Both bike and pedestrian counts are zero AND
    % - Time is before sunrise OR after sunset AND  
    % - Night Total > 0 (indicating nighttime conditions)
    %
    % This provides a more accurate representation of detection capabilities
    % during periods when the computer vision system can actually detect bikes/pedestrians
    
    fprintf('Filtering nighttime data where bike/pedestrian detection is not possible...\n');
    
    filteredLocationData = struct();
    locationNames = fieldnames(locationData);
    
    % Process each location
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        fprintf('Processing location: %s\n', locationInfo.name);
        
        % Apply night filtering to this location's data
        filteredData = applyNightFilter(data.data, weatherData);
        
        % Store filtered data back in the structure
        filteredLocationData.(locationName) = struct();
        filteredLocationData.(locationName).data = filteredData;
        filteredLocationData.(locationName).locationInfo = locationInfo;
        
        % Report filtering statistics
        originalCount = height(data.data);
        filteredCount = height(filteredData);
        removedCount = originalCount - filteredCount;
        
        fprintf('  Original observations: %s\n', num2sepstr(originalCount, '%.0f'));
        fprintf('  Filtered observations: %s\n', num2sepstr(filteredCount, '%.0f'));
        fprintf('  Removed (nighttime zeros): %s (%.1f%%)\n', ...
            num2sepstr(removedCount, '%.0f'), 100 * removedCount / originalCount);
    end
    
    fprintf('Night data filtering complete.\n\n');
end

function filteredData = applyNightFilter(data, weatherData)
    % Apply the night filtering criteria to a single location's data
    
    % Check required dimensions and columns exist
    requiredDimensions = {'Date and Time (Local)'};
    missingDimensions = {};
    requiredColumns = {'Bike Total', 'Pedestrian Total', 'Night Total'};
    missingColumns = {};
    
    for i = 1:length(requiredColumns)
        if ~ismember(requiredColumns{i}, data.Properties.VariableNames)
            missingColumns{end+1} = requiredColumns{i};
        end
    end

    for i = 1:length(requiredDimensions)
        if ~ismember(requiredDimensions{i}, data.Properties.DimensionNames)
            missingDimensions{end+1} = requiredDimensions{i};
        end
    end
    
    if ~isempty(missingColumns)
        warning('Missing required columns for night filtering: %s. Returning original data.', ...
            strjoin(missingColumns, ', '));
        filteredData = data;
        return;
    end

    if ~isempty(missingDimensions)
        warning('Missing required dimensions for night filtering: %s. Returning original data.', ...
            strjoin(missingDimensions, ', '));
        filteredData = data;
        return;
    end
    
    % Extract relevant columns
    dateTimes = data.('Date and Time (Local)');
    bikeCounts = data.('Bike Total');
    pedestrianCounts = data.('Pedestrian Total');
    nightTotal = data.('Night Total');
    
    % Match observation times with sunrise/sunset data
    [sunriseMatched, sunsetMatched] = matchSunriseSunsetTimes(dateTimes, weatherData);
    
    % Identify observations to potentially remove
    % Criteria: bikes = 0 AND pedestrians = 0 AND (before sunrise OR after sunset) AND night total > 0
    
    zeroDetectionCounts = (bikeCounts == 0) & (pedestrianCounts == 0);
    nighttimeConditions = (dateTimes < sunriseMatched) | (dateTimes > sunsetMatched);
    nightSystemActive = nightTotal > 0;
    
    % Combine all criteria - these are the rows to REMOVE
    rowsToRemove = zeroDetectionCounts & (nighttimeConditions | nightSystemActive);
    
    % Apply filter (keep rows that DON'T meet removal criteria)
    filteredData = data(~rowsToRemove, :);
    
    % Report detailed statistics
    totalObs = length(dateTimes);
    zeroCountObs = sum(zeroDetectionCounts);
    nighttimeObs = sum(nighttimeConditions);
    nightActiveObs = sum(nightSystemActive);
    removedObs = sum(rowsToRemove);
    
    fprintf('    Detailed filtering statistics:\n');
    fprintf('      Zero bike & pedestrian counts: %s (%.1f%%)\n', ...
        num2sepstr(zeroCountObs, '%.0f'), 100 * zeroCountObs / totalObs);
    fprintf('      Nighttime observations: %s (%.1f%%)\n', ...
        num2sepstr(nighttimeObs, '%.0f'), 100 * nighttimeObs / totalObs);
    fprintf('      Night system active: %s (%.1f%%)\n', ...
        num2sepstr(nightActiveObs, '%.0f'), 100 * nightActiveObs / totalObs);
    fprintf('      Meeting all removal criteria: %s (%.1f%%)\n', ...
        num2sepstr(removedObs, '%.0f'), 100 * removedObs / totalObs);
    
    % Additional diagnostic: show time distribution of removed observations
    if sum(rowsToRemove) > 0
        removedTimes = dateTimes(rowsToRemove);
        hourOfDay = hour(removedTimes);
        [hourCounts, hourBins] = histcounts(hourOfDay, 0:24);
        
        fprintf('      Removed observations by hour: ');
        peakHours = find(hourCounts > 0);
        for h = peakHours
            fprintf('%02d:00(%d) ', hourBins(h), hourCounts(h));
        end
        fprintf('\n');
    end
end

function [sunriseMatched, sunsetMatched] = matchSunriseSunsetTimes(dateTimes, weatherData)
    % Match each observation time with the appropriate sunrise/sunset times
    
    % Extract just the date part of each observation
    observationDates = dateshift(dateTimes, 'start', 'day');
    
    % Initialize output arrays
    sunriseMatched = NaT(size(dateTimes));
    sunsetMatched = NaT(size(dateTimes));
    
    % Match each unique date with weather data
    uniqueDates = unique(observationDates);
    
    for i = 1:length(uniqueDates)
        currentDate = uniqueDates(i);
        
        % Find matching weather data for this date
        weatherIdx = find(weatherData.dates == currentDate);
        
        if ~isempty(weatherIdx)
            % Get sunrise/sunset for this date
            if iscell(weatherData.sunrise)
                sunrise = weatherData.sunrise{weatherIdx(1)};
                sunset = weatherData.sunset{weatherIdx(1)};
            else
                sunrise = weatherData.sunrise(weatherIdx(1));
                sunset = weatherData.sunset(weatherIdx(1));
            end
            
            % Apply to all observations on this date
            dateMatches = (observationDates == currentDate);
            sunriseMatched(dateMatches) = sunrise;
            sunsetMatched(dateMatches) = sunset;
        else
            % No weather data for this date - use reasonable defaults
            % Default sunrise: 6:00 AM, sunset: 6:00 PM
            defaultSunrise = currentDate + hours(6);
            defaultSunset = currentDate + hours(18);
            
            dateMatches = (observationDates == currentDate);
            sunriseMatched(dateMatches) = defaultSunrise;
            sunsetMatched(dateMatches) = defaultSunset;
        end
    end
    
    % Handle any remaining NaT values with defaults
    stillMissing = isnat(sunriseMatched) | isnat(sunsetMatched);
    if any(stillMissing)
        fprintf('    Warning: Using default sunrise/sunset times for %d observations with missing weather data\n', ...
            sum(stillMissing));
        
        missingDates = dateshift(dateTimes(stillMissing), 'start', 'day');
        sunriseMatched(stillMissing) = missingDates + hours(6);
        sunsetMatched(stillMissing) = missingDates + hours(18);
    end
end

function demonstrateFilterEffect(originalData, filteredData, analysis)
    % Optional function to visualize the effect of night filtering
    % This can be called to show before/after comparison
    
    figure('Position', [408 126 1200 800]);
    
    % Create subplots for comparison
    subplot(2, 1, 1);
    plotHourlyDistribution(originalData, 'Original Data (with nighttime zeros)', analysis);
    
    subplot(2, 1, 2);
    plotHourlyDistribution(filteredData, 'Filtered Data (nighttime zeros removed)', analysis);
    
    sgtitle('Effect of Night Data Filtering on Hourly Count Distribution', 'FontSize', 16, 'FontWeight', 'bold');
end

function plotHourlyDistribution(data, titleStr, analysis)
    % Plot distribution of counts by hour of day
    
    if ismember(analysis.modeString, data.Properties.VariableNames)
        dateTimes = data.('Date and Time (Local)');
        counts = data.(analysis.modeString);
        hours = hour(dateTimes);
        
        % Create box plot or histogram by hour
        boxplot(counts, hours);
        xlabel('Hour of Day');
        ylabel(['Count (' analysis.modeString ')']);
        title(titleStr);
        grid on;
        
        % Add count statistics
        totalObs = length(counts);
        zeroObs = sum(counts == 0);
        text(0.02, 0.98, sprintf('Total obs: %s\nZero counts: %s (%.1f%%)', ...
            num2sepstr(totalObs, '%.0f'), num2sepstr(zeroObs, '%.0f'), ...
            100 * zeroObs / totalObs), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    else
        text(0.5, 0.5, 'Data column not found', 'HorizontalAlignment', 'center');
        title(titleStr);
    end
end

function plotHourlyCountHistogram(locationData, analysis, style)
    % Plot histogram of hourly counts for all locations on the same axes
    
    figure('Position', [408 126 1000 700]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Collect all data first to determine appropriate bin edges
    allCounts = [];
    locationCounts = cell(length(locationNames), 1);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        
        % Extract hourly counts
        hourlyData = extractHourlyData(data, analysis);
        locationCounts{i} = hourlyData.rawCounts;
        allCounts = [allCounts; hourlyData.rawCounts];
    end
    
    % Determine bin edges based on all data
    if ~isempty(allCounts)
        maxCount = max(allCounts);
        
        % Create sensible bin edges
        if maxCount <= 20
            binEdges = 0:1:maxCount+1;  % 1-count bins for small ranges
        elseif maxCount <= 50
            binEdges = 0:2:maxCount+2;  % 2-count bins for medium ranges
        elseif maxCount <= 100
            binEdges = 0:5:maxCount+5;  % 5-count bins for larger ranges
        else
            binEdges = 0:10:maxCount+10; % 10-count bins for very large ranges
        end
    else
        binEdges = 0:1:10;  % Default if no data
    end

    binEdges = 0:1:maxCount+1; % Hard coded
    
    % Plot histogram for each location
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        counts = locationCounts{i};
        
        if ~isempty(counts)
            % Create histogram
            h = histogram(counts, binEdges, ...
                'FaceColor', locationInfo.plotColor, ...
                'EdgeColor', 'white', ...
                'FaceAlpha', 0.6, ...
                'LineWidth', 1, ...
                'DisplayName', sprintf('%s (%s hours; median = %.0f per hour; max = %s per hour)', ...
                    locationInfo.name, ...
                    num2sepstr(length(counts), '%.0f'), ...
                    median(counts), ...
                    num2sepstr(max(counts), '%.0f')));
            
            plotHandles = [plotHandles, h];
        end
    end
    
    % Format the plot
    formatHistogramPlot(analysis, style, binEdges);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function plotHourlyCountHistogramNormalized(locationData, analysis, style)
    % Alternative version: Normalized histograms for better comparison when sample sizes differ
    
    figure('Position', [408 126 1000 700]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Collect all data to determine bin edges
    allCounts = [];
    locationCounts = cell(length(locationNames), 1);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        hourlyData = extractHourlyData(data, analysis);
        locationCounts{i} = hourlyData.rawCounts;
        allCounts = [allCounts; hourlyData.rawCounts];
    end
    
    % Determine bin edges
    if ~isempty(allCounts)
        maxCount = max(allCounts);
        if maxCount <= 20
            binEdges = 0:1:maxCount+1;
        elseif maxCount <= 50
            binEdges = 0:2:maxCount+2;
        elseif maxCount <= 100
            binEdges = 0:5:maxCount+5;
        else
            binEdges = 0:10:maxCount+10;
        end
    else
        binEdges = 0:1:10;
    end
    
    % Plot normalized histogram for each location
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        counts = locationCounts{i};
        
        if ~isempty(counts)
            % Create normalized histogram (shows probability density)
            h = histogram(counts, binEdges, ...
                'Normalization', 'probability', ...
                'FaceColor', locationInfo.plotColor, ...
                'EdgeColor', 'white', ...
                'FaceAlpha', 0.6, ...
                'LineWidth', 1, ...
                'DisplayName', sprintf('%s (n=%s, median=%.1f)', ...
                    locationInfo.name, ...
                    num2sepstr(length(counts), '%.0f'), ...
                    median(counts)));
            
            plotHandles = [plotHandles, h];
        end
    end
    
    % Format the normalized plot
    formatNormalizedHistogramPlot(analysis, style, binEdges);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function formatHistogramPlot(analysis, style, binEdges)
    % Format the histogram plot
    
    xlabel(['Hourly ' analysis.modeDisplayString], ...
        'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel('Number of Hours', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(['Distribution of Hourly ' analysis.modeDisplayString ' (Dawn to Dusk)'], ...
        'FontSize', style.titleFontSize);

    subtitle(sprintf('%s to %s', ...
        datestr(analysis.startTime, 'dd-mmm-yyyy'), ...
        datestr(analysis.endTime, 'dd-mmm-yyyy')), ...
        'FontSize', style.axisFontSize, 'Color', 0.3.*[1 1 1]);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    % Set x-axis limits based on bin edges
    xlim([binEdges(1) binEdges(end)]);
    xticks([binEdges(1):10:binEdges(end)+1])
    
    % Add vertical line at median for reference (position not correct)
    % medianLine = median(binEdges(1:end-1));  % Approximate median bin
    % line([medianLine medianLine], ylim, 'Color', [0.5 0.5 0.5], ...
    %     'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
end

function formatNormalizedHistogramPlot(analysis, style, binEdges)
    % Format the normalized histogram plot
    
    xlabel(['Hourly ' analysis.modeDisplayString], ...
        'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel('Probability', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(['Distribution of Hourly ' analysis.modeDisplayString ' (Normalized)'], ...
        'FontSize', style.titleFontSize);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Set x-axis limits
    xlim([binEdges(1) binEdges(end)]);
    
    % Format y-axis as percentages
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) sprintf('%.1f%%', v*100), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
end

function plotHourlyCountHistogramOverlay(locationData, analysis, style)
    % Alternative version: Overlay with transparency (good for similar distributions)
    
    figure('Position', [408 126 1000 700]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Determine common bin edges
    allCounts = [];
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        hourlyData = extractHourlyData(data, analysis);
        allCounts = [allCounts; hourlyData.rawCounts];
    end
    
    if ~isempty(allCounts)
        maxCount = max(allCounts);
        if maxCount <= 20
            binEdges = 0:1:maxCount+1;
        elseif maxCount <= 50
            binEdges = 0:2:maxCount+2;
        else
            binEdges = 0:5:maxCount+5;
        end
    else
        binEdges = 0:1:10;
    end

    binEdges = 0:1:maxCount+1; % Hard coded
    
    % Plot each location with high transparency for overlay effect
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        hourlyData = extractHourlyData(data, analysis);
        
        if ~isempty(hourlyData.rawCounts)
            h = histogram(hourlyData.rawCounts, binEdges, ...
                'FaceColor', locationInfo.plotColor, ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.5, ...  % High transparency for overlay
                'DisplayName', sprintf('%s (n=%s)', ...
                    locationInfo.name, ...
                    num2sepstr(length(hourlyData.rawCounts), '%.0f')));
            
            plotHandles = [plotHandles, h];
        end
    end
    
    % Format plot
    formatHistogramPlot(analysis, style, binEdges);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function analyzeZeroCountIntervals(locationData, analysis)
    % Find and analyze the longest intervals of consecutive zero daylight counts
    
    fprintf('\n=== Zero Count Interval Analysis for %s ===\n', analysis.modeDisplayString);
    
    locationNames = fieldnames(locationData);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        fprintf('\nLocation: %s\n', locationInfo.name);
        
        % Find zero intervals for this location
        [longestInterval, allIntervals, stats] = findZeroIntervals(data, analysis);
        
        % Report results
        reportZeroIntervals(longestInterval, allIntervals, stats, locationInfo.name, analysis);
    end
    
    fprintf('\n');
end

function [longestInterval, allIntervals, stats] = findZeroIntervals(locationDataStruct, analysis)
    % Find all intervals of consecutive zero counts during daylight hours
    
    % Extract the data
    data = locationDataStruct.data;
    
    % Filter for daylight hours only (where detection is possible)
    if ismember('Daylight', data.Properties.VariableNames)
        daylightData = data(data.Daylight == 1, :);
    else
        % Fallback: assume all data is daylight if Daylight column doesn't exist
        daylightData = data;
        fprintf('  Warning: No Daylight column found, using all data\n');
    end
    
    if isempty(daylightData)
        longestInterval = struct();
        allIntervals = [];
        stats = struct();
        return;
    end
    
    % Get times and counts
    dateTimes = daylightData.('Date and Time (Local)');
    counts = daylightData.(analysis.modeString);
    
    % Remove NaN values
    validIdx = ~isnan(counts);
    dateTimes = dateTimes(validIdx);
    counts = counts(validIdx);
    
    if isempty(counts)
        longestInterval = struct();
        allIntervals = [];
        stats = struct();
        return;
    end
    
    % Find consecutive zero periods
    isZero = (counts == 0);
    
    % Find start and end indices of zero runs
    zeroStarts = find(diff([0; isZero]) == 1);
    zeroEnds = find(diff([isZero; 0]) == -1);
    
    % Calculate intervals
    allIntervals = [];
    
    for i = 1:length(zeroStarts)
        startIdx = zeroStarts(i);
        endIdx = zeroEnds(i);
        
        startTime = dateTimes(startIdx);
        endTime = dateTimes(endIdx);
        duration = hours(endTime - startTime) + 1; % +1 because both endpoints are included
        
        interval = struct();
        interval.startTime = startTime;
        interval.endTime = endTime;
        interval.duration = duration;
        interval.startIdx = startIdx;
        interval.endIdx = endIdx;
        
        allIntervals = [allIntervals; interval];
    end

    % NEW: Filter out intervals that span multiple dates
    % These likely include nighttime periods when detection is impossible
    if ~isempty(allIntervals)
        startDates = dateshift([allIntervals.startTime], 'start', 'day');
        endDates = dateshift([allIntervals.endTime], 'start', 'day');
        sameDayIntervals = (startDates == endDates);
        
        % Keep only same-day intervals
        filteredIntervals = allIntervals(sameDayIntervals);
        
        % Report filtering statistics
        originalCount = length(allIntervals);
        filteredCount = length(filteredIntervals);
        removedCount = originalCount - filteredCount;
        
        fprintf('  Filtered out %d overnight intervals (keeping %d same-day intervals)\n', ...
            removedCount, filteredCount);
        
        % Update allIntervals to contain only same-day intervals
        allIntervals = filteredIntervals;
    end
    
    % Find longest interval
    if ~isempty(allIntervals)
        [~, maxIdx] = max([allIntervals.duration]);
        longestInterval = allIntervals(maxIdx);
    else
        longestInterval = struct();
    end
    
    % Calculate statistics
    stats = struct();
    stats.totalObservations = length(counts);
    stats.zeroObservations = sum(isZero);
    stats.nonZeroObservations = sum(~isZero);
    stats.percentZero = 100 * stats.zeroObservations / stats.totalObservations;
    stats.numZeroIntervals = length(allIntervals);
    
    if ~isempty(allIntervals)
        durations = [allIntervals.duration];
        stats.meanIntervalDuration = mean(durations);
        stats.medianIntervalDuration = median(durations);
        stats.maxIntervalDuration = max(durations);
        stats.totalZeroHours = sum(durations);
    else
        stats.meanIntervalDuration = 0;
        stats.medianIntervalDuration = 0;
        stats.maxIntervalDuration = 0;
        stats.totalZeroHours = 0;
    end
end

function reportZeroIntervals(longestInterval, allIntervals, stats, locationName, analysis)
    % Report the zero interval analysis results
    
    fprintf('  Total daylight observations: %s\n', num2sepstr(stats.totalObservations, '%.0f'));
    fprintf('  Zero count observations: %s (%.1f%%)\n', ...
        num2sepstr(stats.zeroObservations, '%.0f'), stats.percentZero);
    fprintf('  Number of zero intervals: %d\n', stats.numZeroIntervals);
    
    if stats.numZeroIntervals > 0
        fprintf('  Total hours with zero counts: %.1f\n', stats.totalZeroHours);
        fprintf('  Average zero interval duration: %.1f hours\n', stats.meanIntervalDuration);
        fprintf('  Median zero interval duration: %.1f hours\n', stats.medianIntervalDuration);
        
        fprintf('\n  LONGEST ZERO INTERVAL:\n');
        fprintf('    Start: %s\n', datestr(longestInterval.startTime, 'dd-mmm-yyyy HH:MM'));
        fprintf('    End:   %s\n', datestr(longestInterval.endTime, 'dd-mmm-yyyy HH:MM'));
        fprintf('    Duration: %.1f hours (%.1f days)\n', ...
            longestInterval.duration, longestInterval.duration / 24);
        
        % Describe the interval in context
        if longestInterval.duration >= 24 * 7
            fprintf('    Context: More than a week of no activity\n');
        elseif longestInterval.duration >= 24 * 3
            fprintf('    Context: Multiple days of no activity\n');
        elseif longestInterval.duration >= 24
            fprintf('    Context: More than a full day of no activity\n');
        elseif longestInterval.duration >= 12
            fprintf('    Context: More than half a day of no activity\n');
        else
            fprintf('    Context: Several hours of no activity\n');
        end
        
        % Show top 5 longest intervals if there are multiple
        numTop = 20;
        if length(allIntervals) > 1
            fprintf(['\n  TOP ' num2str(numTop) ' LONGEST ZERO INTERVALS:\n']);
            durations = [allIntervals.duration];
            [sortedDurations, sortIdx] = sort(durations, 'descend');
            
            numToShow = min(numTop, length(allIntervals));
            for i = 1:numToShow
                idx = sortIdx(i);
                interval = allIntervals(idx);
                fprintf('    %d. %s to %s (%.1f hours)\n', i, ...
                    datestr(interval.startTime, 'dd-mmm HH:MM'), ...
                    datestr(interval.endTime, 'dd-mmm HH:MM'), ...
                    interval.duration);
            end
        end
    else
        fprintf('  No zero intervals found - there was at least one count in every hour!\n');
    end
end

function plotZeroIntervalAnalysis(locationData, analysis, style)
    % Optional: Create visualization of zero intervals
    
    figure('Position', [408 126 1200 800]);
    
    locationNames = fieldnames(locationData);
    numLocations = length(locationNames);
    
    for i = 1:numLocations
        subplot(numLocations, 1, i);
        
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Get daylight data
        if ismember('Daylight', data.data.Properties.VariableNames)
            daylightData = data.data(data.data.Daylight == 1, :);
        else
            daylightData = data.data;
        end
        
        if ~isempty(daylightData)
            dateTimes = daylightData.('Date and Time (Local)');
            counts = daylightData.(analysis.modeString);
            
            % Remove NaN values
            validIdx = ~isnan(counts);
            dateTimes = dateTimes(validIdx);
            counts = counts(validIdx);
            
            % Plot as line with zero highlighted
            plot(dateTimes, counts, '.', 'MarkerSize', 2, 'Color', locationInfo.plotColor);
            hold on;
            
            % Highlight zero periods
            zeroIdx = (counts == 0);
            if any(zeroIdx)
                plot(dateTimes(zeroIdx), counts(zeroIdx), '.', 'MarkerSize', 4, 'Color', 'red');
            end
            
            % Find and highlight longest zero interval
            [longestInterval, ~, ~] = findZeroIntervals(data, analysis);
            if ~isempty(fieldnames(longestInterval))
                % Add rectangle to highlight longest interval
                yLimits = ylim;
                rectangle('Position', [datenum(longestInterval.startTime), yLimits(1), ...
                    datenum(longestInterval.endTime) - datenum(longestInterval.startTime), ...
                    yLimits(2) - yLimits(1)], ...
                    'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
            
            title(sprintf('%s - Zero Intervals Highlighted', locationInfo.name), ...
                'FontSize', style.titleFontSize * 0.9);
            ylabel(['Hourly ' analysis.modeDisplayString], 'FontSize', style.labelFontSize * 0.9);
            
            if i == numLocations
                xlabel('Date', 'FontSize', style.labelFontSize);
            end
            
            grid on;
            set(gca, 'FontSize', style.axisFontSize * 0.9);
            
            hold off;
        end
    end
    
    sgtitle('Zero Count Interval Analysis', 'FontSize', style.titleFontSize);
end

function plotTemperatureScatterWeekly(locationData, weatherData, analysis, plots, style)
    % Plot scatter plot of weekly counts versus average weekly air temperature for all locations
    % This is the weekly equivalent of plotTemperatureScatter()
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot each location with different colors
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate weekly totals and match with weather data
        [weeklyCounts, weeklyTemperatures, validWeekStarts] = prepareWeeklyTemperatureData(data, weatherData, analysis);
        
        if ~isempty(weeklyCounts)
            % Create scatter plot
            h = scatter(weeklyTemperatures, weeklyCounts, 80, ...
                'MarkerEdgeColor', locationInfo.plotColor, ...
                'MarkerFaceColor', locationInfo.plotColor, ...
                'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeAlpha', 0.8, ...
                'DisplayName', sprintf('%s (%d weeks)', locationInfo.name, length(weeklyCounts)));
            plotHandles = [plotHandles, h];
            
            % Add piecewise linear trend line with breakpoint at 0°C
            if length(weeklyTemperatures) > 5
                [trendParams, trendStats] = fitPiecewiseLinear(weeklyTemperatures, weeklyCounts);
                
                if ~isempty(trendParams)
                    % Plot the piecewise linear trend
                    tempRange = linspace(min(weeklyTemperatures), max(weeklyTemperatures), 100);
                    trendLine = evaluatePiecewiseLinear(tempRange, trendParams);
                    
                    plot(tempRange, trendLine, '-', ...
                        'Color', locationInfo.plotColor, ...
                        'LineWidth', 2, ...
                        'HandleVisibility', 'off');  % Don't show in legend
                    
                    % Store trend stats for later display
                    locationInfo.trendStats = trendStats;
                end
            end
        end
    end
    
    % Format plot
    formatWeeklyTemperatureScatterPlot(analysis, style);
    
    % Add correlation statistics
    addWeeklyCorrelationStats(locationData, weatherData, analysis, style);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function [weeklyCounts, weeklyTemperatures, validWeekStarts] = prepareWeeklyTemperatureData(locationDataStruct, weatherData, analysis)
    % Prepare matched weekly counts and temperature data
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Calculate weekly totals using the existing yearWeekKey approach
    groupedData = groupsummary(data, 'yearWeekKey', 'sum', analysis.modeString);
    
    % Get week start dates for each yearWeekKey
    weekStartGrouped = groupsummary(data, 'yearWeekKey', 'min', 'weekStartDateTimes');
    
    % Create weekly data table
    sumColumnName = ['sum_' analysis.modeString];
    weeklyData = table(weekStartGrouped.min_weekStartDateTimes, ...
        groupedData.(sumColumnName), ...
        'VariableNames', {'WeekStart', 'Count'});
    
    % Match indices between weekly data and week starts
    [~, ia, ib] = intersect(groupedData.yearWeekKey, weekStartGrouped.yearWeekKey);
    weeklyData = weeklyData(ib, :);  % Align the data properly
    
    % Aggregate weather data to weekly averages
    weeklyWeatherData = aggregateWeatherToWeekly(weatherData);
    
    % Match weekly counts with weekly weather data
    [~, ic, id] = intersect(weeklyData.WeekStart, weeklyWeatherData.weekStarts);
    
    if isempty(ic)
        % No matching weeks
        weeklyCounts = [];
        weeklyTemperatures = [];
        validWeekStarts = [];
        return;
    end
    
    % Extract matched data
    weeklyCounts = weeklyData.Count(ic);
    weeklyTemperatures = weeklyWeatherData.avgTemperature(id);
    validWeekStarts = weeklyData.WeekStart(ic);
    
    % Remove any NaN values
    validIdx = ~isnan(weeklyCounts) & ~isnan(weeklyTemperatures);
    weeklyCounts = weeklyCounts(validIdx);
    weeklyTemperatures = weeklyTemperatures(validIdx);
    validWeekStarts = validWeekStarts(validIdx);
end

function formatWeeklyTemperatureScatterPlot(analysis, style)
    % Format the weekly temperature scatter plot
    
    xlabel('Average Weekly Air Temperature (°C)', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel(['Weekly ' analysis.modeDisplayString], 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(['Weekly ' analysis.modeDisplayString ' vs Average Weekly Air Temperature'], ...
        'FontSize', style.titleFontSize);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    % Ensure y-axis starts at 0
    ylim([0 max(ylim) * 1.05]);
    
    % Add temperature reference lines (optional)
    xLimits = xlim;
    
    % Add vertical lines for key temperatures
    line([0 0], ylim, 'Color', [0.7 0.7 1], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    line([10 10], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    line([20 20], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Add temperature labels
    text(0, max(ylim) * 0.95, '0°C', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
    text(10, max(ylim) * 0.95, '10°C', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
    text(20, max(ylim) * 0.95, '20°C', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
end

function addWeeklyCorrelationStats(locationData, weatherData, analysis, style)
    % Add piecewise linear model statistics as text annotation for weekly data
    
    locationNames = fieldnames(locationData);
    statsText = {};
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Get weekly data for analysis
        [weeklyCounts, weeklyTemperatures, ~] = prepareWeeklyTemperatureData(data, weatherData, analysis);
        
        if length(weeklyCounts) > 5
            % Fit piecewise linear model
            [trendParams, trendStats] = fitPiecewiseLinear(weeklyTemperatures, weeklyCounts);
            
            if ~isempty(trendParams)
                % Format statistics for display
                belowFreezing = sprintf('T<0°C: slope=%.1f', trendStats.slopeBelow);
                aboveFreezing = sprintf('T≥0°C: slope=%.1f', trendStats.slopeAbove);
                rSquared = sprintf('R²=%.3f', trendStats.rSquared);
                
                statsText{end+1} = sprintf('%s: %s, %s, %s', ...
                    locationInfo.name, belowFreezing, aboveFreezing, rSquared);
            else
                % Fallback to simple correlation if piecewise fit fails
                corrCoeff = corr(weeklyTemperatures, weeklyCounts, 'type', 'Pearson');
                statsText{end+1} = sprintf('%s: r = %.3f (simple)', locationInfo.name, corrCoeff);
            end
        end
    end
    
    % Display statistics
    if ~isempty(statsText)
        % Position text box in upper left
        xLimits = xlim;
        yLimits = ylim;
        
        textX = xLimits(1) + 0.05 * (xLimits(2) - xLimits(1));
        textY = yLimits(2) - 0.1 * (yLimits(2) - yLimits(1));
        
        % Create text box
        textStr = strjoin(statsText, '\n');
        text(textX, textY, textStr, ...
            'FontSize', style.axisFontSize * 0.9, ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.9], ...
            'EdgeColor', [0.7 0.7 0.7], ...
            'Margin', 5);
    end
end

function plotLocationCorrelation(locationData, analysis, plots, style, aggregationLevel)
    % Plot correlation scatter plot between counts at two locations
    %
    % Inputs:
    %   locationData - structure containing data for all locations
    %   analysis - analysis configuration structure
    %   plots - plotting configuration structure  
    %   style - plotting style configuration structure
    %   aggregationLevel - string: 'hourly' | 'daily' | 'weekly' | 'monthly'
    %
    % This function creates a scatter plot showing the correlation between
    % traffic counts at the two monitoring locations, with data aggregated
    % according to the specified level.
    
    % Validate inputs
    if nargin < 5
        error('aggregationLevel parameter is required');
    end
    
    validLevels = {'hourly', 'daily', 'weekly', 'monthly'};
    if ~ismember(lower(aggregationLevel), validLevels)
        error('aggregationLevel must be one of: %s', strjoin(validLevels, ', '));
    end
    
    locationNames = fieldnames(locationData);
    
    % Ensure we have exactly two locations
    if length(locationNames) ~= 2
        error('This function requires exactly two locations. Found %d locations.', length(locationNames));
    end
    
    % Extract data for both locations
    location1Name = locationNames{1};
    location2Name = locationNames{2};
    
    location1Data = locationData.(location1Name);
    location2Data = locationData.(location2Name);
    
    location1Info = location1Data.locationInfo;
    location2Info = location2Data.locationInfo;
    
    % Get aggregated data for both locations based on specified level
    switch lower(aggregationLevel)
        case 'hourly'
            [counts1, counts2, timePoints, validIdx] = prepareHourlyCorrelationData(location1Data, location2Data, analysis);
            timeAxisLabel = 'Time';
            titleSuffix = 'Hourly';
            
        case 'daily'
            [counts1, counts2, timePoints, validIdx] = prepareDailyCorrelationData(location1Data, location2Data, analysis);
            timeAxisLabel = 'Date';
            titleSuffix = 'Daily';
            
        case 'weekly'
            [counts1, counts2, timePoints, validIdx] = prepareWeeklyCorrelationData(location1Data, location2Data, analysis);
            timeAxisLabel = 'Week Starting';
            titleSuffix = 'Weekly';
            
        case 'monthly'
            [counts1, counts2, timePoints, validIdx] = prepareMonthlyCorrelationData(location1Data, location2Data, analysis);
            timeAxisLabel = 'Month';
            titleSuffix = 'Monthly';
            
        otherwise
            error('Invalid aggregation level: %s', aggregationLevel);
    end
    
    % Check if we have sufficient data
    if isempty(counts1) || isempty(counts2) || length(counts1) < 3
        warning('Insufficient data for correlation analysis at %s level. Need at least 3 data points.', aggregationLevel);
        fprintf('Found %d matching data points between locations.\n', length(counts1));
        return;
    end
    
    % Create the figure
    figure('Position', [408 126 1132 921]);
    
    % Calculate correlation statistics
    [corrStats, trendStats] = calculateCorrelationStats(counts1, counts2);
    
    % Create the scatter plot
    h = scatter(counts1, counts2, 60, ...
        'MarkerFaceColor', [0.2 0.4 0.8], ...
        'MarkerEdgeColor', [0.1 0.2 0.6], ...
        'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeAlpha', 0.8);
    
    hold on;
    
    % Add trend line
    if ~isempty(trendStats)
        xRange = linspace(min(counts1), max(counts1), 100);
        trendLine = trendStats.slope * xRange + trendStats.intercept;
        
        plot(xRange, trendLine, '-', ...
            'Color', [0.8 0.2 0.2], ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('Trend: y = %.3fx + %.1f', trendStats.slope, trendStats.intercept));
    end
    
    % Add 1:1 reference line for comparison
    minVal = min([min(counts1), min(counts2)]);
    maxVal = max([max(counts1), max(counts2)]);
    plot([minVal maxVal], [minVal maxVal], '--', ...
        'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1, ...
        'DisplayName', '1:1 Reference');
    
    % Format the plot
    formatCorrelationPlot(location1Info, location2Info, analysis, style, titleSuffix, corrStats);
    
    % Add correlation statistics text box
    addCorrelationStatsBox(corrStats, trendStats, style, length(counts1));
    
    % Add legend if trend line was plotted
    if ~isempty(trendStats)
        legend('Data Points', 'Trend Line', '1:1 Reference', ...
            'Location', 'northwest', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off;
    
    % Optional: Print summary statistics to console
    printCorrelationSummary(location1Info, location2Info, corrStats, trendStats, titleSuffix, length(counts1));
end

function [counts1, counts2, timePoints, validIdx] = prepareHourlyCorrelationData(location1Data, location2Data, analysis)
    % Prepare hourly data for correlation analysis
    
    % Extract hourly data for both locations
    data1 = location1Data.data;
    data2 = location2Data.data;
    
    % Get time series and counts
    times1 = data1.('Date and Time (Local)');
    times2 = data2.('Date and Time (Local)');
    counts1_raw = data1.(analysis.modeString);
    counts2_raw = data2.(analysis.modeString);
    
    % Find common time points
    [commonTimes, ia, ib] = intersect(times1, times2);
    
    if isempty(commonTimes)
        counts1 = [];
        counts2 = [];
        timePoints = [];
        validIdx = [];
        return;
    end
    
    % Extract counts for common time points
    counts1 = counts1_raw(ia);
    counts2 = counts2_raw(ib);
    timePoints = commonTimes;
    
    % Remove NaN values
    validIdx = ~isnan(counts1) & ~isnan(counts2);
    counts1 = counts1(validIdx);
    counts2 = counts2(validIdx);
    timePoints = timePoints(validIdx);
end

function [counts1, counts2, timePoints, validIdx] = prepareDailyCorrelationData(location1Data, location2Data, analysis)
    % Prepare daily aggregated data for correlation analysis
    
    % Calculate daily totals for both locations
    dailyData1 = calculateDailyTotals(location1Data, analysis);
    dailyData2 = calculateDailyTotals(location2Data, analysis);
    
    % Find common dates
    [commonDates, ia, ib] = intersect(dailyData1.dates, dailyData2.dates);
    
    if isempty(commonDates)
        counts1 = [];
        counts2 = [];
        timePoints = [];
        validIdx = [];
        return;
    end
    
    % Extract counts for common dates
    counts1 = dailyData1.rawCounts(ia);
    counts2 = dailyData2.rawCounts(ib);
    timePoints = commonDates;
    
    % Remove NaN values
    validIdx = ~isnan(counts1) & ~isnan(counts2);
    counts1 = counts1(validIdx);
    counts2 = counts2(validIdx);
    timePoints = timePoints(validIdx);
end

function [counts1, counts2, timePoints, validIdx] = prepareWeeklyCorrelationData(location1Data, location2Data, analysis)
    % Prepare weekly aggregated data for correlation analysis
    
    % Calculate weekly totals for both locations
    weeklyData1 = calculateWeeklyTotals(location1Data, analysis);
    weeklyData2 = calculateWeeklyTotals(location2Data, analysis);
    
    % Find common week starts
    [commonWeeks, ia, ib] = intersect(weeklyData1.weekStarts, weeklyData2.weekStarts);
    
    if isempty(commonWeeks)
        counts1 = [];
        counts2 = [];
        timePoints = [];
        validIdx = [];
        return;
    end
    
    % Extract counts for common weeks
    counts1 = weeklyData1.rawCounts(ia);
    counts2 = weeklyData2.rawCounts(ib);
    timePoints = commonWeeks;
    
    % Remove NaN values
    validIdx = ~isnan(counts1) & ~isnan(counts2);
    counts1 = counts1(validIdx);
    counts2 = counts2(validIdx);
    timePoints = timePoints(validIdx);
end

function [counts1, counts2, timePoints, validIdx] = prepareMonthlyCorrelationData(location1Data, location2Data, analysis)
    % Prepare monthly aggregated data for correlation analysis
    
    % Calculate monthly totals for both locations
    monthlyData1 = calculateMonthlyTotals(location1Data, analysis);
    monthlyData2 = calculateMonthlyTotals(location2Data, analysis);
    
    % Check if either location has no monthly data
    if isempty(monthlyData1.monthStarts) || isempty(monthlyData2.monthStarts)
        counts1 = [];
        counts2 = [];
        timePoints = [];
        validIdx = [];
        return;
    end
    
    % Find common months
    [commonMonths, ia, ib] = intersect(monthlyData1.monthStarts, monthlyData2.monthStarts);
    
    if isempty(commonMonths)
        counts1 = [];
        counts2 = [];
        timePoints = [];
        validIdx = [];
        return;
    end
    
    % Extract counts for common months
    counts1 = monthlyData1.rawCounts(ia);
    counts2 = monthlyData2.rawCounts(ib);
    timePoints = commonMonths;
    
    % Remove NaN values
    validIdx = ~isnan(counts1) & ~isnan(counts2);
    counts1 = counts1(validIdx);
    counts2 = counts2(validIdx);
    timePoints = timePoints(validIdx);
end

function [corrStats, trendStats] = calculateCorrelationStats(counts1, counts2)
    % Calculate correlation and trend statistics
    
    corrStats = struct();
    trendStats = struct();
    
    % Pearson correlation
    [r, p] = corrcoef(counts1, counts2);
    corrStats.pearsonR = r(1,2);
    corrStats.pearsonP = p(1,2);
    
    % Spearman correlation (rank-based, more robust to outliers)
    corrStats.spearmanRho = corr(counts1, counts2, 'type', 'Spearman');
    
    % R-squared
    corrStats.rSquared = corrStats.pearsonR^2;
    
    % Linear regression for trend line
    try
        p = polyfit(counts1, counts2, 1);
        trendStats.slope = p(1);
        trendStats.intercept = p(2);
        
        % Calculate regression statistics
        yPred = polyval(p, counts1);
        SSres = sum((counts2 - yPred).^2);
        SStot = sum((counts2 - mean(counts2)).^2);
        trendStats.rSquared = 1 - SSres/SStot;
        
        % Standard error of regression
        trendStats.standardError = sqrt(SSres / (length(counts1) - 2));
        
    catch
        trendStats = [];
    end
    
    % Basic descriptive statistics
    corrStats.mean1 = mean(counts1);
    corrStats.mean2 = mean(counts2);
    corrStats.std1 = std(counts1);
    corrStats.std2 = std(counts2);
    corrStats.n = length(counts1);
end

function formatCorrelationPlot(location1Info, location2Info, analysis, style, titleSuffix, corrStats)
    % Format the correlation scatter plot
    
    % Create short names for axis labels
    loc1Short = extractLocationShortName(location1Info.name);
    loc2Short = extractLocationShortName(location2Info.name);
    
    xlabel([titleSuffix ' ' analysis.modeDisplayString ' - ' loc1Short], ...
        'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel([titleSuffix ' ' analysis.modeDisplayString ' - ' loc2Short], ...
        'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(sprintf('%s %s Correlation: %s vs %s', ...
        titleSuffix, analysis.modeDisplayString, loc1Short, loc2Short), ...
        'FontSize', style.titleFontSize);
    
    % Add subtitle with correlation coefficient
    subtitle(sprintf('Pearson r = %.3f (R² = %.3f, n = %d)', ...
        corrStats.pearsonR, corrStats.rSquared, corrStats.n), ...
        'FontSize', style.axisFontSize, 'Color', [0.3 0.3 0.3]);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format both axes with separators
    xtick_positions = xticks;
    xtick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), xtick_positions, 'UniformOutput', false);
    xticklabels(xtick_labels);
    
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    % Ensure both axes start at 0 and have equal scaling for easier comparison
    maxVal = max([xlim, ylim]);
    xlim([0 maxVal * 1.05]);
    ylim([0 maxVal * 1.05]);
    
    % Make the plot square for better visual comparison
    axis equal;
    axis([0 maxVal*1.05 0 maxVal*1.05]);
end

function addCorrelationStatsBox(corrStats, trendStats, style, nPoints)
    % Add a text box with correlation statistics
    
    statsText = {};
    statsText{end+1} = sprintf('Sample Size: %d', nPoints);
    statsText{end+1} = sprintf('Pearson r: %.3f', corrStats.pearsonR);
    statsText{end+1} = sprintf('R²: %.3f', corrStats.rSquared);
    statsText{end+1} = sprintf('Spearman ρ: %.3f', corrStats.spearmanRho);
    
    if corrStats.pearsonP < 0.001
        statsText{end+1} = 'p < 0.001';
    else
        statsText{end+1} = sprintf('p = %.3f', corrStats.pearsonP);
    end
    
    if ~isempty(trendStats)
        statsText{end+1} = sprintf('Slope: %.3f', trendStats.slope);
        if abs(trendStats.intercept) < 1
            statsText{end+1} = sprintf('Intercept: %.2f', trendStats.intercept);
        else
            statsText{end+1} = sprintf('Intercept: %.1f', trendStats.intercept);
        end
    end
    
    % Position text box in upper left
    xLimits = xlim;
    yLimits = ylim;
    
    textX = xLimits(1) + 0.05 * (xLimits(2) - xLimits(1));
    textY = yLimits(2) - 0.05 * (yLimits(2) - yLimits(1));
    
    % Create text box
    textStr = strjoin(statsText, '\n');
    text(textX, textY, textStr, ...
        'FontSize', style.axisFontSize * 0.9, ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.9], ...
        'EdgeColor', [0.7 0.7 0.7], ...
        'Margin', 5);
end

function printCorrelationSummary(location1Info, location2Info, corrStats, trendStats, titleSuffix, nPoints)
    % Print correlation summary to console
    
    fprintf('\n=== %s Correlation Analysis ===\n', titleSuffix);
    fprintf('Location 1: %s\n', location1Info.name);
    fprintf('Location 2: %s\n', location2Info.name);
    fprintf('Sample Size: %d\n', nPoints);
    fprintf('Pearson Correlation: r = %.3f (R² = %.3f)\n', corrStats.pearsonR, corrStats.rSquared);
    fprintf('Spearman Correlation: ρ = %.3f\n', corrStats.spearmanRho);
    
    if corrStats.pearsonP < 0.001
        fprintf('Statistical Significance: p < 0.001\n');
    else
        fprintf('Statistical Significance: p = %.3f\n', corrStats.pearsonP);
    end
    
    if ~isempty(trendStats)
        fprintf('Linear Trend: y = %.3fx + %.1f\n', trendStats.slope, trendStats.intercept);
        fprintf('Standard Error: %.2f\n', trendStats.standardError);
    end
    
    % Interpret correlation strength
    absR = abs(corrStats.pearsonR);
    if absR >= 0.9
        strength = 'very strong';
    elseif absR >= 0.7
        strength = 'strong';
    elseif absR >= 0.5
        strength = 'moderate';
    elseif absR >= 0.3
        strength = 'weak';
    else
        strength = 'very weak';
    end
    
    direction = ternary(corrStats.pearsonR > 0, 'positive', 'negative');
    fprintf('Interpretation: %s %s correlation\n', strength, direction);
    fprintf('=====================================\n\n');
end

%% Bike vs other modalities

function plotBikeModalityCorrelation(locationData, analysis, plots, style)
    % Plot correlation scatter plots between weekly bike counts and other modalities
    %
    % This function creates scatter plots showing the correlation between
    % weekly bike counts and weekly counts for other transportation modes
    % (cars and pedestrians) at both monitoring locations.
    %
    % Inputs:
    %   locationData - structure containing data for all locations
    %   analysis - analysis configuration structure
    %   plots - plotting configuration structure  
    %   style - plotting style configuration structure
    
    % Define the modalities to compare against bikes
    comparisonModes = {
        struct('columnName', 'Car Total', 'displayName', 'Car Counts', 'shortName', 'Cars');
        struct('columnName', 'Pedestrian Total', 'displayName', 'Pedestrian Counts', 'shortName', 'Pedestrians');
    };
    
    locationNames = fieldnames(locationData);
    
    % Create a separate plot for each comparison mode
    for modeIdx = 1:length(comparisonModes)
        currentMode = comparisonModes{modeIdx};
        
        % Create figure for this mode comparison
        figure('Position', [408 126 1132 921]);
        hold on;
        
        plotHandles = [];
        allBikeCounts = [];
        allModeCounts = [];
        
        % Process each location
        for i = 1:length(locationNames)
            locationName = locationNames{i};
            data = locationData.(locationName);
            locationInfo = data.locationInfo;
            
            % Get weekly bike counts and weekly counts for the comparison mode
            [bikeCounts, modeCounts, weekStarts, validIdx] = prepareWeeklyModalityData(data, analysis, currentMode.columnName);

            if ~isempty(bikeCounts) && length(bikeCounts) >= 3
                % Calculate trend line first to get slope for legend
                locationSlope = NaN; % Default if trend calculation fails
                locationHorizontalIntercept = NaN; % Default if trend calculation fails
                if length(bikeCounts) > 3
                    [trendParams, trendStats] = calculateLocationTrend(modeCounts, bikeCounts);

                    if ~isempty(trendParams)
                        locationSlope = trendParams.slope;
                        locationHorizontalIntercept = trendParams.horizontalIntercept;

                        % Plot trend line
                        xRange = linspace(min(modeCounts), max(modeCounts), 100);
                        trendLine = trendParams.slope * xRange + trendParams.intercept;

                        plot(xRange, trendLine, '--', ...
                            'Color', locationInfo.plotColor, ...
                            'LineWidth', 1.5, ...
                            'HandleVisibility', 'off');  % Don't show in legend
                    end
                end

                % Create display name with slope and horizontal intercept if available
                if ~isnan(locationSlope) && ~isnan(locationHorizontalIntercept)
                    if locationHorizontalIntercept >= 0
                        displayName = sprintf('%s (%d weeks, slope = %.3f, x-intercept = %s)', ...
                            extractLocationShortName(locationInfo.name), length(bikeCounts), ...
                            locationSlope, num2sepstr(locationHorizontalIntercept,'%.0f'));
                    else
                        % Handle negative intercepts (which don't make physical sense for counts)
                        displayName = sprintf('%s (%d weeks, slope = %.3f, x-intercept < 0)', ...
                            extractLocationShortName(locationInfo.name), length(bikeCounts), ...
                            locationSlope);
                    end
                elseif ~isnan(locationSlope)
                    displayName = sprintf('%s (%d weeks, slope = %.3f)', ...
                        extractLocationShortName(locationInfo.name), length(bikeCounts), locationSlope);
                else
                    displayName = sprintf('%s (%d weeks)', ...
                        extractLocationShortName(locationInfo.name), length(bikeCounts));
                end
                % Create scatter plot for this location
                h = scatter(modeCounts, bikeCounts, 80, ...
                    'MarkerFaceColor', locationInfo.plotColor, ...
                    'MarkerEdgeColor', locationInfo.plotColor * 0.7, ...
                    'MarkerFaceAlpha', 0.6, ...
                    'MarkerEdgeAlpha', 0.8, ...
                    'DisplayName', displayName);
                
                plotHandles = [plotHandles, h];
                
                % Collect data for overall trend analysis
                allBikeCounts = [allBikeCounts; bikeCounts];
                allModeCounts = [allModeCounts; modeCounts];
            end
        end
        
        % Individual location trend lines are already added above
        % No overall trend line needed since locations have different patterns
        
        % Format the plot
        formatModalityCorrelationPlot(currentMode, analysis, style, allModeCounts, allBikeCounts);
        
        % Remove correlation statistics box - slopes are now in legend
        % addModalityCorrelationStatsBox(corrStats, [], style, length(allBikeCounts), currentMode.shortName);
        
        % Add legend
        if ~isempty(plotHandles)
            legend(plotHandles, 'Location', 'northwest', 'Color', style.axisBackgroundColor, ...
                'FontSize', style.legendFontSize);
        end
        
        hold off;
        
        % Print summary to console (optional - can be removed if not needed)
        if ~isempty(allBikeCounts) && length(allBikeCounts) >= 3
            [corrStats, ~] = calculateCorrelationStats(allModeCounts, allBikeCounts);
            printModalityCorrelationSummary(currentMode, corrStats, [], length(allBikeCounts));
        end
    end
end

function [bikeCounts, modeCounts, weekStarts, validIdx] = prepareWeeklyModalityData(locationDataStruct, analysis, modeColumnName)
    % Prepare weekly bike counts and weekly counts for another modality
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Check if the mode column exists
    if ~ismember(modeColumnName, data.Properties.VariableNames)
        warning('Mode column "%s" not found in data', modeColumnName);
        bikeCounts = [];
        modeCounts = [];
        weekStarts = [];
        validIdx = [];
        return;
    end
    
    % Calculate weekly totals for bikes using existing function
    bikeAnalysis = analysis;
    bikeAnalysis.modeString = 'Bike Total';
    weeklyBikeData = calculateWeeklyTotals(locationDataStruct, bikeAnalysis);
    
    % Calculate weekly totals for the comparison mode
    modeAnalysis = analysis;
    modeAnalysis.modeString = modeColumnName;
    weeklyModeData = calculateWeeklyTotals(locationDataStruct, modeAnalysis);
    
    % Find common week starts
    [commonWeeks, ia, ib] = intersect(weeklyBikeData.weekStarts, weeklyModeData.weekStarts);
    
    if isempty(commonWeeks)
        bikeCounts = [];
        modeCounts = [];
        weekStarts = [];
        validIdx = [];
        return;
    end
    
    % Extract counts for common weeks
    bikeCounts = weeklyBikeData.rawCounts(ia);
    modeCounts = weeklyModeData.rawCounts(ib);
    weekStarts = commonWeeks;
    
    % Remove NaN values
    validIdx = ~isnan(bikeCounts) & ~isnan(modeCounts);
    bikeCounts = bikeCounts(validIdx);
    modeCounts = modeCounts(validIdx);
    weekStarts = weekStarts(validIdx);
end

function [trendParams, trendStats] = calculateLocationTrend(xData, yData)
    % Calculate trend line parameters and statistics
    
    trendParams = struct();
    trendStats = struct();
    
    if length(xData) < 3
        trendParams = [];
        trendStats = [];
        return;
    end
    
    try
        % Linear regression
        p = polyfit(xData, yData, 1);
        trendParams.slope = p(1);
        trendParams.intercept = p(2);
        
        % Calculate horizontal intercept (x-intercept: where y = 0)
        % From y = mx + b, when y = 0, x = -b/m
        if abs(trendParams.slope) > 0.001  % Avoid division by very small numbers
            trendParams.horizontalIntercept = -trendParams.intercept / trendParams.slope;
        else
            trendParams.horizontalIntercept = NaN;  % Undefined for near-zero slope
        end
        
        % Calculate regression statistics
        yPred = polyval(p, xData);
        SSres = sum((yData - yPred).^2);
        SStot = sum((yData - mean(yData)).^2);
        
        if SStot > 0
            trendStats.rSquared = 1 - SSres/SStot;
        else
            trendStats.rSquared = 0;
        end
        
        % Standard error of regression
        if length(xData) > 2
            trendStats.standardError = sqrt(SSres / (length(xData) - 2));
        else
            trendStats.standardError = 0;
        end
        
        % Calculate correlation coefficient
        trendStats.correlation = corr(xData, yData);
        
    catch ME
        warning('Error calculating trend: %s', ME.message);
        trendParams = [];
        trendStats = [];
    end
end

function formatModalityCorrelationPlot(currentMode, analysis, style, modeCounts, bikeCounts)
    % Format the modality correlation scatter plot
    
    xlabel(['Weekly ' currentMode.displayName], ...
        'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel('Weekly Bike Counts', ...
        'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(sprintf('Weekly Bike Counts vs %s', currentMode.displayName), ...
        'FontSize', style.titleFontSize);
    
    % Add subtitle with date range
    subtitle(sprintf('%s to %s', ...
        datestr(analysis.startTime, 'mmm dd, yyyy'), ...
        datestr(analysis.endTime, 'mmm dd, yyyy')), ...
        'FontSize', style.axisFontSize, 'Color', [0.3 0.3 0.3]);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format both axes with separators
    if ~isempty(modeCounts)
        xtick_positions = xticks;
        xtick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), xtick_positions, 'UniformOutput', false);
        xticklabels(xtick_labels);
    end
    
    if ~isempty(bikeCounts)
        ytick_positions = yticks;
        ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
        yticklabels(ytick_labels);
    end
    
    % Ensure both axes start at 0
    if ~isempty(modeCounts) && ~isempty(bikeCounts)
        xlim([0 max(xlim) * 1.05]);
        ylim([0 max(ylim) * 1.05]);
    end
end

function addModalityCorrelationStatsBox(corrStats, trendStats, style, nPoints, modeName)
    % Add a text box with correlation statistics
    
    statsText = {};
    statsText{end+1} = sprintf('Sample Size: %d weeks', nPoints);
    statsText{end+1} = sprintf('Pearson r: %.3f', corrStats.pearsonR);
    statsText{end+1} = sprintf('R²: %.3f', corrStats.rSquared);
    statsText{end+1} = sprintf('Spearman ρ: %.3f', corrStats.spearmanRho);
    
    if corrStats.pearsonP < 0.001
        statsText{end+1} = 'p < 0.001';
    else
        statsText{end+1} = sprintf('p = %.3f', corrStats.pearsonP);
    end
    
    if ~isempty(trendStats) && isfield(trendStats, 'slope') && isfield(trendStats, 'intercept')
        statsText{end+1} = sprintf('Slope: %.3f', trendStats.slope);
        if abs(trendStats.intercept) < 1
            statsText{end+1} = sprintf('Intercept: %.2f', trendStats.intercept);
        else
            statsText{end+1} = sprintf('Intercept: %.1f', trendStats.intercept);
        end
    end
    
    % Add interpretation
    absR = abs(corrStats.pearsonR);
    if absR >= 0.7
        interpretation = 'Strong';
    elseif absR >= 0.5
        interpretation = 'Moderate';
    elseif absR >= 0.3
        interpretation = 'Weak';
    else
        interpretation = 'Very weak';
    end
    
    direction = ternary(corrStats.pearsonR > 0, 'positive', 'negative');
    statsText{end+1} = sprintf('%s %s correlation', interpretation, direction);
    
    % Position text box in upper left
    xLimits = xlim;
    yLimits = ylim;
    
    textX = xLimits(1) + 0.05 * (xLimits(2) - xLimits(1));
    textY = yLimits(2) - 0.05 * (yLimits(2) - yLimits(1));
    
    % Create text box
    textStr = strjoin(statsText, '\n');
    text(textX, textY, textStr, ...
        'FontSize', style.axisFontSize * 0.9, ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.9], ...
        'EdgeColor', [0.7 0.7 0.7], ...
        'Margin', 5);
end

function printModalityCorrelationSummary(currentMode, corrStats, trendStats, nPoints)
    % Print correlation summary to console
    
    fprintf('\n=== Weekly Bike vs %s Correlation Analysis ===\n', currentMode.shortName);
    fprintf('Sample Size: %d weeks\n', nPoints);
    fprintf('Pearson Correlation: r = %.3f (R² = %.3f)\n', corrStats.pearsonR, corrStats.rSquared);
    fprintf('Spearman Correlation: ρ = %.3f\n', corrStats.spearmanRho);
    
    if corrStats.pearsonP < 0.001
        fprintf('Statistical Significance: p < 0.001\n');
    else
        fprintf('Statistical Significance: p = %.3f\n', corrStats.pearsonP);
    end
    
    if ~isempty(trendStats) && isfield(trendStats, 'slope') && isfield(trendStats, 'intercept')
        fprintf('Linear Trend: bikes = %.3f × %s + %.1f\n', ...
            trendStats.slope, lower(currentMode.shortName), trendStats.intercept);
        if isfield(trendStats, 'standardError')
            fprintf('Standard Error: %.2f\n', trendStats.standardError);
        end
    end
    
    % Interpret correlation strength
    absR = abs(corrStats.pearsonR);
    if absR >= 0.7
        strength = 'strong';
    elseif absR >= 0.5
        strength = 'moderate';
    elseif absR >= 0.3
        strength = 'weak';
    else
        strength = 'very weak';
    end
    
    direction = ternary(corrStats.pearsonR > 0, 'positive', 'negative');
    fprintf('Interpretation: %s %s correlation\n', strength, direction);
    
    % Provide practical interpretation
    if strcmp(currentMode.shortName, 'Cars') && corrStats.pearsonR > 0.3
        fprintf('Practical meaning: Higher car traffic tends to coincide with higher bike traffic\n');
    elseif strcmp(currentMode.shortName, 'Pedestrians') && corrStats.pearsonR > 0.3
        fprintf('Practical meaning: Higher pedestrian activity tends to coincide with higher bike activity\n');
    end
    
    fprintf('=================================================\n\n');
end

function plotVisibilityScatterWeekly(locationData, weatherData, analysis, plots, style)
    % Plot scatter plot of weekly counts versus average weekly visibility for all locations
    % This is the weekly equivalent of plotTemperatureScatterWeekly() but for visibility
    
    figure('Position', [408 126 1132 921]);
    hold on
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot each location with different colors
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate weekly totals and match with weather data
        [weeklyCounts, weeklyVisibility, validWeekStarts] = prepareWeeklyVisibilityData(data, weatherData, analysis);
        
        if ~isempty(weeklyCounts)
            % Create scatter plot
            h = scatter(weeklyVisibility, weeklyCounts, 80, ...
                'MarkerEdgeColor', locationInfo.plotColor, ...
                'MarkerFaceColor', locationInfo.plotColor, ...
                'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeAlpha', 0.8, ...
                'DisplayName', sprintf('%s (%d weeks)', locationInfo.name, length(weeklyCounts)));
            plotHandles = [plotHandles, h];
            
            % Add simple linear trend line (no breakpoint needed for visibility)
            if length(weeklyVisibility) > 5
                [trendParams, trendStats] = fitLinearTrend(weeklyVisibility, weeklyCounts);
                
                if ~isempty(trendParams)
                    % Plot the linear trend
                    visRange = linspace(min(weeklyVisibility), max(weeklyVisibility), 100);
                    trendLine = trendParams.slope * visRange + trendParams.intercept;
                    
                    plot(visRange, trendLine, '-', ...
                        'Color', locationInfo.plotColor, ...
                        'LineWidth', 2, ...
                        'HandleVisibility', 'off');  % Don't show in legend
                    
                    % Store trend stats for later display
                    locationInfo.trendStats = trendStats;
                end
            end
        end
    end
    
    % Format plot
    formatWeeklyVisibilityScatterPlot(analysis, style);
    
    % Add correlation statistics
    addWeeklyVisibilityCorrelationStats(locationData, weatherData, analysis, style);
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'northeast', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function [weeklyCounts, weeklyVisibility, validWeekStarts] = prepareWeeklyVisibilityData(locationDataStruct, weatherData, analysis)
    % Prepare matched weekly counts and visibility data
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Calculate weekly totals using the existing yearWeekKey approach
    groupedData = groupsummary(data, 'yearWeekKey', 'sum', analysis.modeString);
    
    % Get week start dates for each yearWeekKey
    weekStartGrouped = groupsummary(data, 'yearWeekKey', 'min', 'weekStartDateTimes');
    
    % Create weekly data table
    sumColumnName = ['sum_' analysis.modeString];
    weeklyData = table(weekStartGrouped.min_weekStartDateTimes, ...
        groupedData.(sumColumnName), ...
        'VariableNames', {'WeekStart', 'Count'});
    
    % Match indices between weekly data and week starts
    [~, ia, ib] = intersect(groupedData.yearWeekKey, weekStartGrouped.yearWeekKey);
    weeklyData = weeklyData(ib, :);  % Align the data properly
    
    % Aggregate weather data to weekly averages
    weeklyWeatherData = aggregateWeatherToWeeklyWithVisibility(weatherData);
    
    % Match weekly counts with weekly weather data
    [~, ic, id] = intersect(weeklyData.WeekStart, weeklyWeatherData.weekStarts);
    
    if isempty(ic)
        % No matching weeks
        weeklyCounts = [];
        weeklyVisibility = [];
        validWeekStarts = [];
        return;
    end
    
    % Extract matched data
    weeklyCounts = weeklyData.Count(ic);
    weeklyVisibility = weeklyWeatherData.avgVisibility(id);
    validWeekStarts = weeklyData.WeekStart(ic);
    
    % Remove any NaN values
    validIdx = ~isnan(weeklyCounts) & ~isnan(weeklyVisibility);
    weeklyCounts = weeklyCounts(validIdx);
    weeklyVisibility = weeklyVisibility(validIdx);
    validWeekStarts = validWeekStarts(validIdx);
end

function weeklyWeatherData = aggregateWeatherToWeeklyWithVisibility(weatherData)
    % Create week grouping for weather data including visibility
    tempTable = table(weatherData.dates, weatherData.temperature, weatherData.precipitation, weatherData.visibility, ...
        'VariableNames', {'dates', 'temperature', 'precipitation', 'visibility'});
    
    % Add week start dates
    tempTable.weekStarts = dateshift(dateshift(tempTable.dates,'dayofweek','Monday','previous'),'start','day');
    
    % Group by week
    weeklyGrouped = groupsummary(tempTable, 'weekStarts', {'mean', 'sum'}, {'temperature', 'precipitation', 'visibility'});
    
    weeklyWeatherData = struct();
    weeklyWeatherData.weekStarts = weeklyGrouped.weekStarts;
    weeklyWeatherData.avgTemperature = weeklyGrouped.mean_temperature;
    weeklyWeatherData.totalPrecipitation = weeklyGrouped.sum_precipitation;
    weeklyWeatherData.avgVisibility = weeklyGrouped.mean_visibility;
end

function [trendParams, trendStats] = fitLinearTrend(xData, yData)
    % Fit simple linear model (no breakpoint needed for visibility)
    
    try
        % Linear fit: y = a*x + b
        p = polyfit(xData, yData, 1);
        
        % Store parameters
        trendParams = struct();
        trendParams.slope = p(1);
        trendParams.intercept = p(2);
        
        % Calculate goodness of fit
        predictedCounts = polyval(p, xData);
        
        % Calculate R-squared
        SSres = sum((yData - predictedCounts).^2);
        SStot = sum((yData - mean(yData)).^2);
        rSquared = 1 - SSres/SStot;
        
        % Store statistics
        trendStats = struct();
        trendStats.slope = p(1);
        trendStats.intercept = p(2);
        trendStats.rSquared = rSquared;
        trendStats.nPoints = length(xData);
        
    catch ME
        fprintf('Error fitting linear model: %s\n', ME.message);
        trendParams = [];
        trendStats = [];
    end
end

function formatWeeklyVisibilityScatterPlot(analysis, style)
    % Format the weekly visibility scatter plot
    
    xlabel('Average Weekly Visibility (km)', 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    ylabel(['Weekly ' analysis.modeDisplayString], 'FontSize', style.labelFontSize, 'FontWeight', 'bold');
    
    title(['Weekly ' analysis.modeDisplayString ' vs Average Weekly Visibility'], ...
        'FontSize', style.titleFontSize);
    
    set(gca, 'Color', style.axisBackgroundColor);
    set(gca, 'FontSize', style.axisFontSize);
    grid on;
    
    % Format y-axis with separators
    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);
    
    % Ensure y-axis starts at 0
    ylim([0 max(ylim) * 1.05]);
    
    % Add visibility reference lines (optional)
    xLimits = xlim;
    
    % Add vertical lines for key visibility levels
    % Poor visibility: < 1 km
    % Moderate visibility: 1-5 km  
    % Good visibility: 5-10 km
    % Excellent visibility: > 10 km
    line([1 1], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    line([5 5], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    line([10 10], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Add visibility labels
    text(1, max(ylim) * 0.95, 'Poor', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
    text(5, max(ylim) * 0.95, 'Moderate', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
    text(10, max(ylim) * 0.95, 'Good', 'HorizontalAlignment', 'center', 'FontSize', style.axisFontSize * 0.8, 'Color', [0.6 0.6 0.6]);
end

function addWeeklyVisibilityCorrelationStats(locationData, weatherData, analysis, style)
    % Add linear model statistics as text annotation for weekly visibility data
    
    locationNames = fieldnames(locationData);
    statsText = {};
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Get weekly data for analysis
        [weeklyCounts, weeklyVisibility, ~] = prepareWeeklyVisibilityData(data, weatherData, analysis);
        
        if length(weeklyCounts) > 5
            % Fit linear model
            [trendParams, trendStats] = fitLinearTrend(weeklyVisibility, weeklyCounts);
            
            if ~isempty(trendParams)
                % Format statistics for display
                slopeStr = sprintf('slope=%.1f', trendStats.slope);
                rSquared = sprintf('R²=%.3f', trendStats.rSquared);
                
                statsText{end+1} = sprintf('%s: %s, %s', ...
                    locationInfo.name, slopeStr, rSquared);
            else
                % Fallback to simple correlation
                corrCoeff = corr(weeklyVisibility, weeklyCounts, 'type', 'Pearson');
                statsText{end+1} = sprintf('%s: r = %.3f', locationInfo.name, corrCoeff);
            end
        end
    end
    
    % Display statistics
    if ~isempty(statsText)
        % Position text box in upper left
        xLimits = xlim;
        yLimits = ylim;
        
        textX = xLimits(1) + 0.05 * (xLimits(2) - xLimits(1));
        textY = yLimits(2) - 0.1 * (yLimits(2) - yLimits(1));
        
        % Create text box
        textStr = strjoin(statsText, '\n');
        text(textX, textY, textStr, ...
            'FontSize', style.axisFontSize * 0.9, ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.9], ...
            'EdgeColor', [0.7 0.7 0.7], ...
            'Margin', 5);
    end
end

