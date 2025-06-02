%% Modular Telraam Analysis - Single Script Version
% This version can be run as a single script while maintaining modularity

clear all; close all; clc;

%% Configuration Setup
% Locations to analyze
locations = {
    struct('name', 'Eastern Segment', ...
           'fileStem2024', 'telraam-raw-data-2024East60', ...
           'fileStem2025', 'telraam-raw-data-2025East60', ...
           'plotColor', [0 0 1]);  % Blue
    struct('name', 'Western Segment', ...
           'fileStem2024', 'telraam-raw-data-2024West60', ...
           'fileStem2025', 'telraam-raw-data-2025West60', ...
           'plotColor', [0 0 0]);  % Black
};

% Analysis parameters
analysis = struct( ...
    'startTime', datetime(2024,08,01,00,00,01), ...
    'endTime', datetime(2025,06,01,23,59,59), ...
    'modeString', 'Bike Total', ...
    'modeDisplayString', 'Bike Counts', ...
    'uptimeThreshold', 0.0, ...
    'maxUptimeCorrection', 1.0, ...
    'truncationCutoffTime', timeofday(datetime('today')+hours(15)), ...
    'daylightCorrectionRatioWD', 1.65, ...
    'daylightCorrectionRatioWE', 1.37 ...
);

% Plot configuration - easy to turn elements on/off
plots = struct( ...
    'showRawCounts', true, ...
    'showAdjustedCounts', false, ...  % Turn off for first example
    'showTruncatedCounts', false, ...
    'showWeather', true, ...
    'showPrecipitationBubbles', true, ...  % Turn off for first example
    'plotTypes', {{'daily'}}, ...  % Only daily for now
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
    dates = locationData.(locationNames{i}).data.Date;
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

%% ======================== FUNCTIONS ========================

function inputTable = loadSingleLocationData(location, analysis)
    % Load 2024 data
    inputTable2024 = loadYearData(location.fileStem2024, 2024);
    
    % Load 2025 data
    inputTable2025 = loadYearData(location.fileStem2025, 2025);
    
    % Combine and convert to timetable
    columnsToKeep = {'Uptime','Date','Bike Total','PedestrianTotal',...
                     'Night Total','Speed V85 km/h','Car Total','Large vehicle Total'};
    inputTable2024 = inputTable2024(:,columnsToKeep);
    inputTable2025 = inputTable2025(:,columnsToKeep);
    
    combinedTable = vertcat(inputTable2024, inputTable2025);
    inputTable = table2timetable(combinedTable);
    
    % Apply date range filter
    inputTable = inputTable(inputTable.Date >= analysis.startTime & ...
                           inputTable.Date <= analysis.endTime, :);
end

function yearTable = loadYearData(fileStem, targetYear)
    excelFileName = [fileStem '.xlsx'];
    
    % Read table with preserved variable names to avoid the warning
    opts = detectImportOptions(excelFileName);
    opts.VariableNamingRule = 'preserve';
    yearTable = readtable(excelFileName, opts);
    
    % Debug: Show available column names
    %fprintf('Available columns in %s:\n', excelFileName);
    %disp(yearTable.Properties.VariableNames);
    
    % Find the date column - try multiple possible names
    dateColumns = {'DateAndTime (Local)', 'DateAndTime_Local_', 'Date', 'DateAndTime', 'DateTime','Date and Time (Local)'};
    dateColFound = false;
    
    for i = 1:length(dateColumns)
        if ismember(dateColumns{i}, yearTable.Properties.VariableNames)
            if ~strcmp(dateColumns{i}, 'Date')  % Only rename if it's not already 'Date'
                yearTable = renamevars(yearTable, dateColumns{i}, 'Date');
            end
            dateColFound = true;
            fprintf('Found date column: %s\n', dateColumns{i});
            break;
        end
    end
    
    if ~dateColFound
        error('Could not find a date column in %s. Available columns are: %s', ...
            excelFileName, strjoin(yearTable.Properties.VariableNames, ', '));
    end
    
    % Find and rename speed column
    speedColumns = {'SpeedV85 (Km/h)', 'SpeedV85Km_h', 'SpeedV85', 'Speed V85', 'V85','Speed V85 km/h'};
    speedColFound = false;
    
    for i = 1:length(speedColumns)
        if ismember(speedColumns{i}, yearTable.Properties.VariableNames)
            if ~strcmp(speedColumns{i}, 'SpeedV85')  % Only rename if it's not already 'SpeedV85'
                yearTable = renamevars(yearTable, speedColumns{i}, 'Speed V85 km/h');
            end
            speedColFound = true;
            fprintf('Found speed column: %s\n', speedColumns{i});
            break;
        end
    end
    
    if ~speedColFound
        fprintf('Warning: No speed column found. Creating dummy SpeedV85 column.\n');
        yearTable.SpeedV85 = nan(height(yearTable), 1);
    end
    
    % Compute missing variables
    if ismember('Uptime', yearTable.Properties.VariableNames)
        yearTable.Uptime = ones(size(yearTable.Uptime));
    else
        yearTable.Uptime = ones(height(yearTable), 1);
    end
    
    % Check if the required columns exist before computing - handle original names
    nightColsA = {'ModeNight (A→B)', 'ModeNight_A_B_'};
    nightColsB = {'ModeNight (B→A)', 'ModeNight_B_A_'};
    
    nightColA = '';
    nightColB = '';
    for i = 1:length(nightColsA)
        if ismember(nightColsA{i}, yearTable.Properties.VariableNames)
            nightColA = nightColsA{i};
            break;
        end
    end
    for i = 1:length(nightColsB)
        if ismember(nightColsB{i}, yearTable.Properties.VariableNames)
            nightColB = nightColsB{i};
            break;
        end
    end
    
    if ~isempty(nightColA) && ~isempty(nightColB)
        yearTable.NightTotal_10Modes_ = yearTable.(nightColA) + yearTable.(nightColB);
    else
        yearTable.NightTotal_10Modes_ = zeros(height(yearTable), 1);
    end
    
    % Handle bicycle columns
    bikeColsA = {'ModeBicycle (A→B)', 'ModeBicycle_A_B_'};
    bikeColsB = {'ModeBicycle (B→A)', 'ModeBicycle_B_A_'};
    
    bikeColA = '';
    bikeColB = '';
    for i = 1:length(bikeColsA)
        if ismember(bikeColsA{i}, yearTable.Properties.VariableNames)
            bikeColA = bikeColsA{i};
            break;
        end
    end
    for i = 1:length(bikeColsB)
        if ismember(bikeColsB{i}, yearTable.Properties.VariableNames)
            bikeColB = bikeColsB{i};
            break;
        end
    end
    
    if ~isempty(bikeColA) && ~isempty(bikeColB)
        yearTable.BicycleTotal_10Modes_ = yearTable.(bikeColA) + yearTable.(bikeColB);
    else
        yearTable.BicycleTotal_10Modes_ = zeros(height(yearTable), 1);
    end
    
    % Handle car columns
    carColsA = {'ModeCar (A→B)', 'ModeCar_A_B_', 'CarTotal'};
    carColsB = {'ModeCar (B→A)', 'ModeCar_B_A_'};
    
    carColA = '';
    carColB = '';
    for i = 1:length(carColsA)
        if ismember(carColsA{i}, yearTable.Properties.VariableNames)
            carColA = carColsA{i};
            break;
        end
    end
    for i = 1:length(carColsB)
        if ismember(carColsB{i}, yearTable.Properties.VariableNames)
            carColB = carColsB{i};
            break;
        end
    end
    
    if ~isempty(carColA) && ~isempty(carColB)
        yearTable.CarTotal = yearTable.(carColA) + yearTable.(carColB);
    elseif ismember('CarTotal', yearTable.Properties.VariableNames)
        % CarTotal already exists, keep it
    else
        yearTable.CarTotal = zeros(height(yearTable), 1);
    end
    
    % Handle pedestrian columns
    pedCols = {'PedestrianTotal', 'Pedestrian Total'};
    pedCol = '';
    for i = 1:length(pedCols)
        if ismember(pedCols{i}, yearTable.Properties.VariableNames)
            pedCol = pedCols{i};
            break;
        end
    end
    
    if ~isempty(pedCol)
        yearTable.PedestrianTotal = yearTable.(pedCol);
    else
        yearTable.PedestrianTotal = zeros(height(yearTable), 1);
    end
    
    % Convert dates and filter by year
    yearTable.Date = datetime(char(yearTable.Date),'Format','yy-MM-dd HH:mm');
    yearTable = yearTable(year(yearTable.Date) == targetYear, :);
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
    
    inputTable.dayOfWeek = string(day(inputTable.Date,'name'));
    inputTable.dayOfWeekCat = categorical(inputTable.dayOfWeek);
    inputTable.weekOfYear = week(inputTable.Date,'iso-weekofyear');
    inputTable.isWeekday = ismember(inputTable.dayOfWeekCat, weekdays);
    inputTable.weekStartDateTimes = dateshift(dateshift(inputTable.Date,'dayofweek','Monday','previous'),'start','day');
    inputTable.Daylight = ~((inputTable.('Night Total') > 0) | (isnan(inputTable.('Speed V85 km/h'))));
    inputTable.DaylightUptime = inputTable.Daylight .* inputTable.Uptime;
    
    % Fix week numbering
    inputTable = fixWeekNumbering(inputTable);
end

function inputTable = fixWeekNumbering(inputTable)
    inputTable.weekOfYear((month(inputTable.Date)==12) & (inputTable.weekOfYear==1)) = 53;
    inputTable.weekOfYear((month(inputTable.Date)==1) & (inputTable.weekOfYear==1)) = 53;
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
    inputTable = sortrows([inputTableWD; inputTableWE], 'Date');
end

function plotCombinedDaily(locationData, weatherData, analysis, plots, style)
    figure('Position', [408 126 1132 921]);
    
    locationNames = fieldnames(locationData);
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    if plots.showWeather
        yyaxis right
        weatherHandles = plotWeatherData(weatherData, plots, style);
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
    
    % Add legend
    if ~isempty(plotHandles)
        legend(plotHandles, 'Location', 'southwest', 'Color', style.axisBackgroundColor, ...
            'FontSize', style.legendFontSize);
    end
    
    hold off
end

function dailyData = calculateDailyTotals(locationDataStruct, analysis)
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Create daily grouping
    data.DayOnly = dateshift(data.Date, 'start', 'day');
    
    % Calculate raw daily totals
    groupedData = groupsummary(data, 'DayOnly', 'sum', analysis.modeString);
    
    % Calculate adjusted daily totals (for times up to cutoff)
    truncatedData = data(timeofday(data.Date) <= analysis.truncationCutoffTime, :);
    adjustedGrouped = groupsummary(truncatedData, 'DayOnly', 'sum', 'AdjustedCountsUptimeDaylight');
    
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
    title(['Daily ' analysis.modeDisplayString ' Comparison'], ...
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