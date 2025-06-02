%% Modular Telraam Analysis Framework
% This refactored approach allows for flexible multi-location analysis
% and configurable plotting options

%% Main Analysis Script
function runTelraamAnalysis()
    % Configure analysis
    config = setupAnalysisConfig();
    
    % Load and process data for all locations
    locationData = loadAllLocationData(config);
    
    % Get weather data (once for all locations)
    weatherData = getWeatherForAnalysis(locationData, config);
    
    % Generate plots based on configuration
    generateConfiguredPlots(locationData, weatherData, config);
end

%% Configuration Setup
function config = setupAnalysisConfig()
    config = struct();
    
    % Locations to analyze
    config.locations = {
        struct('name', 'Eastern Segment', ...
               'fileStem2024', 'telraam-raw-data-13690-2024East60-f809348', ...
               'fileStem2025', 'telraam-raw-data-13690-East60-feb4c7', ...
               'plotColor', [0 0 1]);  % Blue
        struct('name', 'Western Segment', ...
               'fileStem2024', 'telraam-raw-data-9000007290-2024West60-a0fb917', ...
               'fileStem2025', 'telraam-raw-data-13374-West60-401b954', ...
               'plotColor', [0 0 0]);  % Black
    };
    
    % Analysis parameters
    config.analysis = struct( ...
        'startTime', datetime(2024,08,01,00,00,01), ...
        'endTime', datetime(2025,06,01,23,59,59), ...
        'modeString', 'BicycleTotal_10Modes_', ...
        'modeDisplayString', 'Bike Counts', ...
        'uptimeThreshold', 0.0, ...
        'maxUptimeCorrection', 1.0, ...
        'truncationCutoffTime', timeofday(datetime('today')+hours(15)), ...
        'computeDaylightCorrection', false, ...
        'daylightCorrectionRatioWD', 1.65, ...
        'daylightCorrectionRatioWE', 1.37 ...
    );
    
    % Plot configuration - easy to turn elements on/off
    config.plots = struct( ...
        'showRawCounts', true, ...
        'showAdjustedCounts', true, ...
        'showTruncatedCounts', false, ...  % "up to 3pm" plots
        'showWeather', true, ...
        'showPrecipitationBubbles', true, ...
        'plotTypes', {{'daily', 'weekly', 'monthly'}}, ...  % Which plot types to generate
        'combinedPlots', true ...  % Plot multiple locations on same axes
    );
    
    % Plotting style parameters
    config.style = struct( ...
        'plotLineWidth', 10.0, ...
        'axisFontSize', 16.0, ...
        'labelFontSize', 20.0, ...
        'titleFontSize', 24.0, ...
        'legendFontSize', 16.0, ...
        'axisBackgroundColor', 0.8.*[1 1 1], ...
        'legendBackgroundAlpha', 0.2 ...
    );
end

%% Data Loading Functions
function locationData = loadAllLocationData(config)
    locationData = struct();
    
    for i = 1:length(config.locations)
        location = config.locations{i};
        fprintf('Loading data for %s...\n', location.name);
        
        % Load raw data
        rawData = loadSingleLocationData(location, config);
        
        % Process data
        processedData = processTelraamData(rawData, config);
        
        % Store in structure
        locationData.(matlab.lang.makeValidName(location.name)) = processedData;
        locationData.(matlab.lang.makeValidName(location.name)).locationInfo = location;
    end
end

function inputTable = loadSingleLocationData(location, config)
    % Load 2024 data
    inputTable2024 = loadYearData(location.fileStem2024, 2024);
    
    % Load 2025 data
    inputTable2025 = loadYearData(location.fileStem2025, 2025);
    
    % Combine and convert to timetable
    columnsToKeep = {'Uptime','Date','BicycleTotal_10Modes_','PedestrianTotal',...
                     'NightTotal_10Modes_','SpeedV85','CarTotal'};
    inputTable2024 = inputTable2024(:,columnsToKeep);
    inputTable2025 = inputTable2025(:,columnsToKeep);
    
    combinedTable = vertcat(inputTable2024, inputTable2025);
    inputTable = table2timetable(combinedTable);
    
    % Apply date range filter
    inputTable = inputTable(inputTable.Date >= config.analysis.startTime & ...
                           inputTable.Date <= config.analysis.endTime, :);
end

function yearTable = loadYearData(fileStem, year)
    excelFileName = [fileStem '.xlsx'];
    yearTable = readtable(excelFileName);
    
    % Rename variables
    yearTable = renamevars(yearTable,'DateAndTime_Local_','Date');
    yearTable = renamevars(yearTable,'SpeedV85Km_h','SpeedV85');
    
    % Compute missing variables
    yearTable.Uptime = ones(size(yearTable.Uptime));
    yearTable.NightTotal_10Modes_ = yearTable.ModeNight_A_B_ + yearTable.ModeNight_B_A_;
    yearTable.BicycleTotal_10Modes_ = yearTable.ModeBicycle_A_B_ + yearTable.ModeBicycle_B_A_;
    
    % Convert dates and filter by year
    yearTable.Date = datetime(char(yearTable.Date),'Format','yy-MM-dd HH:mm');
    yearTable = yearTable(year(yearTable.Date) == year, :);
end

%% Data Processing Functions
function processedTable = processTelraamData(inputTable, config)
    % Add descriptive columns
    processedTable = addTemporalColumns(inputTable);
    
    % Apply uptime corrections
    processedTable = applyUptimeCorrections(processedTable, config);
    
    % Apply daylight corrections
    processedTable = applyDaylightCorrections(processedTable, config);
end

function inputTable = addTemporalColumns(inputTable)
    weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
    
    inputTable.dayOfWeek = string(day(inputTable.Date,'name'));
    inputTable.dayOfWeekCat = categorical(inputTable.dayOfWeek);
    inputTable.weekOfYear = week(inputTable.Date,'iso-weekofyear');
    inputTable.isWeekday = ismember(inputTable.dayOfWeekCat, weekdays);
    inputTable.weekStartDateTimes = dateshift(dateshift(inputTable.Date,'dayofweek','Monday','previous'),'start','day');
    inputTable.Daylight = ~((inputTable.NightTotal_10Modes_ > 0) | (isnan(inputTable.SpeedV85)));
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

function inputTable = applyUptimeCorrections(inputTable, config)
    uptimeCorrection = 1./inputTable.Uptime;
    uptimeCorrection(uptimeCorrection > config.analysis.maxUptimeCorrection) = config.analysis.maxUptimeCorrection;
    inputTable.AdjustedCountsUptime = inputTable.(config.analysis.modeString) .* uptimeCorrection;
end

function inputTable = applyDaylightCorrections(inputTable, config)
    % Split into weekday/weekend
    inputTableWD = inputTable(inputTable.isWeekday,:);
    inputTableWE = inputTable(~inputTable.isWeekday,:);
    
    % Apply daylight corrections
    inputTableWD.AdjustedCountsUptimeDaylight = inputTableWD.AdjustedCountsUptime .* config.analysis.daylightCorrectionRatioWD;
    inputTableWE.AdjustedCountsUptimeDaylight = inputTableWE.AdjustedCountsUptime .* config.analysis.daylightCorrectionRatioWE;
    
    % Recombine
    inputTable = sortrows([inputTableWD; inputTableWE], 'Date');
end

%% Weather Data Functions
function weatherData = getWeatherForAnalysis(locationData, config)
    % Get unique days across all locations
    allDates = [];
    locationNames = fieldnames(locationData);
    for i = 1:length(locationNames)
        dates = locationData.(locationNames{i}).Date;
        allDates = [allDates; dates];
    end
    
    uniqueDays = unique(dateshift(allDates, 'start', 'day'));
    dailyNoonTimes = uniqueDays + hours(12);
    
    % Get weather data
    fprintf('Getting weather data for %d days...\n', length(uniqueDays));
    [precipitationData, temperatureData, sunriseData, sunsetData, sunhoursData, snowData, windspeedData, feelslikeData] = ...
        getWeatherstackData('Montreal', dailyNoonTimes);
    
    weatherData = struct( ...
        'dates', uniqueDays, ...
        'precipitation', precipitationData, ...
        'temperature', temperatureData, ...
        'feelslike', feelslikeData, ...
        'windspeed', windspeedData, ...
        'sunrise', sunriseData, ...
        'sunset', sunsetData, ...
        'sunhours', sunhoursData, ...
        'snow', snowData ...
    );
end

%% Plotting Functions
function generateConfiguredPlots(locationData, weatherData, config)
    if config.plots.combinedPlots
        generateCombinedPlots(locationData, weatherData, config);
    else
        generateSeparatePlots(locationData, weatherData, config);
    end
end

function generateCombinedPlots(locationData, weatherData, config)
    locationNames = fieldnames(locationData);
    
    for plotType = config.plots.plotTypes
        switch plotType{1}
            case 'daily'
                plotCombinedDaily(locationData, weatherData, config);
            case 'weekly'
                plotCombinedWeekly(locationData, weatherData, config);
            case 'monthly'
                plotCombinedMonthly(locationData, weatherData, config);
        end
    end
end

function plotCombinedDaily(locationData, weatherData, config)
    figure('Position', [408 126 1132 921]);
    
    locationNames = fieldnames(locationData);
    legendEntries = {};
    plotHandles = [];
    
    % Plot weather on right axis if enabled
    if config.plots.showWeather
        yyaxis right
        weatherHandle = plotWeatherData(weatherData, config);
    end
    
    % Plot traffic data on left axis
    yyaxis left
    hold on
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        data = locationData.(locationName);
        locationInfo = data.locationInfo;
        
        % Calculate daily totals
        dailyData = calculateDailyTotals(data, config);
        
        % Plot raw counts if enabled
        if config.plots.showRawCounts
            h1 = plot(dailyData.dates, dailyData.rawCounts, '-', ...
                'LineWidth', config.style.plotLineWidth * 0.5, ...
                'Color', locationInfo.plotColor, ...
                'DisplayName', sprintf('%s Raw (Total: %s)', ...
                    locationInfo.name, num2sepstr(sum(dailyData.rawCounts), '%.0f')));
            plotHandles = [plotHandles, h1];
        end
        
        % Plot adjusted counts if enabled
        if config.plots.showAdjustedCounts
            h2 = plot(dailyData.dates, dailyData.adjustedCounts, '--', ...
                'LineWidth', config.style.plotLineWidth * 0.3, ...
                'Color', locationInfo.plotColor * 0.7, ...
                'DisplayName', sprintf('%s Adjusted (Total: %s)', ...
                    locationInfo.name, num2sepstr(sum(dailyData.adjustedCounts), '%.0f')));
            plotHandles = [plotHandles, h2];
        end
    end
    
    % Format plot
    formatCombinedPlot(config, 'Daily', weatherData);
    legend(plotHandles, 'Location', 'southwest');
    
    hold off
end

function dailyData = calculateDailyTotals(data, config)
    % Create daily grouping
    data.DayOnly = dateshift(data.Date, 'start', 'day');
    
    % Calculate raw daily totals
    groupedData = groupsummary(data, 'DayOnly', 'sum', config.analysis.modeString);
    
    % Calculate adjusted daily totals (for times up to cutoff)
    truncatedData = data(timeofday(data.Date) <= config.analysis.truncationCutoffTime, :);
    adjustedGrouped = groupsummary(truncatedData, 'DayOnly', 'sum', 'AdjustedCountsUptimeDaylight');
    
    % Combine results
    dailyData = struct();
    dailyData.dates = groupedData.DayOnly;
    dailyData.rawCounts = groupedData.(['sum_' config.analysis.modeString]);
    
    % Match adjusted counts to raw count dates
    [~, ia, ib] = intersect(groupedData.DayOnly, adjustedGrouped.DayOnly);
    adjustedCounts = nan(size(dailyData.rawCounts));
    adjustedCounts(ia) = adjustedGrouped.sum_AdjustedCountsUptimeDaylight(ib);
    
    % Ensure adjusted is at least as large as raw
    dailyData.adjustedCounts = max(dailyData.rawCounts, adjustedCounts);
end

function weatherHandles = plotWeatherData(weatherData, config)
    if ~config.plots.showWeather
        weatherHandles = [];
        return;
    end
    
    % Plot temperature line
    h1 = plot(weatherData.dates, weatherData.feelslike, '-', ...
        'LineWidth', config.style.plotLineWidth * 0.5, ...
        'Color', [0 0.4471 0.7412 0.3], ...
        'DisplayName', 'Feels Like Temperature (°C)');
    
    weatherHandles = h1;
    
    % Add precipitation bubbles if enabled
    if config.plots.showPrecipitationBubbles
        hold on
        h2 = bubblechart(weatherData.dates, weatherData.feelslike, weatherData.precipitation, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.3, ...
            'DisplayName', 'Precipitation (bubble size)');
        weatherHandles = [weatherHandles, h2];
        hold off
    end
    
    ylabel('Temperature (°C)', 'FontSize', config.style.labelFontSize, ...
        'Color', [1 0 0 0.5] + 0.5.*[0 1 1 0], 'FontWeight', 'bold');
end

function formatCombinedPlot(config, plotType, weatherData)
    % Left axis formatting
    yyaxis left
    ylabel(['Total ' config.analysis.modeDisplayString ' per Day'], ...
        'FontSize', config.style.labelFontSize + 2, 'FontWeight', 'bold');
    
    ax = gca;
    ax.FontSize = config.style.axisFontSize;
    ax.YAxis(1).Color = 'k';
    if config.plots.showWeather
        ax.YAxis(2).Color = [0 0.4471 0.7412];
    end
    
    % Title and formatting
    title([plotType ' ' config.analysis.modeDisplayString ' Comparison'], ...
        'FontSize', config.style.titleFontSize);
    
    xlabel('Date', 'FontSize', config.style.labelFontSize);
    set(gca, 'Color', config.style.axisBackgroundColor);
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

%% Utility Functions (keeping existing functions like num2sepstr, getWeatherstackData, etc.)

% [Include the existing utility functions here - num2sepstr, getWeatherstackData, etc.]