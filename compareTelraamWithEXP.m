%% Traffic Count Comparison with EXP Study - CORRECTED VERSION
% This script generates year-over-year daily traffic count comparisons
% using 15-minute resolution data and compares with EXP Study survey counts

clear all; close all; clc;

%% Configuration Setup

% Location definitions
westernSegmentName = 'rue de Terrebonne @ King Edward';
easternSegmentName = 'rue de Terrebonne @ Draper';

locations = {
    struct('name', easternSegmentName, ...
           'fileStem2024', 'raw-data-2024East15', ...  % 15-minute resolution files
           'fileStem2025', 'raw-data-2025East15', ...
           'shortName', 'Eastern', ...
           'plotColor', [0 0 1]);  % Blue
    struct('name', westernSegmentName, ...
           'fileStem2024', 'raw-data-2024West15', ...  % 15-minute resolution files
           'fileStem2025', 'raw-data-2025West15', ...
           'shortName', 'Western', ...
           'plotColor', [0 0 0]);  % Black
};

% EXP Study reference data
expStudyDate = datetime(2024, 10, 4);  % October 4 (adjust year as needed)

% Define time windows and survey counts for each modality
expStudyData = struct();

% Bikes
expStudyData.bikes = struct();
expStudyData.bikes.timeWindows = struct( ...
    'AM', [6, 30, 9, 30], ...  % [startHour, startMin, endHour, endMin]
    'PM', [15, 45, 16, 45] ...
);
expStudyData.bikes.counts = struct( ...
    'Western', struct('AM', 84, 'PM', 81, 'Total', 165), ...
    'Eastern', struct('AM', 104, 'PM', 72, 'Total', 176) ...
);

% Cars
expStudyData.cars = struct();
expStudyData.cars.timeWindows = struct( ...
    'AM', [7, 45, 8, 45], ...
    'PM', [15, 45, 16, 45] ...
);
expStudyData.cars.counts = struct( ...
    'Western', struct('AM', 50, 'PM', 55, 'Total', 105), ...
    'Eastern', struct('AM', 321, 'PM', 264, 'Total', 585) ...
);

% Analysis parameters
analysis = struct( ...
    'startTime', datetime(2024, 08, 15, 00, 00, 01), ...
    'endTime', datetime(2025, 10, 31, 23, 59, 59), ...
    'uptimeThreshold', 0.0, ...
    'maxUptimeCorrection', 1.0 ...
);

% Plot style parameters
style = struct( ...
    'plotLineWidth', 2.0, ...
    'axisFontSize', 12.0, ...
    'labelFontSize', 14.0, ...
    'titleFontSize', 16.0, ...
    'legendFontSize', 10.0, ...
    'axisBackgroundColor', 0.95.*[1 1 1], ...
    'expLineColor', [0.8 0 0.8], ...  % Purple for EXP count line (distinct)
    'meanLineColor', [0 0.6 0], ...  % Green for mean line
    'markerSize', 8 ...
);

%% Process Each Location and Modality

% Define modalities to analyze - CORRECTED COLUMN NAMES
modalities = {
    struct('name', 'BikeTotal', 'displayName', 'Bikes', 'expData', expStudyData.bikes),
    struct('name', 'CarTotal', 'displayName', 'Cars', 'expData', expStudyData.cars)
};

% Process each location
for locIdx = 1:length(locations)
    location = locations{locIdx};
    fprintf('\n========================================\n');
    fprintf('Processing %s location\n', location.shortName);
    fprintf('========================================\n');
    
    % Load data for this location
    locationData = loadLocationData15Min(location, analysis);
    
    % Process each modality
    for modIdx = 1:length(modalities)
        modality = modalities{modIdx};
        fprintf('\nAnalyzing %s...\n', modality.displayName);
        
        % Generate plots for AM, PM, and Combined
        timeperiods = {'AM', 'PM', 'Combined'};
        
        for tpIdx = 1:length(timeperiods)
            timeperiod = timeperiods{tpIdx};
            
            % Create the year-over-year comparison plot
            plotYearOverYearWithEXP(locationData, location, modality, ...
                timeperiod, expStudyData, expStudyDate, analysis, style);
        end
    end
end

fprintf('\n========================================\n');
fprintf('Analysis complete. All plots generated.\n');
fprintf('========================================\n');

%% Helper Functions

function data = loadLocationData15Min(location, analysis)
    % Load 15-minute resolution data for a location
    
    data = struct();
    
    % Load 2024 data
    filename2024 = [location.fileStem2024 '.xlsx'];
    if exist(filename2024, 'file')
        fprintf('Loading 2024 data from %s...\n', filename2024);
        data.data2024 = readtable(filename2024);
    else
        warning('2024 data file not found: %s', filename2024);
        data.data2024 = table();
    end
    
    % Load 2025 data
    filename2025 = [location.fileStem2025 '.xlsx'];
    if exist(filename2025, 'file')
        fprintf('Loading 2025 data from %s...\n', filename2025);
        data.data2025 = readtable(filename2025);
    else
        warning('2025 data file not found: %s', filename2025);
        data.data2025 = table();
    end
    
    % Combine data
    if ~isempty(data.data2024) && ~isempty(data.data2025)
        data.combined = [data.data2024; data.data2025];
    elseif ~isempty(data.data2024)
        data.combined = data.data2024;
    elseif ~isempty(data.data2025)
        data.combined = data.data2025;
    else
        data.combined = table();
    end
    
    % Filter by analysis time range if data exists
    if ~isempty(data.combined)
        % Convert dates if needed and filter
        dateColumn = 'DateAndTime_Local_';
        if any(strcmp(data.combined.Properties.VariableNames, dateColumn))
            % Convert cell array dates to datetime
            datesRaw = data.combined.(dateColumn);
            if iscell(datesRaw)
                dates = datetime(datesRaw, 'InputFormat', 'yyyy-MM-dd HH:mm');
            else
                dates = datesRaw;
            end
            
            inRange = dates >= analysis.startTime & dates <= analysis.endTime;
            data.combined = data.combined(inRange, :);
            fprintf('Data filtered to analysis period: %d rows\n', sum(inRange));
        end
    end
end

function plotYearOverYearWithEXP(locationData, location, modality, timeperiod, expStudyData, expStudyDate, analysis, style)
    % Create year-over-year comparison plot with EXP study reference lines
    
    % Extract time windows based on modality and time period
    if strcmp(timeperiod, 'AM')
        timeWindow = modality.expData.timeWindows.AM;
        expCount = modality.expData.counts.(location.shortName).AM;
    elseif strcmp(timeperiod, 'PM')
        timeWindow = modality.expData.timeWindows.PM;
        expCount = modality.expData.counts.(location.shortName).PM;
    else % Combined
        % For combined, we'll use both AM and PM windows
        timeWindowAM = modality.expData.timeWindows.AM;
        timeWindowPM = modality.expData.timeWindows.PM;
        expCount = modality.expData.counts.(location.shortName).Total;
    end
    
    % Calculate daily counts for the specified time windows
    if strcmp(timeperiod, 'Combined')
        dailyData = calculateDailyCountsForTimeWindows(locationData.combined, ...
            modality.name, {timeWindowAM, timeWindowPM});
    else
        dailyData = calculateDailyCountsForTimeWindows(locationData.combined, ...
            modality.name, {timeWindow});
    end
    
    if isempty(dailyData.dates)
        warning('No data available for %s %s at %s location', ...
            modality.displayName, timeperiod, location.shortName);
        return;
    end
    
    % Create year periods for comparison
    yearPeriods = createYearPeriods(dailyData);
    
    if isempty(yearPeriods)
        warning('Insufficient data for year-over-year comparison');
        return;
    end
    
    % Create the plot
    fig = figure('Position', [100, 100, 1200, 600]);
    hold on;
    
    % Plot each year period
    legendEntries = {};
    for i = 1:length(yearPeriods)
        period = yearPeriods(i);
        color = getYearColor(i);
        
        plot(period.normalizedDates, period.counts, '-', ...
            'Color', color, 'LineWidth', style.plotLineWidth, ...
            'DisplayName', period.label);
        legendEntries{end+1} = period.label;
    end
    
    % Add EXP Study reference line
    xlims = xlim;
    plot(xlims, [expCount, expCount], '--', ...
        'Color', style.expLineColor, 'LineWidth', 2, ...
        'DisplayName', sprintf('EXP Count (%d)', expCount));
    legendEntries{end+1} = sprintf('EXP Count (%d)', expCount);
    
    % Add marker on EXP Study date
    % Normalize the EXP study date using Aug 15, 2024 as reference
    aug15_2024 = datetime(2024, 8, 15);
    expStudyNormalized = normalizeToReferenceYear(expStudyDate, aug15_2024);
    plot(expStudyNormalized, expCount, 'o', ...
        'MarkerSize', style.markerSize, ...
        'MarkerEdgeColor', style.expLineColor, ...
        'MarkerFaceColor', style.expLineColor, ...
        'HandleVisibility', 'off');
    
    % Calculate and add mean line
    meanCount = mean(dailyData.counts, 'omitnan');
    plot(xlims, [meanCount, meanCount], '--', ...
        'Color', style.meanLineColor, 'LineWidth', 2, ...
        'DisplayName', sprintf('Year Mean (%.1f)', meanCount));
    legendEntries{end+1} = sprintf('Year Mean (%.1f)', meanCount);
    
    % Format the plot
    xlabel('Date', 'FontSize', style.labelFontSize);
    ylabel(sprintf('%s Count', modality.displayName), 'FontSize', style.labelFontSize);
    
    titleStr = sprintf('%s - %s %s Counts\nYear-over-Year Comparison with EXP Study', ...
        location.shortName, timeperiod, modality.displayName);
    title(titleStr, 'FontSize', style.titleFontSize);
    
    % Set x-axis to show months
    ax = gca;
    ax.XAxis.TickLabelFormat = 'MMM';
    ax.FontSize = style.axisFontSize;
    
    % Add grid
    grid on;
    ax.GridAlpha = 0.3;
    ax.Color = style.axisBackgroundColor;
    
    % Add legend
    legend(legendEntries, 'Location', 'best', 'FontSize', style.legendFontSize);
    
    % Add time window information as text
    if strcmp(timeperiod, 'Combined')
        windowText = sprintf('AM: %02d:%02d-%02d:%02d, PM: %02d:%02d-%02d:%02d', ...
            timeWindowAM(1), timeWindowAM(2), timeWindowAM(3), timeWindowAM(4), ...
            timeWindowPM(1), timeWindowPM(2), timeWindowPM(3), timeWindowPM(4));
    else
        windowText = sprintf('%s Window: %02d:%02d-%02d:%02d', ...
            timeperiod, timeWindow(1), timeWindow(2), timeWindow(3), timeWindow(4));
    end
    
    % Add text annotation
    text(0.02, 0.98, windowText, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    hold off;
    
    % Display the figure name in the title bar
    figName = sprintf('%s_%s_%s_YoY_Comparison', ...
        location.shortName, modality.displayName, timeperiod);
    set(fig, 'Name', figName, 'NumberTitle', 'off');
    
    % Note: Figure saving removed - save manually if needed
    fprintf('Generated plot: %s\n', figName);
end

function dailyData = calculateDailyCountsForTimeWindows(data, modeString, timeWindows)
    % Calculate daily counts for specific time windows
    
    dailyData = struct('dates', [], 'counts', []);
    
    if isempty(data)
        return;
    end
    
    % Check if the mode column exists
    if ~any(strcmp(data.Properties.VariableNames, modeString))
        warning('Mode column "%s" not found in data', modeString);
        fprintf('Available columns: %s\n', strjoin(data.Properties.VariableNames, ', '));
        return;
    end
    
    % Get dates and counts
    dateColumn = 'DateAndTime_Local_';
    if ~any(strcmp(data.Properties.VariableNames, dateColumn))
        warning('Date column not found');
        return;
    end
    
    % Extract dates and counts
    datesRaw = data.(dateColumn);
    counts = data.(modeString);
    
    % Convert dates if they are cell arrays (from Excel)
    if iscell(datesRaw)
        try
            dates = datetime(datesRaw, 'InputFormat', 'yyyy-MM-dd HH:mm');
        catch
            try
                dates = datetime(datesRaw);
            catch
                warning('Could not parse date format');
                return;
            end
        end
    else
        dates = datesRaw;
    end
    
    % Remove any NaN or invalid data
    validIdx = ~isnat(dates) & ~isnan(counts);
    dates = dates(validIdx);
    counts = counts(validIdx);
    
    if isempty(dates)
        warning('No valid data after filtering');
        return;
    end
    
    % Get unique days
    uniqueDays = unique(dateshift(dates, 'start', 'day'));
    
    dailyCounts = zeros(length(uniqueDays), 1);
    
    % For each day
    for dayIdx = 1:length(uniqueDays)
        currentDay = uniqueDays(dayIdx);
        dayTotal = 0;
        
        % For each time window
        for winIdx = 1:length(timeWindows)
            window = timeWindows{winIdx};
            
            % Create start and end times for this window
            windowStart = currentDay + hours(window(1)) + minutes(window(2));
            windowEnd = currentDay + hours(window(3)) + minutes(window(4));
            
            % Find data points within this window
            inWindow = dates >= windowStart & dates <= windowEnd;
            
            % Sum counts in this window
            dayTotal = dayTotal + sum(counts(inWindow));
        end
        
        dailyCounts(dayIdx) = dayTotal;
    end
    
    dailyData.dates = uniqueDays;
    dailyData.counts = dailyCounts;
    
    fprintf('Calculated daily counts: %d days, mean=%.1f, max=%d\n', ...
        length(uniqueDays), mean(dailyCounts), max(dailyCounts));
end

function yearPeriods = createYearPeriods(dailyData)
    % Create year periods for comparison (Aug 15 to Aug 14)
    
    yearPeriods = [];
    
    if isempty(dailyData.dates)
        return;
    end
    
    % Sort data
    [sortedDates, sortIdx] = sort(dailyData.dates);
    sortedCounts = dailyData.counts(sortIdx);
    
    % Find data range
    minDate = min(sortedDates);
    maxDate = max(sortedDates);
    
    % Start from August 15, 2024
    currentAug15 = datetime(2024, 8, 15);
    periodCount = 0;
    
    while currentAug15 <= maxDate
        periodStart = currentAug15;
        periodEnd = currentAug15 + years(1) - days(1);
        
        % Find data within this period
        inPeriod = sortedDates >= periodStart & sortedDates <= periodEnd;
        
        if sum(inPeriod) >= 30  % Require at least 30 days of data
            periodCount = periodCount + 1;
            
            periodDates = sortedDates(inPeriod);
            periodCounts = sortedCounts(inPeriod);
            
            % Normalize dates to reference period - pass periodStart as in original
            normalizedDates = normalizeToReferenceYear(periodDates, periodStart);
            
            % No need to sort again since periodDates are already sorted
            % and normalization preserves the order
            
            % Create label
            startYear = year(periodStart);
            endYear = year(periodEnd);
            label = sprintf('%d-%d', startYear, endYear);
            
            yearPeriods(periodCount).originalDates = periodDates;
            yearPeriods(periodCount).normalizedDates = normalizedDates;
            yearPeriods(periodCount).counts = periodCounts;
            yearPeriods(periodCount).label = label;
            yearPeriods(periodCount).periodStart = periodStart;
            yearPeriods(periodCount).periodEnd = periodEnd;
            
            fprintf('Created period %s: %d days of data\n', label, sum(inPeriod));
        end
        
        currentAug15 = currentAug15 + years(1);
    end
end

function normalizedDates = normalizeToReferenceYear(dates, periodStart)
    % Normalize dates to reference period (Aug 15, 2024 to Aug 14, 2025)
    % Using the same logic as the original script
    
    normalizedDates = NaT(size(dates));
    
    % Reference August 15 is always 2024
    refAug15 = datetime(2024, 8, 15);
    
    for i = 1:length(dates)
        % Calculate days since period start
        daysSincePeriodStart = days(dates(i) - periodStart);
        
        % Add these days to reference August 15
        normalizedDates(i) = refAug15 + days(daysSincePeriodStart);
    end
end

function color = getYearColor(index)
    % Get distinct colors for each year
    colors = [
        0 0 1;      % Blue for first year (2024-2025)
        1 0.5 0;    % Orange for second year (2025-2026) - changed from red
        0 0.7 0;    % Green for third year (if needed)
        0.5 0 0.5;  % Dark purple for fourth year (if needed)
        0 0.8 0.8;  % Cyan
        0.8 0.8 0;  % Yellow
    ];
    
    if index <= size(colors, 1)
        color = colors(index, :);
    else
        % Use a colormap for additional colors
        cmap = lines(index);
        color = cmap(index, :);
    end
end