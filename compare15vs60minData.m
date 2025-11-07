%% Compare Daily Totals: 15-minute vs 60-minute Resolution Data
% This script compares daily count totals from both data sources
% to identify and visualize discrepancies in the 15-minute export

clear all; close all; clc;

fprintf('================================================\n');
fprintf('Daily Totals Comparison: 15-min vs 60-min Data\n');
fprintf('================================================\n\n');

%% Configuration
locations = {
    struct('name', 'Eastern', 'file15min_2024', 'raw-data-2024East15.xlsx', ...
           'file60min_2024', 'raw-data-2024East60.xlsx', ...
           'file15min_2025', 'raw-data-2025East15.xlsx', ...
           'file60min_2025', 'raw-data-2025East60.xlsx');
    struct('name', 'Western', 'file15min_2024', 'raw-data-2024West15.xlsx', ...
           'file60min_2024', 'raw-data-2024West60.xlsx', ...
           'file15min_2025', 'raw-data-2025West15.xlsx', ...
           'file60min_2025', 'raw-data-2025West60.xlsx');
};

% Process each location
for locIdx = 1:length(locations)
    location = locations{locIdx};
    
    fprintf('\n========================================\n');
    fprintf('Processing %s Location\n', location.name);
    fprintf('========================================\n');
    
    %% Load and process 15-minute data
    fprintf('\nLoading 15-minute resolution data...\n');
    data15min_2024 = loadAndProcessData(location.file15min_2024, '15-min 2024');
    data15min_2025 = loadAndProcessData(location.file15min_2025, '15-min 2025');
    
    % Combine 15-minute data
    if ~isempty(data15min_2024) && ~isempty(data15min_2025)
        data15min = [data15min_2024; data15min_2025];
    elseif ~isempty(data15min_2024)
        data15min = data15min_2024;
    elseif ~isempty(data15min_2025)
        data15min = data15min_2025;
    else
        data15min = table();
    end
    
    %% Load and process 60-minute data
    fprintf('\nLoading 60-minute resolution data...\n');
    data60min_2024 = loadAndProcessData(location.file60min_2024, '60-min 2024');
    data60min_2025 = loadAndProcessData(location.file60min_2025, '60-min 2025');
    
    % Combine 60-minute data
    if ~isempty(data60min_2024) && ~isempty(data60min_2025)
        data60min = [data60min_2024; data60min_2025];
    elseif ~isempty(data60min_2024)
        data60min = data60min_2024;
    elseif ~isempty(data60min_2025)
        data60min = data60min_2025;
    else
        data60min = table();
    end
    
    %% Calculate daily totals for both datasets
    fprintf('\nCalculating daily totals...\n');
    
    % For bikes
    dailyBikes15min = calculateDailyTotals(data15min, 'BikeTotal');
    dailyBikes60min = calculateDailyTotals(data60min, 'BikeTotal');
    
    % For cars
    dailyCars15min = calculateDailyTotals(data15min, 'CarTotal');
    dailyCars60min = calculateDailyTotals(data60min, 'CarTotal');
    
    %% Create comparison plots
    
    % BIKES comparison
    fig1 = figure('Position', [100, 100, 1400, 800]);
    
    % Plot both datasets
    subplot(3, 1, 1);
    hold on;
    plot(dailyBikes60min.dates, dailyBikes60min.counts, 'b-', 'LineWidth', 1.5, 'DisplayName', '60-min data');
    plot(dailyBikes15min.dates, dailyBikes15min.counts, 'r-', 'LineWidth', 1.5, 'DisplayName', '15-min data');
    hold off;
    
    xlabel('Date');
    ylabel('Daily Bike Count');
    title(sprintf('%s Location - Daily Bike Totals Comparison', location.name));
    legend('Location', 'best');
    grid on;
    
    % Plot difference
    subplot(3, 1, 2);
    [commonDates, idx15, idx60] = intersect(dailyBikes15min.dates, dailyBikes60min.dates);
    if ~isempty(commonDates)
        difference = dailyBikes60min.counts(idx60) - dailyBikes15min.counts(idx15);
        bar(commonDates, difference, 'FaceColor', [0.8 0.2 0.2]);
        xlabel('Date');
        ylabel('Difference (60min - 15min)');
        title('Daily Count Differences');
        grid on;
        
        % Add statistics text
        meanDiff = mean(difference);
        maxDiff = max(abs(difference));
        percentDaysWithDiff = 100 * sum(abs(difference) > 5) / length(difference);
        
        text(0.02, 0.98, sprintf('Mean difference: %.1f bikes\nMax difference: %.0f bikes\n%.1f%% of days differ by >5 bikes', ...
            meanDiff, maxDiff, percentDaysWithDiff), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    % Plot scatter comparison
    subplot(3, 1, 3);
    if ~isempty(commonDates)
        scatter(dailyBikes15min.counts(idx15), dailyBikes60min.counts(idx60), 20, 'filled');
        hold on;
        % Add 1:1 reference line
        maxVal = max([dailyBikes15min.counts(idx15); dailyBikes60min.counts(idx60)]);
        plot([0 maxVal], [0 maxVal], 'r--', 'LineWidth', 1);
        hold off;
        
        xlabel('15-min Data Daily Count');
        ylabel('60-min Data Daily Count');
        title('Correlation Plot');
        grid on;
        
        % Calculate correlation
        validIdx = ~isnan(dailyBikes15min.counts(idx15)) & ~isnan(dailyBikes60min.counts(idx60));
        if sum(validIdx) > 2
            R = corrcoef(dailyBikes15min.counts(idx15(validIdx)), dailyBikes60min.counts(idx60(validIdx)));
            text(0.02, 0.98, sprintf('Correlation: %.3f', R(1,2)), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
    end
    
    sgtitle(sprintf('%s Location - BIKES - 15-min vs 60-min Data Comparison', location.name));
    
    % Save figure name
    set(fig1, 'Name', sprintf('%s_Bikes_Comparison', location.name), 'NumberTitle', 'off');
    
    % CARS comparison
    fig2 = figure('Position', [100, 100, 1400, 800]);
    
    % Plot both datasets
    subplot(3, 1, 1);
    hold on;
    plot(dailyCars60min.dates, dailyCars60min.counts, 'b-', 'LineWidth', 1.5, 'DisplayName', '60-min data');
    plot(dailyCars15min.dates, dailyCars15min.counts, 'r-', 'LineWidth', 1.5, 'DisplayName', '15-min data');
    hold off;
    
    xlabel('Date');
    ylabel('Daily Car Count');
    title(sprintf('%s Location - Daily Car Totals Comparison', location.name));
    legend('Location', 'best');
    grid on;
    
    % Plot difference
    subplot(3, 1, 2);
    [commonDates, idx15, idx60] = intersect(dailyCars15min.dates, dailyCars60min.dates);
    if ~isempty(commonDates)
        difference = dailyCars60min.counts(idx60) - dailyCars15min.counts(idx15);
        bar(commonDates, difference, 'FaceColor', [0.8 0.2 0.2]);
        xlabel('Date');
        ylabel('Difference (60min - 15min)');
        title('Daily Count Differences');
        grid on;
        
        % Add statistics text
        meanDiff = mean(difference);
        maxDiff = max(abs(difference));
        percentDaysWithDiff = 100 * sum(abs(difference) > 10) / length(difference);
        
        text(0.02, 0.98, sprintf('Mean difference: %.1f cars\nMax difference: %.0f cars\n%.1f%% of days differ by >10 cars', ...
            meanDiff, maxDiff, percentDaysWithDiff), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    % Plot scatter comparison
    subplot(3, 1, 3);
    if ~isempty(commonDates)
        scatter(dailyCars15min.counts(idx15), dailyCars60min.counts(idx60), 20, 'filled');
        hold on;
        % Add 1:1 reference line
        maxVal = max([dailyCars15min.counts(idx15); dailyCars60min.counts(idx60)]);
        plot([0 maxVal], [0 maxVal], 'r--', 'LineWidth', 1);
        hold off;
        
        xlabel('15-min Data Daily Count');
        ylabel('60-min Data Daily Count');
        title('Correlation Plot');
        grid on;
        
        % Calculate correlation
        validIdx = ~isnan(dailyCars15min.counts(idx15)) & ~isnan(dailyCars60min.counts(idx60));
        if sum(validIdx) > 2
            R = corrcoef(dailyCars15min.counts(idx15(validIdx)), dailyCars60min.counts(idx60(validIdx)));
            text(0.02, 0.98, sprintf('Correlation: %.3f', R(1,2)), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
    end
    
    sgtitle(sprintf('%s Location - CARS - 15-min vs 60-min Data Comparison', location.name));
    
    % Save figure name
    set(fig2, 'Name', sprintf('%s_Cars_Comparison', location.name), 'NumberTitle', 'off');
    
    %% Generate discrepancy report
    fprintf('\n--- Discrepancy Report for %s Location ---\n', location.name);
    
    % Find days with large discrepancies for bikes
    [commonDates, idx15, idx60] = intersect(dailyBikes15min.dates, dailyBikes60min.dates);
    if ~isempty(commonDates)
        bikeDiff = dailyBikes60min.counts(idx60) - dailyBikes15min.counts(idx15);
        largeDiffIdx = find(abs(bikeDiff) > 50);
        
        if ~isempty(largeDiffIdx)
            fprintf('\nDays with >50 bike count difference:\n');
            fprintf('Date\t\t\t15-min\t60-min\tDifference\n');
            for i = 1:min(10, length(largeDiffIdx))
                idx = largeDiffIdx(i);
                fprintf('%s\t%d\t%d\t%+d\n', ...
                    datestr(commonDates(idx), 'yyyy-mm-dd'), ...
                    dailyBikes15min.counts(idx15(idx)), ...
                    dailyBikes60min.counts(idx60(idx)), ...
                    bikeDiff(idx));
            end
            if length(largeDiffIdx) > 10
                fprintf('... and %d more days\n', length(largeDiffIdx) - 10);
            end
        end
        
        % Summary statistics
        fprintf('\nBIKES Summary:\n');
        fprintf('Total days compared: %d\n', length(commonDates));
        fprintf('Mean absolute difference: %.1f bikes\n', mean(abs(bikeDiff)));
        fprintf('Days with >20 bike difference: %d (%.1f%%)\n', ...
            sum(abs(bikeDiff) > 20), 100*sum(abs(bikeDiff) > 20)/length(bikeDiff));
        fprintf('Days with >50 bike difference: %d (%.1f%%)\n', ...
            sum(abs(bikeDiff) > 50), 100*sum(abs(bikeDiff) > 50)/length(bikeDiff));
        fprintf('Days with >100 bike difference: %d (%.1f%%)\n', ...
            sum(abs(bikeDiff) > 100), 100*sum(abs(bikeDiff) > 100)/length(bikeDiff));
    end
    
    % Similar analysis for cars
    [commonDates, idx15, idx60] = intersect(dailyCars15min.dates, dailyCars60min.dates);
    if ~isempty(commonDates)
        carDiff = dailyCars60min.counts(idx60) - dailyCars15min.counts(idx15);
        
        fprintf('\nCARS Summary:\n');
        fprintf('Total days compared: %d\n', length(commonDates));
        fprintf('Mean absolute difference: %.1f cars\n', mean(abs(carDiff)));
        fprintf('Days with >50 car difference: %d (%.1f%%)\n', ...
            sum(abs(carDiff) > 50), 100*sum(abs(carDiff) > 50)/length(carDiff));
        fprintf('Days with >100 car difference: %d (%.1f%%)\n', ...
            sum(abs(carDiff) > 100), 100*sum(abs(carDiff) > 100)/length(carDiff));
    end
end

fprintf('\n================================================\n');
fprintf('Comparison Complete\n');
fprintf('================================================\n');

%% Helper Functions

function data = loadAndProcessData(filename, description)
    if exist(filename, 'file')
        fprintf('  Loading %s: %s\n', description, filename);
        data = readtable(filename);
    else
        fprintf('  File not found: %s\n', filename);
        data = table();
    end
end

function dailyData = calculateDailyTotals(data, modeColumn)
    dailyData = struct('dates', [], 'counts', []);
    
    if isempty(data)
        return;
    end
    
    % Check if mode column exists
    if ~any(strcmp(data.Properties.VariableNames, modeColumn))
        return;
    end
    
    % Get dates - handle different possible column names
    dateColumn = '';
    possibleDateCols = {'DateAndTime_Local_', 'Date_x0020_and_x0020_Time_x0020__x0028_Local_x0029_'};
    for i = 1:length(possibleDateCols)
        if any(strcmp(data.Properties.VariableNames, possibleDateCols{i}))
            dateColumn = possibleDateCols{i};
            break;
        end
    end
    
    if isempty(dateColumn)
        warning('Date column not found');
        return;
    end
    
    % Extract dates and counts
    datesRaw = data.(dateColumn);
    counts = data.(modeColumn);
    
    % Convert dates if needed
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
    
    % Calculate daily totals
    uniqueDays = unique(dateshift(dates, 'start', 'day'));
    dailyCounts = zeros(length(uniqueDays), 1);
    
    for i = 1:length(uniqueDays)
        dayMask = dateshift(dates, 'start', 'day') == uniqueDays(i);
        dailyCounts(i) = sum(counts(dayMask), 'omitnan');
    end
    
    dailyData.dates = uniqueDays;
    dailyData.counts = dailyCounts;
end