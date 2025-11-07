%% Diagnostic Script - Analyze Data Gaps in 15-minute Resolution Files
% This script investigates gaps in the data to understand missing periods

clear all; clc;

fprintf('================================================\n');
fprintf('Data Gap Analysis for 15-minute Resolution Files\n');
fprintf('================================================\n\n');

%% Load all data files
files = {
    'raw-data-2024East15.xlsx', '2024 Eastern';
    'raw-data-2025East15.xlsx', '2025 Eastern';
    'raw-data-2024West15.xlsx', '2024 Western';
    'raw-data-2025West15.xlsx', '2025 Western';
};

for fileIdx = 1:size(files, 1)
    filename = files{fileIdx, 1};
    description = files{fileIdx, 2};
    
    if ~exist(filename, 'file')
        fprintf('File not found: %s\n\n', filename);
        continue;
    end
    
    fprintf('========================================\n');
    fprintf('Analyzing: %s\n', description);
    fprintf('========================================\n\n');
    
    % Load data
    data = readtable(filename);
    
    % Get dates and bike counts
    datesRaw = data.DateAndTime_Local_;
    bikeCounts = data.BikeTotal;
    
    % Convert dates
    if iscell(datesRaw)
        dates = datetime(datesRaw, 'InputFormat', 'yyyy-MM-dd HH:mm');
    else
        dates = datesRaw;
    end
    
    % Sort by date
    [dates, sortIdx] = sort(dates);
    bikeCounts = bikeCounts(sortIdx);
    
    % Basic statistics
    fprintf('Date range: %s to %s\n', datestr(dates(1)), datestr(dates(end)));
    fprintf('Total records: %d\n', length(dates));
    fprintf('Expected 15-min intervals in range: %d\n', ...
        round(minutes(dates(end) - dates(1)) / 15));
    
    % Find gaps (where time difference > 15 minutes)
    timeDiffs = diff(dates);
    gapThreshold = minutes(16); % Allow 1 minute tolerance
    gaps = find(timeDiffs > gapThreshold);
    
    fprintf('\nData Gaps (> 15 minutes):\n');
    fprintf('-------------------------\n');
    
    if isempty(gaps)
        fprintf('No gaps found!\n');
    else
        fprintf('Found %d gaps:\n\n', length(gaps));
        
        % Analyze gaps
        totalMissingTime = 0;
        gapSizes = [];
        
        for i = 1:min(10, length(gaps)) % Show first 10 gaps
            gapIdx = gaps(i);
            gapStart = dates(gapIdx);
            gapEnd = dates(gapIdx + 1);
            gapDuration = gapEnd - gapStart;
            gapSizes(end+1) = hours(gapDuration);
            totalMissingTime = totalMissingTime + hours(gapDuration);
            
            fprintf('Gap %d: %s to %s (%.1f hours)\n', ...
                i, datestr(gapStart), datestr(gapEnd), hours(gapDuration));
        end
        
        if length(gaps) > 10
            fprintf('... and %d more gaps\n', length(gaps) - 10);
        end
        
        fprintf('\nGap Statistics:\n');
        fprintf('Total gaps: %d\n', length(gaps));
        fprintf('Total missing time: %.1f hours (%.1f days)\n', ...
            sum(hours(timeDiffs(gaps))), sum(days(timeDiffs(gaps))));
        fprintf('Longest gap: %.1f hours\n', max(hours(timeDiffs(gaps))));
        fprintf('Average gap: %.1f hours\n', mean(hours(timeDiffs(gaps))));
    end
    
    % Check for specific time window coverage
    fprintf('\nTime Window Coverage Analysis:\n');
    fprintf('-------------------------------\n');
    
    % Define time windows
    windows = {
        'Bike AM (6:30-9:30)', [6, 30, 9, 30];
        'Bike PM (15:45-16:45)', [15, 45, 16, 45];
        'Car AM (7:45-8:45)', [7, 45, 8, 45];
    };
    
    % Check coverage for each day
    uniqueDays = unique(dateshift(dates, 'start', 'day'));
    
    for winIdx = 1:size(windows, 1)
        winName = windows{winIdx, 1};
        winTime = windows{winIdx, 2};
        
        daysWithData = 0;
        daysWithoutData = 0;
        
        for dayIdx = 1:length(uniqueDays)
            currentDay = uniqueDays(dayIdx);
            winStart = currentDay + hours(winTime(1)) + minutes(winTime(2));
            winEnd = currentDay + hours(winTime(3)) + minutes(winTime(4));
            
            % Check if we have any data in this window
            inWindow = dates >= winStart & dates <= winEnd;
            
            if any(inWindow)
                daysWithData = daysWithData + 1;
            else
                daysWithoutData = daysWithoutData + 1;
            end
        end
        
        fprintf('%s: %d days with data, %d days without (%.1f%% coverage)\n', ...
            winName, daysWithData, daysWithoutData, ...
            100 * daysWithData / (daysWithData + daysWithoutData));
    end
    
    % Check for systematic patterns in missing data
    fprintf('\nSystematic Patterns in Missing Data:\n');
    fprintf('------------------------------------\n');
    
    % Check if certain hours are consistently missing
    hours_of_day = hour(dates);
    for h = 0:23
        count = sum(hours_of_day == h);
        expected = length(uniqueDays) * 4; % 4 15-min periods per hour
        if count < expected * 0.5 % Less than 50% of expected
            fprintf('Hour %02d:00 has low coverage: %d records (%.1f%% of expected)\n', ...
                h, count, 100 * count / expected);
        end
    end
    
    fprintf('\n');
end

fprintf('================================================\n');
fprintf('Analysis Complete\n');
fprintf('================================================\n');