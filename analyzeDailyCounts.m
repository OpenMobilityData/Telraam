%% Analyze Daily Count Patterns
% Check what's actually happening with daily counts

clear all; clc;

% Load Eastern 2025 data as example
filename = 'raw-data-2025East15.xlsx';
data = readtable(filename);

% Convert dates
datesRaw = data.DateAndTime_Local_;
if iscell(datesRaw)
    dates = datetime(datesRaw, 'InputFormat', 'yyyy-MM-dd HH:mm');
else
    dates = datesRaw;
end

bikeCounts = data.BikeTotal;

% Define AM window (6:30-9:30)
windowAM = [6, 30, 9, 30];

% Get unique days
uniqueDays = unique(dateshift(dates, 'start', 'day'));

% Calculate for each day
fprintf('Day-by-day analysis of AM bike counts:\n');
fprintf('=======================================\n');
fprintf('Date\t\t\tIntervals\tTotal\tExpected\tCoverage\n');

zeroCountDays = [];
lowCountDays = [];

for dayIdx = 1:min(30, length(uniqueDays))  % First 30 days
    currentDay = uniqueDays(dayIdx);
    
    % Create window times
    windowStart = currentDay + hours(windowAM(1)) + minutes(windowAM(2));
    windowEnd = currentDay + hours(windowAM(3)) + minutes(windowAM(4));
    
    % Find data in window
    inWindow = dates >= windowStart & dates <= windowEnd;
    
    % Calculate stats
    intervalsFound = sum(inWindow);
    expectedIntervals = 12;  % 3 hours * 4 intervals per hour
    totalCount = sum(bikeCounts(inWindow));
    coverage = 100 * intervalsFound / expectedIntervals;
    
    fprintf('%s\t%d\t\t%d\t\t%d\t\t%.0f%%\n', ...
        datestr(currentDay, 'yyyy-mm-dd'), intervalsFound, totalCount, ...
        expectedIntervals, coverage);
    
    if totalCount == 0
        zeroCountDays(end+1) = dayIdx;
    elseif totalCount < 10
        lowCountDays(end+1) = dayIdx;
    end
end

fprintf('\nZero-count days: %d\n', length(zeroCountDays));
fprintf('Low-count days (<10): %d\n', length(lowCountDays));

% Now check a specific problematic period (around the March gap)
fprintf('\n\nChecking around March 2025 gap:\n');
fprintf('================================\n');

marchStart = datetime(2025, 3, 1);
marchEnd = datetime(2025, 4, 15);
marchDays = uniqueDays(uniqueDays >= marchStart & uniqueDays <= marchEnd);

for dayIdx = 1:length(marchDays)
    currentDay = marchDays(dayIdx);
    
    % Check if ANY data exists for this day
    dayData = dates >= currentDay & dates < currentDay + days(1);
    
    if sum(dayData) == 0
        fprintf('%s: NO DATA AT ALL\n', datestr(currentDay, 'yyyy-mm-dd'));
    else
        % Check AM window
        windowStart = currentDay + hours(windowAM(1)) + minutes(windowAM(2));
        windowEnd = currentDay + hours(windowAM(3)) + minutes(windowAM(4));
        inWindow = dates >= windowStart & dates <= windowEnd;
        
        fprintf('%s: %d intervals in AM window, total count = %d\n', ...
            datestr(currentDay, 'yyyy-mm-dd'), sum(inWindow), ...
            sum(bikeCounts(inWindow)));
    end
end