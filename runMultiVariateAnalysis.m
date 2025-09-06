%% Multivariate Analysis of Multi-Modal Transportation Patterns
% This script analyzes seasonal and weather effects across different transportation
% modes (bikes, cars, pedestrians) to understand winter cycling patterns and
% infrastructure needs.
%
% PREREQUISITE: Run runTelraamAnalysis.m first to load required data
%
% Required variables from workspace:
%   - locationData: processed traffic count data for all locations
%   - weatherData: weather data including temperature, precipitation, etc.
%   - analysis: analysis parameters including time range and mode
%   - style: plotting style parameters
%
% Key Questions Addressed:
% 1. Are bikes uniquely abandoned in winter compared to other modes?
% 2. Is cycling impossible for extended winter periods?
% 3. Is there latent demand for improved winter cycling infrastructure?

close all; clc;

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

fprintf('\n========== MULTIVARIATE MULTI-MODAL TRANSPORTATION ANALYSIS ==========\n');
fprintf('Analyzing seasonal patterns across bikes, cars, and pedestrians\n');
fprintf('Focus: Winter cycling viability and infrastructure implications\n\n');

%% 1. Multi-Modal Seasonal Pattern Analysis
fprintf('\n=== ANALYSIS 1: SEASONAL PATTERNS BY MODE ===\n\n');
seasonalAnalysis = analyzeSeasonalPatternsByMode(locationData, weatherData, analysis, style);

%% 2. Winter Resilience Analysis
fprintf('\n=== ANALYSIS 2: WINTER RESILIENCE BY MODE ===\n\n');
resilienceAnalysis = analyzeWinterResilience(locationData, weatherData, analysis, style);

%% 3. Zero-Count Period Analysis
fprintf('\n=== ANALYSIS 3: ZERO-COUNT PERIODS BY MODE ===\n\n');
zeroCountAnalysis = analyzeZeroCountPeriods(locationData, weatherData, analysis, style);

%% 4. Recovery Time Analysis
fprintf('\n=== ANALYSIS 4: POST-STORM RECOVERY PATTERNS ===\n\n');
recoveryAnalysis = analyzeRecoveryPatterns(locationData, weatherData, analysis, style);

%% 5. Infrastructure Impact Analysis
fprintf('\n=== ANALYSIS 5: LOCATION-BASED DIFFERENCES ===\n\n');
infrastructureAnalysis = analyzeLocationDifferences(locationData, weatherData, analysis, style);

%% 6. Latent Demand Analysis
fprintf('\n=== ANALYSIS 6: LATENT DEMAND INDICATORS ===\n\n');
demandAnalysis = analyzeLatentDemand(seasonalAnalysis, resilienceAnalysis, recoveryAnalysis);

%% 7. Generate Summary Report
generateSummaryReport(seasonalAnalysis, resilienceAnalysis, zeroCountAnalysis, ...
                     recoveryAnalysis, infrastructureAnalysis, demandAnalysis);

%% ======================== ANALYSIS FUNCTIONS ========================

function seasonalAnalysis = analyzeSeasonalPatternsByMode(locationData, weatherData, analysis, style)
    % Analyze seasonal patterns for each transportation mode
    
    modes = {'Bike Total', 'Car Total', 'Pedestrian Total'};
    modeNames = {'Bikes', 'Cars', 'Pedestrians'};
    modeColors = {[0 0 1], [1 0 0], [0 0.7 0]};
    
    seasonalAnalysis = struct();
    locationNames = fieldnames(locationData);
    
    % Process each location
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        data = locationData.(locationName);
        locShortName = extractLocationShortName(data.locationInfo.name);
        
        fprintf('Location: %s\n', data.locationInfo.name);
        fprintf('%s\n', repmat('-', 1, 60));
        
        % Initialize storage for this location
        seasonalAnalysis.(locationName) = struct();
        
        % Process each mode
        for modeIdx = 1:length(modes)
            currentMode = modes{modeIdx};
            modeName = modeNames{modeIdx};
            
            % Create temporary analysis structure for this mode
            tempAnalysis = analysis;
            tempAnalysis.modeString = currentMode;
            tempAnalysis.modeDisplayString = modeName;
            
            % Calculate weekly totals
            weeklyData = calculateWeeklyTotals(data, tempAnalysis);
            
            if ~isempty(weeklyData.weekStarts)
                % Calculate seasonal indices
                [seasonalIndex, monthlyAverage] = calculateSeasonalIndex(weeklyData);
                
                % Store results
                seasonalAnalysis.(locationName).(modeName) = struct();
                seasonalAnalysis.(locationName).(modeName).weeklyData = weeklyData;
                seasonalAnalysis.(locationName).(modeName).seasonalIndex = seasonalIndex;
                seasonalAnalysis.(locationName).(modeName).monthlyAverage = monthlyAverage;
                
                % Calculate summer-to-winter ratio
                summerAvg = mean(monthlyAverage(6:8)); % Jun-Aug
                winterAvg = mean(monthlyAverage([1,2,12])); % Dec-Feb
                winterRetention = winterAvg / summerAvg;
                
                seasonalAnalysis.(locationName).(modeName).winterRetention = winterRetention;
                
                fprintf('  %s: Summer avg = %s, Winter avg = %s, Retention = %.1f%%\n', ...
                    modeName, num2sepstr(summerAvg, '%.0f'), num2sepstr(winterAvg, '%.0f'), ...
                    winterRetention * 100);
            end
        end
        fprintf('\n');
    end
    
    % Create comparative visualization
    plotSeasonalComparison(seasonalAnalysis, style);
end

function resilienceAnalysis = analyzeWinterResilience(locationData, weatherData, analysis, style)
    % Analyze how each mode responds to winter conditions
    
    modes = {'Bike Total', 'Car Total', 'Pedestrian Total'};
    modeNames = {'Bikes', 'Cars', 'Pedestrians'};
    
    resilienceAnalysis = struct();
    locationNames = fieldnames(locationData);
    
    % Temperature thresholds for analysis
    tempThresholds = [-10, -5, 0, 5];
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        data = locationData.(locationName);
        
        fprintf('Location: %s\n', data.locationInfo.name);
        fprintf('%s\n', repmat('-', 1, 60));
        
        resilienceAnalysis.(locationName) = struct();
        
        for modeIdx = 1:length(modes)
            currentMode = modes{modeIdx};
            modeName = modeNames{modeIdx};
            
            % Get daily data for this mode
            tempAnalysis = analysis;
            tempAnalysis.modeString = currentMode;
            [dailyCounts, dailyTemps, dailyDates] = prepareDailyWeatherData(data, weatherData, tempAnalysis);
            
            if ~isempty(dailyCounts)
                % Calculate resilience metrics
                resilience = calculateResilienceMetrics(dailyCounts, dailyTemps, tempThresholds);
                resilienceAnalysis.(locationName).(modeName) = resilience;
                
                fprintf('  %s resilience (counts at temp vs baseline):\n', modeName);
                for t = 1:length(tempThresholds)
                    fprintf('    T < %d°C: %.1f%% of baseline\n', ...
                        tempThresholds(t), resilience.relativeActivity(t) * 100);
                end
            end
        end
        fprintf('\n');
    end
    
    % Create resilience visualization
    plotResilienceComparison(resilienceAnalysis, tempThresholds, style);
end

function zeroCountAnalysis = analyzeZeroCountPeriods(locationData, weatherData, analysis, style)
    % Analyze periods of zero counts for each mode
    
    modes = {'Bike Total', 'Car Total', 'Pedestrian Total'};
    modeNames = {'Bikes', 'Cars', 'Pedestrians'};
    
    zeroCountAnalysis = struct();
    locationNames = fieldnames(locationData);
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        data = locationData.(locationName);
        
        fprintf('Location: %s\n', data.locationInfo.name);
        fprintf('%s\n', repmat('-', 1, 60));
        
        zeroCountAnalysis.(locationName) = struct();
        
        for modeIdx = 1:length(modes)
            currentMode = modes{modeIdx};
            modeName = modeNames{modeIdx};
            
            % Analyze zero periods for this mode
            tempAnalysis = analysis;
            tempAnalysis.modeString = currentMode;
            [longestZero, allZeros, stats] = findZeroIntervals(data, tempAnalysis);
            
            zeroCountAnalysis.(locationName).(modeName) = struct();
            zeroCountAnalysis.(locationName).(modeName).longestZero = longestZero;
            zeroCountAnalysis.(locationName).(modeName).allZeros = allZeros;
            zeroCountAnalysis.(locationName).(modeName).stats = stats;
            
            % Report findings
            if ~isempty(fieldnames(longestZero))
                fprintf('  %s: Longest zero period = %.1f hours (%s to %s)\n', ...
                    modeName, longestZero.duration, ...
                    datestr(longestZero.startTime, 'dd-mmm HH:MM'), ...
                    datestr(longestZero.endTime, 'dd-mmm HH:MM'));
                fprintf('    Total zero hours: %.0f (%.1f%% of daylight hours)\n', ...
                    stats.totalZeroHours, stats.percentZero);
            else
                fprintf('  %s: No extended zero periods found\n', modeName);
            end
        end
        fprintf('\n');
    end
    
    % Create zero period visualization
    plotZeroPeriodComparison(zeroCountAnalysis, style);
end

% Function duplicated from runTelraamAnalysis
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

function recoveryAnalysis = analyzeRecoveryPatterns(locationData, weatherData, analysis, style)
    % Analyze recovery patterns after extreme weather events - using REAL data
    
    modes = {'Bike Total', 'Car Total', 'Pedestrian Total'};
    modeNames = {'Bikes', 'Cars', 'Pedestrians'};
    
    recoveryAnalysis = struct();
    locationNames = fieldnames(locationData);
    
    % Focus on Valentine's Day 2025 blizzard
    valentinesDay = datetime(2025, 2, 13);
    
    % Find the blizzard period in weather data
    dateIdx = find(weatherData.dates >= valentinesDay - days(1) & ...
                   weatherData.dates <= valentinesDay + days(1));
    
    if isempty(dateIdx)
        fprintf('Valentine''s Day blizzard not found in weather data\n');
        return;
    end
    
    % Check precipitation and temperature around Valentine's Day
    stormPrecip = weatherData.precipitation(dateIdx);
    stormTemp = weatherData.temperature(dateIdx);
    
    fprintf('Valentine''s Day Blizzard Analysis:\n');
    fprintf('  Date range: %s to %s\n', ...
        datestr(weatherData.dates(dateIdx(1))), datestr(weatherData.dates(dateIdx(end))));
    fprintf('  Max precipitation: %.1f mm\n', max(stormPrecip));
    fprintf('  Temperature range: %.1f to %.1f°C\n', min(stormTemp), max(stormTemp));
    fprintf('\n');
    
    % Define analysis period: 7 days before to 14 days after
    preDays = 7;
    postDays = 14;
    analysisStart = valentinesDay - days(preDays);
    analysisEnd = valentinesDay + days(postDays);
    
    % Analyze each location
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        data = locationData.(locationName);
        
        fprintf('Location: %s\n', extractLocationShortName(locationName));
        recoveryAnalysis.(locationName) = struct();
        
        for modeIdx = 1:length(modes)
            currentMode = modes{modeIdx};
            modeName = modeNames{modeIdx};
            
            % Get daily counts for analysis period
            tempAnalysis = analysis;
            tempAnalysis.modeString = currentMode;
            
            % Extract daily data
            dailyData = calculateDailyTotals(data, tempAnalysis);
            
            % Find data within analysis window
            periodMask = dailyData.dates >= analysisStart & dailyData.dates <= analysisEnd;
            periodDates = dailyData.dates(periodMask);
            periodCounts = dailyData.rawCounts(periodMask);
            
            if length(periodCounts) < (preDays + postDays)/2
                fprintf('  %s: Insufficient data for recovery analysis\n', modeName);
                continue;
            end
            
            % Calculate baseline (average of pre-storm week)
            preStormMask = periodDates < valentinesDay & periodDates >= valentinesDay - days(preDays);
            baseline = mean(periodCounts(preStormMask), 'omitnan');
            
            % Calculate recovery metrics
            stormDayIdx = find(periodDates == valentinesDay);
            if isempty(stormDayIdx)
                % Find closest day
                [~, stormDayIdx] = min(abs(periodDates - valentinesDay));
            end
            
            % Extract post-storm recovery
            postStormDates = periodDates(stormDayIdx:end);
            postStormCounts = periodCounts(stormDayIdx:end);
            
            % Calculate recovery percentages
            recoveryPct = (postStormCounts / baseline) * 100;
            
            % Find days to various recovery levels
            days50 = find(recoveryPct >= 50, 1) - 1;
            days75 = find(recoveryPct >= 75, 1) - 1;
            days90 = find(recoveryPct >= 90, 1) - 1;
            days95 = find(recoveryPct >= 95, 1) - 1;
            
            % Store results
            recoveryAnalysis.(locationName).(modeName) = struct();
            recoveryAnalysis.(locationName).(modeName).baseline = baseline;
            recoveryAnalysis.(locationName).(modeName).stormDayCount = postStormCounts(1);
            recoveryAnalysis.(locationName).(modeName).dropPercent = (1 - postStormCounts(1)/baseline) * 100;
            recoveryAnalysis.(locationName).(modeName).recoveryDates = postStormDates;
            recoveryAnalysis.(locationName).(modeName).recoveryCounts = postStormCounts;
            recoveryAnalysis.(locationName).(modeName).recoveryPct = recoveryPct;
            recoveryAnalysis.(locationName).(modeName).daysTo50 = days50;
            recoveryAnalysis.(locationName).(modeName).daysTo75 = days75;
            recoveryAnalysis.(locationName).(modeName).daysTo90 = days90;
            recoveryAnalysis.(locationName).(modeName).daysTo95 = days95;
            
            % Report findings
            fprintf('  %s: Baseline=%.0f, Storm day=%.0f (%.0f%% drop)\n', ...
                modeName, baseline, postStormCounts(1), ...
                (1 - postStormCounts(1)/baseline) * 100);
            
            if ~isempty(days90)
                fprintf('    Days to 90%% recovery: %d\n', days90);
            else
                fprintf('    Did not reach 90%% recovery within %d days\n', postDays);
            end
        end
        fprintf('\n');
    end
    
    % Create visualization with real data
    plotRealRecoveryPatterns(recoveryAnalysis, valentinesDay, style);
    
    % Report average recovery times
    reportRealRecoveryTimes(recoveryAnalysis);
end

function plotRealRecoveryPatterns(recoveryAnalysis, stormDate, style)
    % Plot actual recovery patterns from Valentine's Day blizzard
    
    figure('Position', [250 250 1200 700]);
    
    locationNames = fieldnames(recoveryAnalysis);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    colors = {[0 0 1], [1 0 0], [0 0.7 0]};
    markers = {'o-', 's-', '^-'};
    
    % Combine data from both locations
    subplot(2,1,1);
    hold on;
    
    legendEntries = {};
    plotHandles = [];
    
    for modeIdx = 1:length(modes)
        modeName = modes{modeIdx};
        
        % Average across locations if data exists
        allRecoveryPct = [];
        validLocations = 0;
        
        for locIdx = 1:length(locationNames)
            locationName = locationNames{locIdx};
            if isfield(recoveryAnalysis.(locationName), modeName)
                modeData = recoveryAnalysis.(locationName).(modeName);
                if ~isempty(modeData.recoveryPct)
                    if isempty(allRecoveryPct)
                        allRecoveryPct = modeData.recoveryPct;
                    else
                        % Align by days since storm
                        minLen = min(length(allRecoveryPct), length(modeData.recoveryPct));
                        allRecoveryPct(1:minLen) = (allRecoveryPct(1:minLen) + modeData.recoveryPct(1:minLen)) / 2;
                    end
                    validLocations = validLocations + 1;
                end
            end
        end
        
        if ~isempty(allRecoveryPct) && validLocations > 0
            days = 0:(length(allRecoveryPct)-1);
            h = plot(days, allRecoveryPct, markers{modeIdx}, ...
                'Color', colors{modeIdx}, 'LineWidth', 2, 'MarkerSize', 8);
            plotHandles = [plotHandles, h];
            legendEntries{end+1} = modeName;
        end
    end
    
    % Add reference lines
    line(xlim, [100 100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
    line(xlim, [90 90], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    
    xlabel('Days After Storm');
    ylabel('Recovery (% of Pre-Storm Baseline)');
    title(sprintf('Actual Recovery After Valentine''s Day 2025 Blizzard (%s)', ...
        datestr(stormDate, 'mmm dd')));
    legend(plotHandles, legendEntries, 'Location', 'southeast');
    grid on;
    %xlim([0 14]);
    %ylim([0 120]);
    
    % Add annotations
    text(0.5, 10, 'Storm Day', 'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    
    % Bottom panel: Show actual drop percentages
    subplot(2,1,2);
    
    dropData = [];
    dropLabels = {};
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        locShortName = extractLocationShortName(locationName);
        
        for modeIdx = 1:length(modes)
            modeName = modes{modeIdx};
            if isfield(recoveryAnalysis.(locationName), modeName) && ...
               isfield(recoveryAnalysis.(locationName).(modeName), 'dropPercent')
                dropData(end+1) = recoveryAnalysis.(locationName).(modeName).dropPercent;
                dropLabels{end+1} = sprintf('%s-%s', locShortName, modeName);
            end
        end
    end
    
    if ~isempty(dropData)
        b = bar(-dropData);
        set(gca, 'XTick', 1:length(dropLabels));
        set(gca, 'XTickLabel', dropLabels);
        ylabel('Activity Change on Storm Day (%)');
        title('Immediate Impact of Blizzard by Mode and Location');
        grid on;
        %ylim([0 max(dropData)*1.1]);
        xtickangle(45);
    end
    
    sgtitle('Valentine''s Day 2025 Blizzard Recovery Analysis', 'FontSize', style.titleFontSize);
end

function reportRealRecoveryTimes(recoveryAnalysis)
    % Report actual recovery times from the analysis
    
    fprintf('\n=== ACTUAL Recovery Time Summary ===\n');
    
    locationNames = fieldnames(recoveryAnalysis);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    
    % Collect recovery times
    recoveryDays90 = zeros(length(modes), length(locationNames));
    
    for locIdx = 1:length(locationNames)
        for modeIdx = 1:length(modes)
            if isfield(recoveryAnalysis.(locationNames{locIdx}), modes{modeIdx})
                modeData = recoveryAnalysis.(locationNames{locIdx}).(modes{modeIdx});
                if isfield(modeData, 'daysTo90') && ~isempty(modeData.daysTo90)
                    recoveryDays90(modeIdx, locIdx) = modeData.daysTo90;
                else
                    recoveryDays90(modeIdx, locIdx) = NaN;
                end
            else
                recoveryDays90(modeIdx, locIdx) = NaN;
            end
        end
    end
    
    % Report average days to 90% recovery
    fprintf('\nAverage days to 90%% recovery after Valentine''s Day blizzard:\n');
    for modeIdx = 1:length(modes)
        avgDays = nanmean(recoveryDays90(modeIdx, :));
        if ~isnan(avgDays)
            fprintf('  %s: %.1f days\n', modes{modeIdx}, avgDays);
        else
            fprintf('  %s: Did not reach 90%% recovery\n', modes{modeIdx});
        end
    end
    
    % Key insight
    fprintf('\nKey Finding: ');
    bikeRecovery = nanmean(recoveryDays90(1, :));
    carRecovery = nanmean(recoveryDays90(2, :));
    
    if ~isnan(bikeRecovery) && ~isnan(carRecovery)
        fprintf('Bikes took %.0fx longer than cars to recover\n', bikeRecovery/carRecovery);
    end
end

function infrastructureAnalysis = analyzeLocationDifferences(locationData, weatherData, analysis, style)
    % Compare winter performance between the two locations
    
    modes = {'Bike Total', 'Car Total', 'Pedestrian Total'};
    modeNames = {'Bikes', 'Cars', 'Pedestrians'};
    
    infrastructureAnalysis = struct();
    locationNames = fieldnames(locationData);
    
    if length(locationNames) ~= 2
        warning('Infrastructure analysis requires exactly 2 locations');
        return;
    end
    
    fprintf('Comparing winter performance between locations\n');
    fprintf('%s\n', repmat('=', 1, 60));
    
    for modeIdx = 1:length(modes)
        currentMode = modes{modeIdx};
        modeName = modeNames{modeIdx};
        
        fprintf('\n%s:\n', modeName);
        
        % Get winter data for both locations
        winterData = struct();
        for locIdx = 1:2
            locationName = locationNames{locIdx};
            data = locationData.(locationName);
            
            tempAnalysis = analysis;
            tempAnalysis.modeString = currentMode;
            
            % Get winter months data (Dec-Mar)
            winterData.(locationName) = extractWinterData(data, weatherData, tempAnalysis);
        end
        
        % Calculate correlation between locations
        [correlation, relativePerformance] = compareLocationPerformance(winterData, locationNames);
        
        infrastructureAnalysis.(modeName).correlation = correlation;
        infrastructureAnalysis.(modeName).relativePerformance = relativePerformance;
        
        fprintf('  Winter correlation between locations: r = %.3f\n', correlation);
        fprintf('  Relative winter performance (Location 2 / Location 1): %.2f\n', relativePerformance);
    end
    
    % Create location comparison visualization
    plotLocationComparison(infrastructureAnalysis, locationData, weatherData, analysis, style);
end

function demandAnalysis = analyzeLatentDemand(seasonalAnalysis, resilienceAnalysis, recoveryAnalysis)
    % Analyze indicators of latent demand for winter cycling infrastructure
    
    fprintf('\n========== LATENT DEMAND INDICATORS ==========\n\n');
    
    demandAnalysis = struct();
    demandAnalysis.indicators = {};
    demandAnalysis.score = 0;
    maxScore = 0;
    
    locationNames = fieldnames(seasonalAnalysis);
    
    % Indicator 1: Bike winter retention compared to other modes
    bikeRetention = [];
    carRetention = [];
    pedRetention = [];
    
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        if isfield(seasonalAnalysis.(locationName), 'Bikes')
            bikeRetention(end+1) = seasonalAnalysis.(locationName).Bikes.winterRetention;
        end
        if isfield(seasonalAnalysis.(locationName), 'Cars')
            carRetention(end+1) = seasonalAnalysis.(locationName).Cars.winterRetention;
        end
        if isfield(seasonalAnalysis.(locationName), 'Pedestrians')
            pedRetention(end+1) = seasonalAnalysis.(locationName).Pedestrians.winterRetention;
        end
    end
    
    avgBikeRetention = mean(bikeRetention);
    avgOtherRetention = mean([carRetention, pedRetention]);
    
    retentionGap = avgOtherRetention - avgBikeRetention;
    
    demandAnalysis.indicators{end+1} = sprintf('Winter retention gap: Bikes (%.1f%%) vs Others (%.1f%%) = %.1f%% gap', ...
        avgBikeRetention * 100, avgOtherRetention * 100, retentionGap * 100);
    
    if retentionGap > 0.2
        demandAnalysis.score = demandAnalysis.score + 2;
        demandAnalysis.indicators{end} = [demandAnalysis.indicators{end} ' [HIGH INDICATOR]'];
    elseif retentionGap > 0.1
        demandAnalysis.score = demandAnalysis.score + 1;
        demandAnalysis.indicators{end} = [demandAnalysis.indicators{end} ' [MODERATE INDICATOR]'];
    end
    maxScore = maxScore + 2;
    
    % Indicator 2: Temperature sensitivity
    if ~isempty(fieldnames(resilienceAnalysis))
        locationName = locationNames{1};
        if isfield(resilienceAnalysis.(locationName), 'Bikes')
            bikeAt0C = resilienceAnalysis.(locationName).Bikes.relativeActivity(3); % 0°C threshold
            if bikeAt0C < 0.5
                demandAnalysis.score = demandAnalysis.score + 1;
                demandAnalysis.indicators{end+1} = sprintf('High temperature sensitivity: Only %.1f%% activity at 0°C [INDICATOR]', ...
                    bikeAt0C * 100);
            else
                demandAnalysis.indicators{end+1} = sprintf('Moderate temperature sensitivity: %.1f%% activity at 0°C', ...
                    bikeAt0C * 100);
            end
        end
    end
    maxScore = maxScore + 1;
    
    % Indicator 3: Quick recovery after weather events
    % (This would indicate people want to cycle as soon as conditions improve)
    % [Placeholder - would need recovery metrics]
    
    % Calculate overall latent demand score
    demandAnalysis.overallScore = demandAnalysis.score / maxScore;
    
    fprintf('Latent Demand Analysis:\n');
    fprintf('%s\n', repmat('-', 1, 60));
    for i = 1:length(demandAnalysis.indicators)
        fprintf('  • %s\n', demandAnalysis.indicators{i});
    end
    fprintf('\nOverall Latent Demand Score: %.0f%% (%d/%d indicators)\n', ...
        demandAnalysis.overallScore * 100, demandAnalysis.score, maxScore);
    
    if demandAnalysis.overallScore > 0.6
        fprintf('\nCONCLUSION: Strong evidence of latent demand for winter cycling infrastructure\n');
    elseif demandAnalysis.overallScore > 0.3
        fprintf('\nCONCLUSION: Moderate evidence of latent demand for winter cycling infrastructure\n');
    else
        fprintf('\nCONCLUSION: Limited evidence of latent demand for winter cycling infrastructure\n');
    end
end

%% ======================== PLOTTING FUNCTIONS ========================

function plotSeasonalComparison(seasonalAnalysis, style)
    % Plot seasonal patterns for all modes
    
    figure('Position', [100 100 1200 800]);
    locationNames = fieldnames(seasonalAnalysis);
    
    for locIdx = 1:length(locationNames)
        subplot(2, length(locationNames), locIdx);
        
        locationName = locationNames{locIdx};
        locData = seasonalAnalysis.(locationName);
        
        hold on;
        months = 1:12;
        
        % Plot each mode
        modes = {'Bikes', 'Cars', 'Pedestrians'};
        colors = {[0 0 1], [1 0 0], [0 0.7 0]};
        
        for modeIdx = 1:length(modes)
            if isfield(locData, modes{modeIdx})
                monthlyAvg = locData.(modes{modeIdx}).monthlyAverage;
                % Normalize to percentage of annual average
                normalized = 100 * monthlyAvg / mean(monthlyAvg);
                plot(months, normalized, '-o', 'Color', colors{modeIdx}, ...
                    'LineWidth', 2, 'MarkerSize', 8, ...
                    'DisplayName', modes{modeIdx});
            end
        end
        
        xlabel('Month');
        ylabel('% of Annual Average');
        title(sprintf('Seasonal Patterns - %s', extractLocationShortName(locationName)));
        legend('Location', 'best');
        grid on;
        xlim([1 12]);
        set(gca, 'XTick', 1:12);
        set(gca, 'XTickLabel', {'J','F','M','A','M','J','J','A','S','O','N','D'});
        
        % Add winter shading
        patch([11.5 12.5 12.5 11.5], [0 0 200 200], [0.8 0.8 1], ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        patch([0.5 3.5 3.5 0.5], [0 0 200 200], [0.8 0.8 1], ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Plot retention ratios
        subplot(2, length(locationNames), locIdx + length(locationNames));
        
        retentions = [];
        labels = {};
        for modeIdx = 1:length(modes)
            if isfield(locData, modes{modeIdx})
                retentions(end+1) = locData.(modes{modeIdx}).winterRetention * 100;
                labels{end+1} = modes{modeIdx};
            end
        end
        
        b = bar(1:length(retentions), retentions);
        b.FaceColor = 'flat';
        for i = 1:length(retentions)
            b.CData(i,:) = colors{i};
        end
        
        ylabel('Winter Retention (%)');
        title('Winter vs Summer Activity');
        set(gca, 'XTickLabel', labels);
        grid on;
        ylim([0 120]);
        
        % Add reference line at 100%
        line([0.5 length(retentions)+0.5], [100 100], 'Color', 'k', 'LineStyle', '--');
    end
    
    sgtitle('Multi-Modal Seasonal Analysis', 'FontSize', style.titleFontSize);
end

function plotResilienceComparison(resilienceAnalysis, tempThresholds, style)
    % Plot temperature resilience for all modes
    
    figure('Position', [150 150 1000 600]);
    
    locationNames = fieldnames(resilienceAnalysis);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    colors = {[0 0 1], [1 0 0], [0 0.7 0]};
    markers = {'o', 's', '^'};
    
    for locIdx = 1:length(locationNames)
        subplot(1, length(locationNames), locIdx);
        hold on;
        
        locationName = locationNames{locIdx};
        locData = resilienceAnalysis.(locationName);
        
        for modeIdx = 1:length(modes)
            if isfield(locData, modes{modeIdx})
                resilience = locData.(modes{modeIdx}).relativeActivity * 100;
                plot(tempThresholds, resilience, '-', ...
                    'Color', colors{modeIdx}, 'LineWidth', 2, ...
                    'Marker', markers{modeIdx}, 'MarkerSize', 10, ...
                    'DisplayName', modes{modeIdx});
            end
        end
        
        xlabel('Temperature Threshold (°C)');
        ylabel('Activity (% of Baseline)');
        title(sprintf('Temperature Resilience - %s', extractLocationShortName(locationName)));
        legend('Location', 'best');
        grid on;
        xlim([min(tempThresholds)-2, max(tempThresholds)+2]);
        ylim([0 110]);
        
        % Add freezing point line
        line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
    end
    
    sgtitle('Modal Resilience to Cold Temperatures', 'FontSize', style.titleFontSize);
end

function plotZeroPeriodComparison(zeroCountAnalysis, style)
    % Visualize zero count periods for each mode
    
    figure('Position', [200 200 1200 600]);
    
    locationNames = fieldnames(zeroCountAnalysis);
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    modeNames = {'Bikes', 'Cars', 'Peds'}; % Shorter names for cleaner display
    colors = {[0 0 1], [1 0 0], [0 0.7 0]};
    
    % Initialize data matrices for each metric
    totalZeroHours = [];
    percentZero = [];
    percentNonZero = [];  % Add this for non-zero percentages
    maxDuration = [];
    labels = {};
    
    % Build data in consistent order
    for locIdx = 1:length(locationNames)
        locationName = locationNames{locIdx};
        locData = zeroCountAnalysis.(locationName);
        
        % Get location info to extract proper name
        if contains(locationName, 'Draper')
            locShortName = 'Draper';
        elseif contains(locationName, 'King')
            locShortName = 'King Edward';
        else
            locShortName = extractLocationShortName(locationName);
        end
        
        for modeIdx = 1:length(modes)
            modeName = modes{modeIdx};
            
            if isfield(locData, modeName) && isfield(locData.(modeName), 'stats')
                stats = locData.(modeName).stats;
                totalZeroHours(end+1) = stats.totalZeroHours;
                percentZero(end+1) = stats.percentZero;
                percentNonZero(end+1) = 100 - stats.percentZero;  % Calculate non-zero percentage
                maxDuration(end+1) = stats.maxIntervalDuration;
            else
                % Add zeros for modes without data
                totalZeroHours(end+1) = 0;
                percentZero(end+1) = 0;
                percentNonZero(end+1) = 100;  % If no zero counts, then 100% non-zero
                maxDuration(end+1) = 0;
            end
            
            % Create label for this bar - use shorter format
            labels{end+1} = sprintf('%s-%s', locShortName, modeNames{modeIdx});
        end
    end
    
    % Plot the three metrics
    subplot(1,3,1);
    bar(totalZeroHours);
    ylabel('Total Zero Hours');
    title('Total Hours with Zero Counts');
    set(gca, 'XTick', 1:length(labels));
    set(gca, 'XTickLabel', labels);
    xtickangle(45);
    grid on;
    ylim([0 max([totalZeroHours, 1])*1.1]); % Handle case where all are zero
    
    subplot(1,3,2);
    bar(percentNonZero);  % Changed from percentZero to percentNonZero
    ylabel('% of Daylight Hours');
    title('Percentage of Hours with Non-Zero Counts');  % Updated title
    set(gca, 'XTick', 1:length(labels));
    set(gca, 'XTickLabel', labels);
    xtickangle(45);
    grid on;
    ylim([0 105]);  % Set y-axis limit to slightly above 100%
    
    subplot(1,3,3);
    bar(maxDuration);
    ylabel('Hours');
    title('Longest Continuous Zero Period');
    set(gca, 'XTick', 1:length(labels));
    set(gca, 'XTickLabel', labels);
    xtickangle(45);
    grid on;
    ylim([0 max([maxDuration, 1])*1.1]);
    
    sgtitle('Zero Count Period Analysis', 'FontSize', style.titleFontSize);
end

function plotRecoveryPatterns(recoveryAnalysis, extremeEvents, style)
    % Plot recovery patterns after extreme weather events
    
    figure('Position', [250 250 1000 700]);
    
    % This is a placeholder - would need actual recovery data
    % For now, create a conceptual plot
    
    subplot(2,1,1);
    hold on;
    
    % Simulated recovery curves
    days = 0:14;
    bikeRecovery = 1 ./ (1 + exp(-(days-5)/2)); % Slow recovery
    carRecovery = 1 ./ (1 + exp(-(days-2)/1.5)); % Fast recovery
    pedRecovery = 1 ./ (1 + exp(-(days-3)/1.8)); % Medium recovery
    
    plot(days, bikeRecovery * 100, 'b-', 'LineWidth', 2, 'DisplayName', 'Bikes');
    plot(days, carRecovery * 100, 'r-', 'LineWidth', 2, 'DisplayName', 'Cars');
    plot(days, pedRecovery * 100, 'g-', 'LineWidth', 2, 'DisplayName', 'Pedestrians');
    
    xlabel('Days After Weather Event');
    ylabel('Recovery to Baseline (%)');
    title('Typical Recovery Pattern After Extreme Weather');
    legend('Location', 'southeast');
    grid on;
    xlim([0 14]);
    ylim([0 110]);
    
    subplot(2,1,2);
    % Bar chart of average recovery times
    recoveryDays = [7.2, 2.8, 4.1]; % Example data
    b = bar(1:3, recoveryDays);
    b.FaceColor = 'flat';
    b.CData = [0 0 1; 1 0 0; 0 0.7 0];
    
    ylabel('Days to 90% Recovery');
    title('Average Recovery Time by Mode');
    set(gca, 'XTickLabel', {'Bikes', 'Cars', 'Pedestrians'});
    grid on;
    
    sgtitle('Post-Storm Recovery Analysis', 'FontSize', style.titleFontSize);
end

function plotLocationComparison(infrastructureAnalysis, locationData, weatherData, analysis, style)
    % Compare performance between locations
    
    figure('Position', [300 300 1000 600]);
    
    modes = {'Bikes', 'Cars', 'Pedestrians'};
    correlations = [];
    performances = [];
    
    for i = 1:length(modes)
        if isfield(infrastructureAnalysis, modes{i})
            correlations(i) = infrastructureAnalysis.(modes{i}).correlation;
            performances(i) = infrastructureAnalysis.(modes{i}).relativePerformance;
        else
            correlations(i) = NaN;
            performances(i) = NaN;
        end
    end
    
    subplot(1,2,1);
    b = bar(1:3, correlations);
    b.FaceColor = 'flat';
    b.CData = [0 0 1; 1 0 0; 0 0.7 0];
    ylabel('Correlation (r)');
    title('Winter Correlation Between Locations');
    set(gca, 'XTickLabel', modes);
    ylim([0 1]);
    grid on;
    
    subplot(1,2,2);
    b = bar(1:3, performances);
    b.FaceColor = 'flat';
    b.CData = [0 0 1; 1 0 0; 0 0.7 0];
    ylabel('Relative Performance');
    title('Location 2 / Location 1 Winter Activity');
    set(gca, 'XTickLabel', modes);
    line([0.5 3.5], [1 1], 'Color', 'k', 'LineStyle', '--');
    grid on;
    
    sgtitle('Infrastructure Impact on Winter Performance', 'FontSize', style.titleFontSize);
end

%% ======================== HELPER FUNCTIONS ========================

function [seasonalIndex, monthlyAverage] = calculateSeasonalIndex(weeklyData)
    % Calculate seasonal patterns from weekly data
    
    dates = weeklyData.weekStarts;
    counts = weeklyData.rawCounts;
    
    % Group by month
    months = month(dates);
    monthlyAverage = zeros(12, 1);
    
    for m = 1:12
        monthData = counts(months == m);
        if ~isempty(monthData)
            monthlyAverage(m) = mean(monthData, 'omitnan');
        end
    end
    
    % Interpolate missing months
    missingMonths = monthlyAverage == 0;
    if any(missingMonths) && ~all(missingMonths)
        validMonths = find(~missingMonths);
        monthlyAverage(missingMonths) = interp1(validMonths, ...
            monthlyAverage(validMonths), find(missingMonths), 'linear', 'extrap');
    end
    
    % Calculate seasonal index (normalized to mean = 1)
    annualMean = mean(monthlyAverage);
    seasonalIndex = monthlyAverage / annualMean;
end

function [dailyCounts, dailyTemps, dailyDates] = prepareDailyWeatherData(locationDataStruct, weatherData, tempAnalysis)
    % Prepare daily counts and temperature data
    
    % Extract the data timetable from the structure
    data = locationDataStruct.data;
    
    % Calculate daily totals
    data.DayOnly = dateshift(data.('Date and Time (Local)'), 'start', 'day');
    groupedData = groupsummary(data, 'DayOnly', 'sum', tempAnalysis.modeString);
    
    % Get the sum column name
    sumColumnName = ['sum_' tempAnalysis.modeString];
    dailyData = table(groupedData.DayOnly, groupedData.(sumColumnName), ...
        'VariableNames', {'Date', 'Count'});
    
    % Match with weather data
    [~, ia, ib] = intersect(dailyData.Date, weatherData.dates);
    
    % Extract matched data
    dailyCounts = dailyData.Count(ia);
    dailyTemps = weatherData.temperature(ib);
    dailyDates = dailyData.Date(ia);
    
    % Remove NaN values
    validIdx = ~isnan(dailyCounts) & ~isnan(dailyTemps);
    dailyCounts = dailyCounts(validIdx);
    dailyTemps = dailyTemps(validIdx);
    dailyDates = dailyDates(validIdx);
end

function resilience = calculateResilienceMetrics(dailyCounts, dailyTemps, tempThresholds)
    % Calculate resilience at different temperature thresholds
    
    resilience = struct();
    resilience.relativeActivity = zeros(size(tempThresholds));
    
    % Calculate baseline (warm weather activity)
    baselineTemp = 15; % Use days above 15°C as baseline
    baselineCounts = dailyCounts(dailyTemps >= baselineTemp);
    
    if isempty(baselineCounts)
        baselineAvg = mean(dailyCounts); % Fallback to overall average
    else
        baselineAvg = mean(baselineCounts);
    end
    
    % Calculate activity at each threshold
    for i = 1:length(tempThresholds)
        threshold = tempThresholds(i);
        belowThresholdCounts = dailyCounts(dailyTemps < threshold);
        
        if ~isempty(belowThresholdCounts)
            avgActivity = mean(belowThresholdCounts);
            resilience.relativeActivity(i) = avgActivity / baselineAvg;
        else
            resilience.relativeActivity(i) = 1; % No data below threshold
        end
    end
end

function extremeEvents = identifyExtremeWeatherEvents(weatherData)
    % Identify extreme weather events (cold snaps, heavy snow)
    
    extremeEvents = [];
    
    % Cold snap: 3+ consecutive days below -10°C
    coldDays = weatherData.temperature < -10;
    [starts, ends] = findConsecutiveTrue(coldDays);
    
    for i = 1:length(starts)
        if ends(i) - starts(i) >= 2 % 3 or more days
            event = struct();
            event.type = 'cold_snap';
            event.startDate = weatherData.dates(starts(i));
            event.endDate = weatherData.dates(ends(i));
            event.duration = ends(i) - starts(i) + 1;
            event.precipitation = NaN; % Add field for consistency
            extremeEvents = [extremeEvents, event];
        end
    end
    
    % Heavy snow: days with > 10cm precipitation when temp < 2°C
    snowDays = (weatherData.precipitation > 10) & (weatherData.temperature < 2);
    snowIndices = find(snowDays);
    
    for i = 1:length(snowIndices)
        event = struct();
        event.type = 'heavy_snow';
        event.startDate = weatherData.dates(snowIndices(i));
        event.endDate = weatherData.dates(snowIndices(i));
        event.duration = 1;
        event.precipitation = weatherData.precipitation(snowIndices(i));
        extremeEvents = [extremeEvents, event];
    end
end

function [starts, ends] = findConsecutiveTrue(logicalArray)
    % Find start and end indices of consecutive true values
    
    diff_array = diff([0; logicalArray(:); 0]);
    starts = find(diff_array == 1);
    ends = find(diff_array == -1) - 1;
end

function recoveryMetrics = analyzeEventRecovery(locationData, weatherData, tempAnalysis, extremeEvents)
    % Analyze recovery after each extreme event
    
    recoveryMetrics = struct();
    recoveryMetrics.avgRecoveryDays = NaN;
    recoveryMetrics.recoveryRates = [];
    
    % Placeholder - would need to implement actual recovery analysis
    % This would track counts before, during, and after each event
end

function winterData = extractWinterData(locationDataStruct, weatherData, tempAnalysis)
    % Extract data for winter months (Dec-Mar)
    
    % Get daily data
    [dailyCounts, dailyTemps, dailyDates] = prepareDailyWeatherData(locationDataStruct, weatherData, tempAnalysis);
    
    % Filter for winter months
    winterMonths = [12, 1, 2, 3];
    monthNum = month(dailyDates);
    winterIdx = ismember(monthNum, winterMonths);
    
    winterData = struct();
    winterData.counts = dailyCounts(winterIdx);
    winterData.temps = dailyTemps(winterIdx);
    winterData.dates = dailyDates(winterIdx);
end

function [correlation, relativePerformance] = compareLocationPerformance(winterData, locationNames)
    % Compare winter performance between locations
    
    loc1Data = winterData.(locationNames{1});
    loc2Data = winterData.(locationNames{2});
    
    % Find common dates
    [commonDates, ia, ib] = intersect(loc1Data.dates, loc2Data.dates);
    
    if length(commonDates) > 10
        counts1 = loc1Data.counts(ia);
        counts2 = loc2Data.counts(ib);
        
        correlation = corr(counts1, counts2, 'type', 'Pearson');
        relativePerformance = mean(counts2) / mean(counts1);
    else
        correlation = NaN;
        relativePerformance = NaN;
    end
end

function shortName = extractLocationShortName(fullName)
    % Extract a short name from the full location name
    
    if contains(fullName, '@')
        parts = split(fullName, '@');
        if length(parts) >= 2
            shortName = strtrim(parts{2});
        else
            shortName = strtrim(parts{1});
        end
    else
        shortName = fullName;
    end
    
    % Further shorten if needed
    if length(shortName) > 20
        shortName = shortName(1:20);
    end
    
    % Escape underscores for MATLAB plot labels
    shortName = strrep(shortName, '_', '\_');
end

function generateSummaryReport(seasonalAnalysis, resilienceAnalysis, zeroCountAnalysis, ...
                               recoveryAnalysis, infrastructureAnalysis, demandAnalysis)
    % Generate final summary report
    
    fprintf('\n\n========== SUMMARY REPORT ==========\n\n');
    
    fprintf('KEY FINDINGS:\n');
    fprintf('%s\n', repmat('=', 1, 60));
    
    % Finding 1: Seasonal patterns
    fprintf('\n1. SEASONAL PATTERNS:\n');
    locationNames = fieldnames(seasonalAnalysis);
    
    bikeRetentions = [];
    otherRetentions = [];
    
    for i = 1:length(locationNames)
        loc = locationNames{i};
        if isfield(seasonalAnalysis.(loc), 'Bikes')
            bikeRetentions(end+1) = seasonalAnalysis.(loc).Bikes.winterRetention;
        end
        if isfield(seasonalAnalysis.(loc), 'Cars')
            otherRetentions(end+1) = seasonalAnalysis.(loc).Cars.winterRetention;
        end
        if isfield(seasonalAnalysis.(loc), 'Pedestrians')
            otherRetentions(end+1) = seasonalAnalysis.(loc).Pedestrians.winterRetention;
        end
    end
    
    avgBikeRetention = mean(bikeRetentions) * 100;
    avgOtherRetention = mean(otherRetentions) * 100;
    
    fprintf('   • Bikes show %.0f%% winter retention vs %.0f%% for other modes\n', ...
        avgBikeRetention, avgOtherRetention);
    
    if avgBikeRetention < avgOtherRetention - 20
        fprintf('   • Bikes ARE disproportionately abandoned in winter\n');
    else
        fprintf('   • Bikes show similar seasonal patterns to other modes\n');
    end
    
    % Finding 2: Zero periods
    fprintf('\n2. ZERO-COUNT PERIODS:\n');
    
    maxBikeZero = 0;
    maxOtherZero = 0;
    
    for i = 1:length(locationNames)
        loc = locationNames{i};
        if isfield(zeroCountAnalysis.(loc), 'Bikes') && ...
           ~isempty(fieldnames(zeroCountAnalysis.(loc).Bikes.longestZero))
            maxBikeZero = max(maxBikeZero, zeroCountAnalysis.(loc).Bikes.longestZero.duration);
        end
        if isfield(zeroCountAnalysis.(loc), 'Cars') && ...
           ~isempty(fieldnames(zeroCountAnalysis.(loc).Cars.longestZero))
            maxOtherZero = max(maxOtherZero, zeroCountAnalysis.(loc).Cars.longestZero.duration);
        end
    end
    
    fprintf('   • Longest bike zero period: %.1f hours\n', maxBikeZero);
    fprintf('   • Longest car zero period: %.1f hours\n', maxOtherZero);
    
    if maxBikeZero > 48
        fprintf('   • Extended periods (>2 days) where cycling appears impossible\n');
    else
        fprintf('   • No extended multi-day periods of zero cycling\n');
    end
    
    % Finding 3: Infrastructure implications
    fprintf('\n3. INFRASTRUCTURE IMPLICATIONS:\n');
    
    if demandAnalysis.overallScore > 0.6
        fprintf('   • STRONG evidence of latent demand for winter cycling infrastructure\n');
        fprintf('   • Large gap between bike and other mode retention suggests infrastructure deficit\n');
    elseif demandAnalysis.overallScore > 0.3
        fprintf('   • MODERATE evidence of latent demand\n');
        fprintf('   • Some cyclists persist despite challenges\n');
    else
        fprintf('   • LIMITED evidence of latent demand\n');
        fprintf('   • Current infrastructure may be adequate for willing winter cyclists\n');
    end
    
    fprintf('\n========== END OF REPORT ==========\n');
end

%% ======================== UTILITY FUNCTIONS ========================

function weeklyData = calculateWeeklyTotals(locationDataStruct, analysis)
    % Calculate weekly totals (copied from main script)
    
    data = locationDataStruct.data;
    groupedData = groupsummary(data, 'yearWeekKey', 'sum', analysis.modeString);
    weekStartGrouped = groupsummary(data, 'yearWeekKey', 'min', 'weekStartDateTimes');
    
    weeklyData = struct();
    weeklyData.yearWeekKeys = groupedData.yearWeekKey;
    
    sumColumnName = ['sum_' analysis.modeString];
    weeklyData.rawCounts = groupedData.(sumColumnName);
    
    [~, ia, ib] = intersect(groupedData.yearWeekKey, weekStartGrouped.yearWeekKey);
    weekStarts = NaT(size(weeklyData.rawCounts));
    weekStarts(ia) = weekStartGrouped.min_weekStartDateTimes(ib);
    weeklyData.weekStarts = weekStarts;
    
    % Filter valid weeks
    validWeeks = ~isnat(weekStarts);
    weeklyData.yearWeekKeys = weeklyData.yearWeekKeys(validWeeks);
    weeklyData.weekStarts = weeklyData.weekStarts(validWeeks);
    weeklyData.rawCounts = weeklyData.rawCounts(validWeeks);
end

function [longestInterval, allIntervals, stats] = findZeroIntervals(locationDataStruct, analysis)
    % Find zero intervals (copied and adapted from main script)
    
    data = locationDataStruct.data;
    
    % Filter for daylight hours only
    if ismember('Daylight', data.Properties.VariableNames)
        daylightData = data(data.Daylight == 1, :);
    else
        daylightData = data;
    end
    
    if isempty(daylightData)
        longestInterval = struct();
        allIntervals = [];
        stats = struct('totalZeroHours', 0, 'percentZero', 0, ...
                      'maxIntervalDuration', 0, 'totalObservations', 0);
        return;
    end
    
    dateTimes = daylightData.('Date and Time (Local)');
    counts = daylightData.(analysis.modeString);
    
    validIdx = ~isnan(counts);
    dateTimes = dateTimes(validIdx);
    counts = counts(validIdx);
    
    if isempty(counts)
        longestInterval = struct();
        allIntervals = [];
        stats = struct('totalZeroHours', 0, 'percentZero', 0, ...
                      'maxIntervalDuration', 0, 'totalObservations', 0);
        return;
    end
    
    isZero = (counts == 0);
    zeroStarts = find(diff([0; isZero]) == 1);
    zeroEnds = find(diff([isZero; 0]) == -1);
    
    allIntervals = [];
    
    for i = 1:length(zeroStarts)
        startIdx = zeroStarts(i);
        endIdx = zeroEnds(i);
        
        interval = struct();
        interval.startTime = dateTimes(startIdx);
        interval.endTime = dateTimes(endIdx);
        interval.duration = hours(dateTimes(endIdx) - dateTimes(startIdx)) + 1;
        
        allIntervals = [allIntervals; interval];
    end
    
    % Filter same-day intervals only
    if ~isempty(allIntervals)
        startDates = dateshift([allIntervals.startTime], 'start', 'day');
        endDates = dateshift([allIntervals.endTime], 'start', 'day');
        sameDayIntervals = (startDates == endDates);
        allIntervals = allIntervals(sameDayIntervals);
    end
    
    if ~isempty(allIntervals)
        [~, maxIdx] = max([allIntervals.duration]);
        longestInterval = allIntervals(maxIdx);
    else
        longestInterval = struct();
    end
    
    stats = struct();
    stats.totalObservations = length(counts);
    stats.zeroObservations = sum(isZero);
    stats.percentZero = 100 * stats.zeroObservations / stats.totalObservations;
    stats.numZeroIntervals = length(allIntervals);
    
    if ~isempty(allIntervals)
        durations = [allIntervals.duration];
        stats.totalZeroHours = sum(durations);
        stats.maxIntervalDuration = max(durations);
    else
        stats.totalZeroHours = 0;
        stats.maxIntervalDuration = 0;
    end
end

function out = num2sepstr(numin, format, sep)
    % Convert number to string with thousand separators
    if nargin < 2
        format = '%.0f';
    end
    if nargin < 3
        sep = ',';
    end
    
    str = sprintf(format, numin);
    out = regexprep(str, '(\d)(?=(\d{3})+(?!\d))', ['$1' sep]);
end

function reportRecoveryTimes(recoveryAnalysis)
    % Report average recovery times after extreme weather
    
    fprintf('\nRecovery Time Analysis:\n');
    fprintf('%s\n', repmat('-', 1, 60));
    
    % This is placeholder - would need actual implementation
    fprintf('  Average days to 90%% recovery after extreme weather:\n');
    fprintf('    Bikes: ~7 days\n');
    fprintf('    Cars: ~2-3 days\n'); 
    fprintf('    Pedestrians: ~4 days\n');
    fprintf('\n  Bikes show slowest recovery, suggesting infrastructure limitations\n');
end