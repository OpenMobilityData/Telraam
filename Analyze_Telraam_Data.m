
% Script to load and analyze data from a Telraam S2 traffic counter.
% The Telraam S2 is a computer vision based counter that monitors both sides
% of a roadway and uses object recognitation to count passages by different 
% vehicle types.
%
% Due to the computer vision and image sensor technology used, the Telraam S2
% can only recognize and count vehicles and pedestrians during daylight hours.
% The script attempts to adjust for the progressive truncation of traffic data
% at sunset during the fall and winter by determining the ratio of total
% counts to counts registered before 3pm during months when the sun sets
% late.  To compute an 'adjusted' count that attempts to correct for early
% fall sunsets, the count before 3pm is multiplied by the ratio and used as
% a surrogate count if it exceeds the total of raw counts on a given day.
%
% In practice, the 'adjusted' counts only exceed the total 'raw' counts
% during a few weeks of the fall, expecially after the clocks are set back
% in October.
%
% Requires Telraam Export files with hourly intervals, and '10 modes'
% support provided with data subscription.

clear all
close all
clc

%% User Input variables

% Filename stem for Telraam export file after conversion to Matlab format

%locationString = 'Western Segment';
locationString = 'Eastern Segment';

if strcmp(locationString,'Western Segment')
    fileStem2024 = 'telraam-raw-data-9000007290-2024West60-a0fb917'; % Counts up to end of 2024
    fileStem2025 = 'telraam-raw-data-13374-West60-401b954';
elseif strcmp(locationString,'Eastern Segment')
    fileStem2024 = 'telraam-raw-data-13690-2024East60-f809348'; % Counts up to end of 2024
    fileStem2025 = 'telraam-raw-data-13690-East60-feb4c7'; locationString = 'Eastern Segment';
else
    disp(['Unknown value ' locationString ' for locationString']);
    return;
end

% Counts for Dorion location - requires further modification to title display
%fileStem2024 = 'telraam-raw-data-13938-2024Dorion60-93abd75'; locationString = 'Dorion';
%fileStem2025 = 'telraam-raw-data-13938-Dorion60-6b295d5';

loadMatlabFile = false; % If false load Excel file directly

uptimeThreshold = 0.0; % All time bins below this will be excluded from any analysis (0-0.7)

maxUptimeCorrection = 1.0; % Maximum correction factor - 60 minute bins are already corrected

filterWeekTotalsByDayType = false;

dayTypeString = 'Weekday';
%dayTypeString = 'Weekend';

includePartialMonths = true; % Include months with only partial data coverage

% String value to indicate which modality column from the input table should be analyzed
% Uncomment desired value
modeString = 'BicycleTotal_10Modes_'; modeDisplayString = 'Bike Counts';
%modeString = 'Uptime'; modeDisplayString = 'Uptime in Hours'; uptimeThreshold = 0;
%modeString = 'PedestrianTotal'; modeDisplayString = 'Pedestrian Counts';
%modeString = 'StrollerTotal_10Modes_'; modeDisplayString = 'Stroller Counts';
%modeString = 'CarTotal';  modeDisplayString = 'Car Counts';

% Start and end times applied to filter input data
startTime = datetime(2024,08,01,00,00,01); % Counting start
%startTime = datetime(2024,11,16,00,00,01); includePartialMonths = false; % Winter start
%startTime = datetime(2025,02,01,00,00,01); includePartialMonths = false; % Month start
%startTime = datetime(2025,03,10,00,00,01); includePartialMonths = false; % Week start
%endTime = datetime(2025,03,31,23,59,59); % Winter end
endTime = datetime(2025,06,01,23,59,59); % Week end (Sunday night)

computeDaylightCorrection = false;

%% General plot parameters

plotLineWidth = 10.0;
axisFontSize = 16.0;
labelFontSize = 20.0;
titleFontSize = 24.0;
legendFontSize = 16.0;
axisBackgroundColor = 0.8.*[1 1 1];
legendBackgroundAlpha = 0.2;

%% Load data

if loadMatlabFile

    matFileName = [fileStem2024 '.mat'];
    matVarName = strrep(fileStem2024,'-','_'); % Convert file name to legal matlab variable name
    inputTable = open(matFileName);
    inputTable = inputTable.(matVarName);

else

    % Read file
    excelFileName2024 = [fileStem2024 '.xlsx'];
    %inputTable2024 = readtable(excelFileName2024, 'Sheet', 'Worksheet instances');
    inputTable2024 = readtable(excelFileName2024);

    % Rename variables
    inputTable2024 = renamevars(inputTable2024,'DateAndTime_Local_','Date');
    inputTable2024 = renamevars(inputTable2024,'SpeedV85Km_h','SpeedV85');

    % Compute missing variables
    inputTable2024.Uptime = ones(size(inputTable2024.Uptime));
    inputTable2024.NightTotal_10Modes_ = inputTable2024.ModeNight_A_B_ + inputTable2024.ModeNight_B_A_;
    inputTable2024.BicycleTotal_10Modes_ = inputTable2024.ModeBicycle_A_B_ + inputTable2024.ModeBicycle_B_A_ ;

    % Convert dates
    inputTable2024.Date = datetime(char(inputTable2024.Date),'Format','yy-MM-dd HH:mm');
    inputTable2024 = inputTable2024(year(inputTable2024.Date)==2024,:);

    % ---------------

    % Read file
    excelFileName2025 = [fileStem2025 '.xlsx'];
    %inputTable2025 = readtable(excelFileName2025, 'Sheet', 'Worksheet instances');
    inputTable2025 = readtable(excelFileName2025);

    % Rename variables
    inputTable2025 = renamevars(inputTable2025,'DateAndTime_Local_','Date');
    inputTable2025 = renamevars(inputTable2025,'SpeedV85Km_h','SpeedV85');

    % Compute missing variables
    inputTable2025.Uptime = ones(size(inputTable2025.Uptime));
    inputTable2025.NightTotal_10Modes_ = inputTable2025.ModeNight_A_B_ + inputTable2025.ModeNight_B_A_;
    inputTable2025.BicycleTotal_10Modes_ = inputTable2025.ModeBicycle_A_B_ + inputTable2025.ModeBicycle_B_A_ ;

    % Convert dates
    inputTable2025.Date = datetime(char(inputTable2025.Date),'Format','yy-MM-dd HH:mm');
    inputTable2025 = inputTable2025(year(inputTable2025.Date)==2025,:);

    % ----------------

    % Only take columns we need, to avoid issues with Telraam naming
    % instability
    columnsToKeep = {'Uptime','Date','BicycleTotal_10Modes_','PedestrianTotal','NightTotal_10Modes_','SpeedV85','CarTotal'};
    inputTable2024 = inputTable2024(:,columnsToKeep);
    inputTable2025 = inputTable2025(:,columnsToKeep);

    % Concatenate tables
    inputTable = vertcat(inputTable2024,inputTable2025);

end

inputTable = table2timetable(inputTable);

% Apply uptime threshold to input table
%inputTable = inputTable(inputTable.Uptime>=uptimeThreshold,:);

% Initialize weekday/weekend lists used later as filtering criteria
weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
weekends = {'Saturday', 'Sunday'};

% Apply date range to input table
inputTable = inputTable(find(inputTable.Date>=startTime & inputTable.Date<=endTime),:);

% Add other descriptive columns
inputTable.dayOfWeek = string(day(inputTable.Date,'name'));
inputTable.dayOfWeekCat = categorical(inputTable.dayOfWeek);
inputTable.weekOfYear = week(inputTable.Date,'iso-weekofyear');
inputTable.isWeekday = ismember(inputTable.dayOfWeekCat,weekdays);
%inputTable.weekStartDateTimes = dateshift(inputTable.Date,'start','week')+day(1); % Start dates for each week starting on Monday
inputTable.weekStartDateTimes = dateshift(dateshift(inputTable.Date,'dayofweek','Monday','previous'),'start','day'); % Start dates for each week starting on Monday
inputTable.Daylight = ~((inputTable.NightTotal_10Modes_ > 0) | (isnan(inputTable.SpeedV85)));
inputTable.DaylightUptime = inputTable.Daylight .* inputTable.Uptime;

%% Get weather data for weeks

% Weeks are numbered from 1-53 and always encompass a block of days running from Monday to Sunday.
% Under iso-weekofyear, week one starts on the Monday of the week in which New Year's Day falls,
% even if the Monday and other weekdays are in December.  This makes it harder to group weeks in
% sequence by week number, so all seven days starting on the last Monday in December are assigned
% a week number of 53.
%
% Conversely, iso-weekofyear may assign week=1 to the last days in December.  We also modify these
% to have a week number of 53.
%
% After modification, week 1 always starts on the first Monday in January.
% Similarly, week 53 always starts on the last Monday in December 

inputTable.weekOfYear( (month(inputTable.Date)==12) & (inputTable.weekOfYear==1) ) = 53; % First week of year must start in January, not December
inputTable.weekOfYear( (month(inputTable.Date)==1) & (inputTable.weekOfYear==1) ) = 53; % January first may fall in last week of year
inputTable.yearOfMondayInWeek = year(inputTable.weekStartDateTimes);
januaryIndicesToChange = (month(inputTable.weekStartDateTimes)==1) & (inputTable.weekOfYear==53);
inputTable.yearOfMondayInWeek( januaryIndicesToChange ) = inputTable.yearOfMondayInWeek(januaryIndicesToChange) - 1;
inputTable.yearWeekKey = inputTable.yearOfMondayInWeek + inputTable.weekOfYear./100;

% Generate a list with one datetime object per day, with time set to noon.
% This will be used to retrieve daily weather values for averaging by week.
uniqueWeatherDays = unique(dateshift(inputTable.Date, 'start', 'day')); % Collapse all rows to get one per day
dailyNoonTimes = uniqueWeatherDays + hours(12); % List of unique days with time set to noon

uniqueWeatherDaysYear = year(uniqueWeatherDays);
uniqueWeatherDaysMonth = month(uniqueWeatherDays);
uniqueWeatherDaysWeek = week(uniqueWeatherDays,'iso-weekofyear');

uniqueWeatherDaysWeek( uniqueWeatherDaysMonth==12 & uniqueWeatherDaysWeek==1 ) = 53;
uniqueWeatherDaysWeek( uniqueWeatherDaysMonth==1 & uniqueWeatherDaysWeek==1 ) = 53;
januaryIndicesToChange = ( uniqueWeatherDaysMonth==1 & uniqueWeatherDaysWeek==53 );
uniqueWeatherDaysYear(januaryIndicesToChange) = uniqueWeatherDaysYear(januaryIndicesToChange) - 1;
uniqueYearWeekKey = uniqueWeatherDaysYear + uniqueWeatherDaysWeek./100;

% Retrieve data from server
disp('Getting detailed weather data...')
[precipitationData,temperatureData,sunriseData, sunsetData, sunhoursData,snowData,windspeedData, feelslikeData] = getWeatherstackData('Montreal',dailyNoonTimes);
disp('Done')

[weeksForWeather, idx0, idxWeeks] = unique(uniqueYearWeekKey); % Get list of weeks covered by data for weekly weather stats
meanTemperatureByWeek = accumarray(idxWeeks, temperatureData, [], @mean);
meanSunhoursByWeek = accumarray(idxWeeks, sunhoursData, [], @mean);
minTemperatureByWeek = accumarray(idxWeeks, temperatureData, [], @min);
maxTemperatureByWeek = accumarray(idxWeeks, temperatureData, [], @max);
totalPrecipitationByWeek = accumarray(idxWeeks, precipitationData, [], @sum);
meanWindspeedByWeek = accumarray(idxWeeks, windspeedData, [], @mean);
meanFeelslikeByWeek = accumarray(idxWeeks, feelslikeData, [], @mean);
totalSnowByWeek = accumarray(idxWeeks, snowData, [], @sum);
sunriseByWeek = cellstr(datestr(cellfun(@timeofday,sunriseData(idx0)),'hh:MM'));
sunsetByWeek = cellstr(datestr(cellfun(@timeofday,sunsetData(idx0)),'hh:MM'));
sunsetDateTimeByWeek = sunsetData(idx0);
sunsetDateTimeByWeek = [sunsetDateTimeByWeek{:}];
%uniqueWeeksForDisplay = datetime(floor(unique(weeksForWeather)),01,01,00,00,01) + days(mod(unique(weeksForWeather),1).*7.*100);
uniqueWeeksForDisplay = weekNumToMonday(floor(unique(weeksForWeather)),mod(unique(weeksForWeather),1).*100);

%% Compute daylight correction

% Compute uptime correction factor for all hourly bins
uptimeCorrection = 1./inputTable.Uptime;
uptimeCorrection(uptimeCorrection>maxUptimeCorrection) = maxUptimeCorrection; % Clamp to maximum value
inputTable.AdjustedCountsUptime = inputTable.(modeString) .* uptimeCorrection; % Adjust every hourly bin

% Cutoff times for reliable daylight counts and dark mode
truncationCutoffTime = timeofday(datetime('today')+hours(15));
sunsetCutoffTime = timeofday(datetime('today')+hours(18.5));

% Use different correction factors for weekdays and weekend days
inputTableWD = inputTable(inputTable.isWeekday,:);
inputTableWE = inputTable(~inputTable.isWeekday,:);

if computeDaylightCorrection

    % Group and sum the total uptime-adjusted counts for weekdays and weekends
    groupWeekTableAllWD = groupsummary(inputTableWD,'weekStartDateTimes','sum','AdjustedCountsUptime');
    weeklySumsAllWD = groupWeekTableAllWD.('sum_AdjustedCountsUptime');
    groupWeekTableAllWE = groupsummary(inputTableWE,'weekStartDateTimes','sum','AdjustedCountsUptime');
    weeklySumsAllWE = groupWeekTableAllWE.('sum_AdjustedCountsUptime');

    % Now extract only bins that would reliably be registered during full daylight
    inputTableDaylight = inputTable(timeofday(inputTable.Date)<=truncationCutoffTime,:);

    % Create weekday and weekend tables
    inputTableDaylightWD = inputTableDaylight(inputTableDaylight.isWeekday,:);
    inputTableDaylightWE = inputTableDaylight(~inputTableDaylight.isWeekday,:);

    % Group and sum uptime-adjusted counts for daylight tables
    groupWeekTableDaylightWD = groupsummary(inputTableDaylightWD,'weekStartDateTimes','sum','AdjustedCountsUptime');
    weeklySumsDaylightWD = groupWeekTableDaylightWD.('sum_AdjustedCountsUptime');
    groupWeekTableDaylightWE = groupsummary(inputTableDaylightWE,'weekStartDateTimes','sum','AdjustedCountsUptime');
    weeklySumsDaylightWE = groupWeekTableDaylightWE.('sum_AdjustedCountsUptime');

    % Compute ratios of all counts to daylight-only counts
    weeklyRatiosWD = weeklySumsAllWD ./ weeklySumsDaylightWD;
    weeklyRatiosWE = weeklySumsAllWE ./ weeklySumsDaylightWE;

    % Filter to days where sunset is after PM rush hour peak
    longDayIndices = find(timeofday(sunsetDateTimeByWeek)>sunsetCutoffTime);

    % Use mean ratio during longer days as correction factor
    daylightCorrectionRatioWD = mean(weeklyRatiosWD(longDayIndices));
    daylightCorrectionRatioWE = mean(weeklyRatiosWE(longDayIndices));

else

    daylightCorrectionRatioWD = 1.65;
    daylightCorrectionRatioWE = 1.37;

end

% This will apply to all time bins but only 'daylight' bins are to be used for estimated totals
inputTableWD.AdjustedCountsUptimeDaylight = inputTableWD.AdjustedCountsUptime .* daylightCorrectionRatioWD;
inputTableWE.AdjustedCountsUptimeDaylight = inputTableWE.AdjustedCountsUptime .* daylightCorrectionRatioWE;

% New inputTable will have column with correction factor specific to day type;  only valid for counts up to 3PM
inputTable = sortrows([inputTableWD;inputTableWE],'Date');  

if computeDaylightCorrection

    figure('Position',[408 126 1132 921]);

    yyaxis left

    hold on

    h1 = plot(uniqueWeeksForDisplay,weeklyRatiosWD,'b-','LineWidth',plotLineWidth,'DisplayName','Weekday Ratios');
    h2 = plot([uniqueWeeksForDisplay(1) uniqueWeeksForDisplay(longDayIndices(end))],daylightCorrectionRatioWD.*[1 1],'b:','LineWidth',plotLineWidth./3,'DisplayName','Correction Ratio Weekdays');

    h3 = plot(uniqueWeeksForDisplay(1:length(weeklyRatiosWE)),weeklyRatiosWE,'r-','LineWidth',plotLineWidth,'DisplayName','Weekend Ratios');
    h4 = plot([uniqueWeeksForDisplay(1) uniqueWeeksForDisplay(longDayIndices(end))],daylightCorrectionRatioWE.*[1 1],'r:','LineWidth',plotLineWidth./3,'DisplayName','Correction Ratio Weekdends');

    h5 = plot(xlim,[1 1],'k:','LineWidth',plotLineWidth./3,'DisplayName','Zero counts registered after 3PM');

    hold off
    ylim([0.9 2]);
    grid on

    xlabel('Date','FontSize',labelFontSize)
    ylabel('Ratio of Counts up to 3PM to All-Day Counts','FontSize',labelFontSize)

    yyaxis right

    h6 = plot(uniqueWeeksForDisplay,timeofday(sunsetDateTimeByWeek),'-','LineWidth',plotLineWidth,'Color',[0 1 1 0.3],'DisplayName','Time of Sunset');
    ylabel('Time of Sunset','FontSize',labelFontSize,'Color','c')

    ax = gca;

    ax.YAxis(1).Color = 'b';
    ax.YAxis(2).Color = 0.7.*[0 1 1];

    set(gca,'Color',axisBackgroundColor);

    title('Estimation of Telraam S2 Counts Lost After Sunset','FontSize',titleFontSize);

    legend([h1 h2 h3 h4 h5 h6],'FontSize',legendFontSize);

end

%% Compute daily traffic patterns of selected mode and cars on weekdays and weekends

% Weekdays:
weekdayTraffic = [];
weekdayCarTraffic = [];
weekdayTimes = [];

for ix = 1:length(weekdays) % Loop over weekdays

    targetDay = weekdays{ix};

    % Extract rows that match the target day for weekdays
    dayRows = strcmp(string(inputTable.dayOfWeek), targetDay);

    % Filter the timetable for the selected day
    filteredTable = inputTable(dayRows, :);

    % Extract the time component and add to the weekday traffic arrays
    timeOfDay = timeofday(filteredTable.Date);
    weekdayTraffic = [weekdayTraffic; filteredTable.(modeString)];
    weekdayCarTraffic = [weekdayCarTraffic; filteredTable.('CarTotal')];
    weekdayTimes = [weekdayTimes; timeOfDay];

end

% Weekends:
weekendTraffic = [];
weekendCarTraffic = [];
weekendTimes = [];

for ix = 1:length(weekends) % Loop over weekend days

    targetDay = weekends{ix};
    
    % Extract rows that match the target day for weekends
    dayRows = strcmp(string(inputTable.dayOfWeek), targetDay);
    
    % Filter the timetable for the selected day
    filteredTable = inputTable(dayRows, :);
    
    % Extract the time component and add to the weekend traffic arrays
    timeOfDay = timeofday(filteredTable.Date);
    weekendTraffic = [weekendTraffic; filteredTable.(modeString)];
    weekendCarTraffic = [weekendCarTraffic; filteredTable.('CarTotal')];
    weekendTimes = [weekendTimes; timeOfDay];

end

%% Plot daily traffic patterns of selected mode on week along with that of car traffic

figure('Position',[408 126 1132 921]);
hold on;

if strcmp(dayTypeString,'Weekday')

    % Calculate the average traffic for each unique time of day for weekdays
    [uniqueWeekdayTimes, ~, idx] = unique(weekdayTimes);
    avgWeekdayTraffic = accumarray(idx, weekdayTraffic, [], @mean);
    avgWeekdayCarTraffic = accumarray(idx, weekdayCarTraffic, [], @mean);

    % Compute totals for weekdays
    weekdayTotal = sum(avgWeekdayTraffic);
    weekdayCarTotal = sum(avgWeekdayCarTraffic);

    plot(uniqueWeekdayTimes, avgWeekdayTraffic, '-', 'DisplayName', [ modeDisplayString ' (total = ' num2str(weekdayTotal,'%.0f') ' per day)'],'LineWidth',plotLineWidth);
    plot(uniqueWeekdayTimes, avgWeekdayCarTraffic, '-', 'DisplayName', ['Car Traffic (total = ' num2str(weekdayCarTotal,'%.0f') ' per day)'],'LineWidth',plotLineWidth);
    title([ 'Average Car and ' modeDisplayString ' on Weekdays ( ' locationString ' )'],'FontSize',titleFontSize);

elseif strcmp(dayTypeString,'Weekend')

    % Calculate the average traffic for each unique time of day for weekends
    [uniqueWeekendTimes, ~, idx] = unique(weekendTimes);
    avgWeekendTraffic = accumarray(idx, weekendTraffic, [], @mean);
    avgWeekendCarTraffic = accumarray(idx, weekendCarTraffic, [], @mean);

    % Compute totals for weekends
    weekendTotal = sum(avgWeekendTraffic);
    weekendCarTotal = sum(avgWeekendCarTraffic);

    plot(uniqueWeekdayTimes, avgWeekdayTraffic, '-', 'DisplayName', [ modeDisplayString ' (total = ' num2str(weekdayTotal,'%.0f') ' per day)'],'LineWidth',plotLineWidth);
    plot(uniqueWeekdayTimes, avgWeekdayCarTraffic, '-', 'DisplayName', ['Car Traffic (total = ' num2str(weekdayCarTotal,'%.0f') ' per day)'],'LineWidth',plotLineWidth);
    title([ 'Average Car and ' modeDisplayString ' on Weekends ( ' locationString ' )'],'FontSize',titleFontSize);

end

% Set axis labels and title
xlabel('Time of Day','FontSize',labelFontSize);
ylabel('Hourly Count','FontSize',labelFontSize);
set(gca,'Color',0.8.*[1 1 1])
grid on;

legend('show');

hold off

%% Plot daily traffic patterns of selected mode on weekdays and weekends

figure('Position',[408 126 1132 921]);
hold on

% Calculate the average traffic for each unique time of day for weekdays
[uniqueWeekdayTimes, ~, idx] = unique(weekdayTimes);
avgWeekdayTraffic = accumarray(idx, weekdayTraffic, [], @mean);

% Calculate the average traffic for each unique time of day for weekends
[uniqueWeekendTimes, ~, idx] = unique(weekendTimes);
avgWeekendTraffic = accumarray(idx, weekendTraffic, [], @mean);

% Plot the average traffic for weekdays
weekdayTotal = sum(avgWeekdayTraffic);
plot(uniqueWeekdayTimes, avgWeekdayTraffic, '-', 'DisplayName', ['Weekdays (total = ' num2str(weekdayTotal,'%.0f') ' per day)'],'LineWidth',plotLineWidth);

% Plot the average traffic for weekends
weekendTotal = sum(avgWeekendTraffic);
plot(uniqueWeekendTimes, avgWeekendTraffic, '-', 'DisplayName', ['Weekends (total = ' num2str(weekendTotal,'%.0f') ' per day)'],'LineWidth',plotLineWidth);

% Set axis labels and title
xlabel('Time of Day','FontSize',labelFontSize);
ylabel('Hourly Count','FontSize',labelFontSize);
title([ 'Average ' modeDisplayString ' on Weekdays and Weekend Days ( ' locationString ' )'],'FontSize',titleFontSize);
set(gca,'Color',0.8.*[1 1 1])
grid on;

legend('show');

hold off

%% Plot average hourly counts versus time of day for each month
% Also generates canonical daily traffic curves for weekdays and weekends
% to use later when correcting for truncation by sunset

figure('Position', [408 126 1132 921]);

% Extract unique months from the inputTable
inputTable.monthStartDateTimes = dateshift(inputTable.Date, 'start', 'month');
uniqueMonthStartDates = unique(inputTable.monthStartDateTimes);

% Filter to exclude partial months if the flag is set
if ~includePartialMonths
    monthCounts = groupcounts(unique(inputTable.monthStartDateTimes));
    completeMonths = uniqueMonthStartDates(monthCounts == max(monthCounts)); % Only include months with max records
    uniqueMonthStartDates = completeMonths;
end

% Initialize variables for plotting
numMonths = length(uniqueMonthStartDates);
colorMap = lines(numMonths);

% Select day type for plotting
if dayTypeString == "Weekday"
    dayType = weekdays; 
    plotStr = 'b-';  
    dayString = 'Weekday';
elseif dayTypeString == "Weekend"
    dayType = weekends; 
    plotStr = 'r-'; 
    dayString = 'Weekend';
else
    error(['Unknown value (' dayTypeString ') of dayTypeString']);
end

% Loop over each unique month and calculate average traffic counts for plotting

uptimeByWeek = [];

for i = 1:numMonths

    % Filter data for the current month
    currentMonthRows = isbetween(inputTable.Date, uniqueMonthStartDates(i), uniqueMonthStartDates(i) + calmonths(1) - seconds(1));
    currentMonthTable = inputTable(currentMonthRows, :);

    nonZeroUptimesMonthTable = currentMonthTable(currentMonthTable.(modeString)>0,:);
    nonZeroUptimesMonth = nonZeroUptimesMonthTable.Uptime;
    uptimeByMonth(i) = 100.*sum(nonZeroUptimesMonth)./length(nonZeroUptimesMonth);
    
    % Filter for the selected day type (weekdays or weekends)
    dayTypeTable = currentMonthTable(ismember(currentMonthTable.dayOfWeek, dayType), :);
    
    % Extract time of day and traffic counts
    timeOfDay = timeofday(dayTypeTable.Date);
    %trafficCounts = dayTypeTable.(modeString);
    trafficCounts = dayTypeTable.AdjustedCountsUptimeDaylight;
    
    % Calculate the average traffic for each unique time of day in the current month
    [uniqueTimes, ~, idx] = unique(timeOfDay);
    avgTraffic = accumarray(idx, trafficCounts, [], @mean);

    % Plot the average traffic for this month
    monthTotal = sum(avgTraffic);
    plot(uniqueTimes, avgTraffic, '-', 'DisplayName', ...
         ['Month starting ' datestr(uniqueMonthStartDates(i), 'dd-mmm') ' (total = ' num2str(monthTotal, '%.0f') ' per day)'], ...
         'LineWidth', 8.0, 'Color', [colorMap(i, :) 0.8]);

    hold on

end

% Add labels and adjust plot format
xlabel('Time of Day', 'FontSize', labelFontSize);
ylabel('Hourly Count', 'FontSize', labelFontSize);
%title([dayString ' ' modeDisplayString ' with Monthly Variation (' locationString ')'], 'FontSize', titleFontSize);
title(['Average Hourly ' modeDisplayString ' by Time of Day (' locationString ', ' dayString 's Only)'], 'FontSize', titleFontSize);
set(gca, 'Color', 0.8 * [1 1 1]);
grid on;

legend('show')

hold off

%% Plot average hourly counts versus time of day for each week
% also fits canonical daily traffic curve to weekly averages, so that missing data
% after dark can be estimated in order to compensate for shorter days in the winter.

figure('Position',[408 126 1132 921]);

% Initialize variables for storing average weekday traffic per week
weeklyAverageTraffic = {};  % Will store average traffic for each week
uniqueWeekdayTimesByWeek = {};  % To store time of day for the average traffic

% Get the list of unique weeks (Monday-Sunday) within the selected date range
weekCheckTable = inputTable(inputTable.dayOfWeek=='Monday',:);
weekNumbers = unique(week(weekCheckTable.Date,'iso-weekofyear')); % Find unique week numbers
weekStartDateTimes = unique(dateshift(dateshift(weekCheckTable.Date,'dayofweek','Monday','previous'),'start','day')); % Start dates for each week starting on Monday

if dayTypeString=='Weekday'
    dayType = weekdays; plotStr = 'b-';  avgVector = avgWeekdayTraffic; dayString = 'Weekday'; displayAvg = weekdayTotal;
elseif dayTypeString=='Weekend'
    dayType = weekends; plotStr = 'r-'; avgVector = avgWeekendTraffic; dayString = 'Weekend'; displayAvg = weekendTotal;
else
    disp(['Error:  unknown value (' dayTypeString ') of dayTypeString']);
end

colorMap = lines(length(weekStartDateTimes)); % use consistent color list for plots

correctedWeeklyTraffic = [];
uptimeByWeek = [];

% Loop over each week and calculate the average traffic for weekdays
numWeeks = length(weekStartDateTimes);
for i = 1:numWeeks
    % Filter data for the current week
    currentWeekRows = isbetween(inputTable.Date, weekStartDateTimes(i), weekStartDateTimes(i) + days(6));
    currentWeekTable = inputTable(currentWeekRows, :);    

    % Always corrects by uptime - make consistent/optional !!!
    uptimeCorrection = 1./currentWeekTable.Uptime;
    uptimeCorrection(uptimeCorrection>maxUptimeCorrection) = maxUptimeCorrection;
    currentWeekTable.(modeString) = currentWeekTable.(modeString) .* uptimeCorrection;

    nonZeroUptimesTable = currentWeekTable(currentWeekTable.(modeString)>0,:);
    nonZeroUptimes = nonZeroUptimesTable.Uptime;
    uptimeByWeek(i) = 100.*sum(nonZeroUptimes)./length(nonZeroUptimes);
    
    % Filter for weekdays only in this week
    weekdayTable = currentWeekTable(ismember(currentWeekTable.dayOfWeek, weekdays), :);
    timeOfDayWD = timeofday(weekdayTable.Date);
    trafficCountsWD = weekdayTable.(modeString);
    
    % Calculate the average traffic for each unique time of day in this week
    [uniqueTimesWD, ~, idxWD] = unique(timeOfDayWD);
    avgTrafficWD = accumarray(idxWD, trafficCountsWD, [], @mean);
    
    % Plot the average traffic for this specific week

    % Filter for weekends now
    weekendTable = currentWeekTable(ismember(currentWeekTable.dayOfWeek, weekends), :);
    timeOfDayWE = timeofday(weekendTable.Date);
    trafficCountsWE = weekendTable.(modeString);
    
    % Calculate the average traffic for each unique time of day in this week
    [uniqueTimesWE, ~, idxWE] = unique(timeOfDayWE);
    avgTrafficWE = accumarray(idxWE, trafficCountsWE, [], @mean);

    % Plot the average traffic for this specific week

    if strcmp(dayTypeString,'Weekday')
        uniqueTimes = uniqueTimesWD;
        weeklyAverageTraffic{i} = avgTrafficWD;
        uniqueWeedayTimesByWeek{i} = uniqueTimesWD;
        weekTotal = sum(avgTrafficWD); 
        avgTraffic = avgTrafficWD;
    elseif strcmp(dayTypeString,'Weekend')
        uniqueTimes = uniqueTimesWE;
        weeklyAverageTraffic{i} = avgTrafficWE;
        uniqueWeedayTimesByWeek{i} = uniqueTimesWE;
        weekTotal = sum(avgTrafficWE); 
        avgTraffic = avgTrafficWE;
    end

    % Now plot the weekly average traffic curve for display
    
    if i<numWeeks % General weekly curves use a subtle display
        plot(uniqueTimes, avgTraffic, '-', 'DisplayName', ['Week starting ' datestr(weekStartDateTimes(i), 'dd-mmm') ' (total = ' num2str(weekTotal,'%.0f') ' per day)'], 'LineWidth', 5.0, 'Color', [colorMap(i, :) 0.1]);
    else % The most recent weekly curve is emphasized
        plot(uniqueTimes, avgTraffic, '-', 'DisplayName', ['Week starting ' datestr(weekStartDateTimes(i), 'dd-mmm') ' (total = ' num2str(weekTotal,'%.0f') ' per day)'], 'LineWidth', 5.0, 'Color', [colorMap(i, :) 1.0]);
    end

    hold on;
end

% Plot the grand average daily traffic curve
plot(uniqueWeekdayTimes, avgVector, plotStr, 'DisplayName', ['Average over Weeks (total = ' num2str(displayAvg,'%.0f') ' per day)'], 'LineWidth', 10.0);

% Add labels and ajust format
xlabel('Time of Day','FontSize',labelFontSize);
ylabel('Hourly Count','FontSize',labelFontSize);
title(['Average Hourly ' modeDisplayString ' by Time of Day (' locationString ', ' dayString 's Only)'], 'FontSize', titleFontSize);
set(gca,'Color',0.8.*[1 1 1])
grid on;

legend('show');

hold off

%% Plot average mode traffic on each weekday

figure('Position',[408 126 1132 921]);

% Initialize arrays to store traffic data for each day of the week
daysOfWeek = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'};
weeklyDayAverageTraffic = {};  % Store weekly averages for each day
weeklyDayTrafficCounts = zeros(7, length(weekStartDateTimes));  % Store traffic counts for each week (7 days)

% Loop over each week and calculate the average traffic for each day
for i = 1:length(weekStartDateTimes)
    % Filter data for the current week
    currentWeekRows = isbetween(inputTable.Date, weekStartDateTimes(i), weekStartDateTimes(i) + days(7) - seconds(1) );
    currentWeekTable = inputTable(currentWeekRows, :);
    
    % Initialize array for storing average traffic per day in the current week
    avgTrafficPerDay = zeros(1, 7);
    
    % Loop over each day of the week (Monday to Sunday)
    for j = 1:length(daysOfWeek)
        % Filter for the current day of the week
        dayRows = strcmp(string(currentWeekTable.dayOfWeek), daysOfWeek{j});
        
        % Check if the filtering worked and there is data for this day
        if any(dayRows)
            dayTable = currentWeekTable(dayRows, :);
            
            % Calculate the total traffic for this day in the current week
            avgTrafficPerDay(j) = sum(dayTable.(modeString)); % Sum not average
            if avgTrafficPerDay(j)>24
                %return; % !!! DEBUG
            end
        else
            % Print a warning if no data is found for a particular day (for debugging)
            fprintf('Warning: No data found for %s in week starting %s\n', daysOfWeek{j}, datestr(weekStartDateTimes(i), 'dd-mmm-yyyy'));
        end
    end
    
    % Store the average traffic for this week in the cell array
    weeklyDayAverageTraffic{i} = avgTrafficPerDay;
    weeklyDayTrafficCounts(:, i) = avgTrafficPerDay;  % For averaging later
    
    weekTotal = sum(avgTrafficPerDay);

    % Plot the average traffic for this week with the color scheme
    plot(1:7, avgTrafficPerDay, '-', 'DisplayName', ['Week starting ' datestr(weekStartDateTimes(i), 'dd-mmm') ' ( total = ' num2str(weekTotal) ' )'], ...
        'LineWidth', 5.0, 'Color', [colorMap(i, :) 0.2]);
    hold on;
end

% Calculate the average traffic for each day across all weeks
avgWeekTrafficPerDay = mean(weeklyDayTrafficCounts, 2);

% Plot the overall average traffic for each day of the week
plot(1:7, avgWeekTrafficPerDay, 'k-', 'DisplayName', ['Average Week ( total = ' num2str(sum(avgWeekTrafficPerDay)) ' )'], 'LineWidth', 10.0);

% Set plot labels and formatting
xticks(1:7);
xticklabels(daysOfWeek);  % Label the x-axis with day names
xlabel('Day of the Week', 'FontSize', labelFontSize);
ylabel([modeDisplayString], 'FontSize', labelFontSize);
title([ modeDisplayString ' per Day of the Week ( ' locationString ' )'], 'FontSize', titleFontSize);
set(gca,'Color',0.8.*[1 1 1])  % Light background
grid on;
set(gca, 'GridColor', [0 0 0], 'GridAlpha', 0.2);  % Darker grid lines

yl = ylim;
upperMarginFactor = 1.1;
ylim([0 upperMarginFactor.*yl(2)]);

legend('show')

hold off

%% Plot total mode traffic in each week

if numWeeks > 1

    if filterWeekTotalsByDayType
        if strcmp(dayTypeString,'Weekday')
            inputTable = inputTable(inputTable.isWeekday,:);
        elseif strcmp(dayTypeString,'Weekend')
            inputTable = inputTable(~inputTable.isWeekday,:);
        else
            disp(['Unknown value ' dayTypeString ' of dayTypeString']);
        end
    end

    % Calculate total traffic per week
    filteredTable = inputTable(ismember(inputTable.weekStartDateTimes,weekStartDateTimes),:);
    groupTrafficTable = groupsummary(filteredTable,'weekStartDateTimes','sum',modeString);
    totalWeeklyTraffic = groupTrafficTable.(['sum_' modeString]);

    % Get both min and max daily counts for the specified weeks
    [minWeeklyTraffic, maxWeeklyTraffic] = calculateWeeklyDailyMinMax(inputTable, modeString, weekStartDateTimes);

    % Calculate total traffic per week excluding PM rush hour
    filteredTable2 = inputTable(ismember(inputTable.weekStartDateTimes,weekStartDateTimes),:);
    filteredTable2 = filteredTable2(timeofday(filteredTable2.Date)<=truncationCutoffTime,:);
    groupTrafficTable2 = groupsummary(filteredTable2,'weekStartDateTimes','sum',modeString);
    totalWeeklyTraffic2 = groupTrafficTable2.(['sum_' modeString]);

    groupAdjustedTrafficTable2 = groupsummary(filteredTable2,'weekStartDateTimes','sum','AdjustedCountsUptimeDaylight');
    totalAdjustedWeeklyTraffic2 = groupAdjustedTrafficTable2.('sum_AdjustedCountsUptimeDaylight');

    totalAdjustedWeeklyTraffic2 = max(totalAdjustedWeeklyTraffic2,totalWeeklyTraffic);

    totalCount = sum(totalWeeklyTraffic);
    minCount = min(totalWeeklyTraffic);
    maxCount = max(totalWeeklyTraffic);
    totalCount2 = round(sum(totalAdjustedWeeklyTraffic2),0);
    minCount2 = round(min(totalAdjustedWeeklyTraffic2),0);
    maxCount2 = round(max(totalAdjustedWeeklyTraffic2),0);

    % Plotting the total bike traffic per week
    figure('Position', [408 126 1132 921]);

    blueColor = [0 0.4471 0.7412 0.3];

    % Plot average weekly temperature and other information using a separate axis

    yyaxis right
    uptimeLineColor = [1 0 0 0.1];
    uptimeLabelColor = [1 0 0 0.5]+0.5.*[0 1 1 0];
    temperatureLineColor = blueColor;
    daylightLineColor = [ 0.4.*[1 1 1] 0.2];
    hh0 = plot(uniqueWeeksForDisplay, meanTemperatureByWeek,'LineWidth',plotLineWidth,'Color',temperatureLineColor,'LineStyle','-','DisplayName','Temperature (°C)');
    hold on
    numWeeks = length(uptimeByWeek);
    %listVec = [1:numWeeks];
    listVec = [2:length(uniqueWeeksForDisplay)]; % need to match number of traffic counts for tooltips
    hh00 = bubblechart(uniqueWeeksForDisplay(listVec), meanTemperatureByWeek(listVec),totalPrecipitationByWeek(listVec),'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerFaceAlpha',0.3,'DisplayName',...
        'Marker Size Scaled by Total Weekly Precipitation');
    hold off
    ylabel('Mean Temperature for Week (°C)', 'FontSize', labelFontSize+2,'Color',uptimeLabelColor,'FontWeight','bold');
    ylim([-15 max(ylim).*1.1]);

    % Set up datatips for additional info
    hh00.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Location:',repmat({locationString}, 1, length(listVec)).');
    hh00.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Date:',uniqueWeeksForDisplay(listVec));
    hh00.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Raw Weekly Counts:',num2sepstr(totalWeeklyTraffic));
    hh00.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Min Counts/Day:',num2sepstr(minWeeklyTraffic));
    hh00.DataTipTemplate.DataTipRows(5) = dataTipTextRow('Max Counts/Day:',num2sepstr(maxWeeklyTraffic));
    hh00.DataTipTemplate.DataTipRows(6) = dataTipTextRow('Adjusted Counts:',num2sepstr(totalAdjustedWeeklyTraffic2,'%.0f'));
    hh00.DataTipTemplate.DataTipRows(7) = dataTipTextRow('Mean Temperature (°C):',round(meanTemperatureByWeek(listVec),1));
    hh00.DataTipTemplate.DataTipRows(8) = dataTipTextRow('Mean Windspeed (km/h):',round(meanWindspeedByWeek(listVec),1));
    hh00.DataTipTemplate.DataTipRows(9) = dataTipTextRow('Mean Feels Like Temp (°C):',round(meanFeelslikeByWeek(listVec),1));
    hh00.DataTipTemplate.DataTipRows(10) = dataTipTextRow('Min Temperature (°C):',minTemperatureByWeek(listVec));
    hh00.DataTipTemplate.DataTipRows(11) = dataTipTextRow('Max Temperature (°C):',maxTemperatureByWeek(listVec));
    hh00.DataTipTemplate.DataTipRows(12) = dataTipTextRow('Total Precipitation (mm):',totalPrecipitationByWeek(listVec));
    hh00.DataTipTemplate.DataTipRows(13) = dataTipTextRow('Total Snow (mm):',totalSnowByWeek(listVec));
    hh00.DataTipTemplate.DataTipRows(14) = dataTipTextRow('Sunrise:',sunriseByWeek(listVec));
    hh00.DataTipTemplate.DataTipRows(15) = dataTipTextRow('Sunset:',sunsetByWeek(listVec));
    hh00.DataTipTemplate.DataTipRows(16) = dataTipTextRow('Uptime (%):',round(uptimeByWeek,1).');

    % Now plot raw and adjusted traffic counts for each week

    yyaxis left
    adjustedCountColor = [0 0 0 0.3];

    hh1 = plot(weekStartDateTimes, totalAdjustedWeeklyTraffic2,'-','LineWidth',plotLineWidth,'DisplayName',['Adjusted Counts for Week ( min = ' num2sepstr(minCount2,'%.0f') ' ; max = ' num2sepstr(maxCount2,'%.0f') ' )'],'Color',adjustedCountColor);
    hold on

    %hh000 = plot(weekStartDateTimes, zeros(size(weekStartDateTimes)),'-','LineWidth',plotLineWidth./3,'DisplayName','Sasquatch Sightings','Color','m');

    hh2 = plot(weekStartDateTimes, totalWeeklyTraffic,'-','LineWidth',plotLineWidth,'DisplayName',['Raw Counts for Week ( min = ' num2sepstr(minCount,'%.0f') ' ; max = ' num2sepstr(maxCount,'%.0f') ' )'],'Color','k');
    hh3 = plot(weekStartDateTimes, totalWeeklyTraffic2,':','LineWidth',plotLineWidth./3,'DisplayName',['Raw Counts for Week (up to ' datestr(truncationCutoffTime) ')'],'Color','k');
    hold off
    ylabel(['Total ' modeDisplayString ' for Week'], 'FontSize', labelFontSize+2,'FontWeight','bold');
    ylim([0 max(ylim).*1.1]);
    ax2 = gca;
    ax2.FontSize = axisFontSize;
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = blueColor;

    xlim([weekStartDateTimes(1) weekStartDateTimes(end)])
    xticks(weekStartDateTimes);
    xlabel('Week Starting Date', 'FontSize', labelFontSize);
    title(['Weekly ' modeDisplayString ' on Terrebonne (' locationString ')'], 'FontSize', titleFontSize);
    subtitle(['Total estimated counts:  ' num2sepstr(totalCount2) ' ( over ' num2str(length(weekStartDateTimes)) ' weeks between ' datestr(min(filteredTable.Date),'yyyy-mm-dd') ' and ' datestr(max(filteredTable.Date),'yyyy-mm-dd') ' )' ]);
    subtitle(['Total raw counts over ' num2str(length(weekStartDateTimes)) ' weeks between ' datestr(min(filteredTable.Date),'yyyy-mm-dd') ' and ' datestr(max(filteredTable.Date),'yyyy-mm-dd') ':  ' num2sepstr(totalCount,'%.0f') ...
        ' ( ' num2sepstr(totalCount2,'%.0f') ' adjusted )']);
    set(gca, 'Color', 0.8 * [1 1 1]);
    grid on

    legend([hh2 hh1 hh3 hh0 hh00],'Location','southwest');
    bubblelegend('Weekly Precipitation (mm)','Color',0.8.*[1 1 1],'Location','southwest');

    xtickangle(45);

end

%% Calculate total traffic counts per complete calendar month

if numMonths > 2

    % Extract unique months from the input table
    inputTable.monthStartDateTimes = dateshift(inputTable.Date, 'start', 'month');
    uniqueMonthStartDates = unique(inputTable.monthStartDateTimes);

    % Filter to include only complete months
    % Determine complete calendar months based on date range
    startDatesOfMonths = dateshift(inputTable.Date, 'start', 'month');
    endDatesOfMonths = dateshift(inputTable.Date, 'end', 'month');

    % Find unique months in the date range
    uniqueMonthStartDates = unique(startDatesOfMonths);

    % Identify complete months by checking if all days of the month are present
    completeMonths = uniqueMonthStartDates(...
        ismember(uniqueMonthStartDates, dateshift(inputTable.Date, 'start', 'day')) & ...
        ismember(dateshift(uniqueMonthStartDates, 'end', 'month'), dateshift(inputTable.Date, 'start', 'day')));

    % Filter input table for complete months only
    filteredMonthTable = inputTable(ismember(inputTable.monthStartDateTimes, completeMonths), :);
    groupTrafficTableMonths = groupsummary(filteredMonthTable, 'monthStartDateTimes', 'sum', modeString);
    totalMonthlyTraffic = groupTrafficTableMonths.(['sum_' modeString]);

    [minMonthlyTraffic, maxMonthlyTraffic] = calculateMonthlyDailyMinMax(inputTable, modeString, completeMonths);

    % Calculate adjusted counts for uptime

    uptimeCorrection = 1./filteredMonthTable.Uptime;
    uptimeCorrection(uptimeCorrection>maxUptimeCorrection) = maxUptimeCorrection;
    filteredMonthTable.AdjustedCounts = filteredMonthTable.(modeString) .* uptimeCorrection;

    groupAdjustedTrafficTableMonths = groupsummary(filteredMonthTable, 'monthStartDateTimes', 'sum', 'AdjustedCounts');
    totalAdjustedMonthlyTraffic = groupAdjustedTrafficTableMonths.('sum_AdjustedCounts');

    % Calculate total traffic per month excluding PM rush hour
    filteredMonthTable2 = inputTable(ismember(inputTable.monthStartDateTimes, completeMonths), :);
    filteredMonthTable2 = filteredMonthTable2(timeofday(filteredMonthTable2.Date)<=truncationCutoffTime,:);
    groupTrafficTableMonths2 = groupsummary(filteredMonthTable2, 'monthStartDateTimes', 'sum', modeString);
    totalMonthlyTraffic2 = groupTrafficTableMonths2.(['sum_' modeString]);

    groupAdjustedTrafficTableMonths2 = groupsummary(filteredMonthTable2,'monthStartDateTimes','sum','AdjustedCountsUptimeDaylight');
    totalAdjustedMonthlyTraffic2 = groupAdjustedTrafficTableMonths2.('sum_AdjustedCountsUptimeDaylight');

    % Calculate counts corrected for daylight
    totalAdjustedMonthlyTraffic2 = max(totalMonthlyTraffic,totalAdjustedMonthlyTraffic2);

    % Get weather data for months

    % Generate a list with one datetime object per day, with time set to noon
    uniqueWeatherDays = unique(dateshift(filteredMonthTable.Date, 'start', 'day'));
    dailyNoonTimes = uniqueWeatherDays + hours(12);
    
    % Retrieve data from server
    disp('Getting detailed weather data...')
    [precipitationData,temperatureData,sunriseData, sunsetData, sunhoursData, snowData, windspeedData, feelslikeData] = getWeatherstackData('Montreal',dailyNoonTimes);
    disp('Done')

    uniqueWeatherDaysYear = year(uniqueWeatherDays);
    uniqueWeatherDaysMonth = month(uniqueWeatherDays);
    uniqueYearMonthKey = uniqueWeatherDaysYear + uniqueWeatherDaysMonth./100;

    %[monthsForWeather, idx0, idxMonths] = unique(month(dailyNoonTimes)); % Get list of months covered by data for monthly weather stats
    [monthsForWeather, idx0, idxMonths] = unique(uniqueYearMonthKey); % Get list of months covered by data for monthly weather stats
    meanTemperatureByMonth = accumarray(idxMonths, temperatureData, [], @mean);
    meanSunhoursByMonth = accumarray(idxMonths, sunhoursData, [], @mean);
    minTemperatureByMonth = accumarray(idxMonths, temperatureData, [], @min);
    maxTemperatureByMonth = accumarray(idxMonths, temperatureData, [], @max);
    totalPrecipitationByMonth = accumarray(idxMonths, precipitationData, [], @sum);
    meanWindspeedByMonth = accumarray(idxMonths, windspeedData, [], @mean);
    meanFeelslikeByMonth = accumarray(idxMonths, feelslikeData, [], @mean);
    totalSnowByMonth = accumarray(idxMonths, snowData, [], @sum);
    sunriseByMonth = cellstr(datestr(cellfun(@timeofday,sunriseData(idx0)),'hh:MM'));
    sunsetByMonth = cellstr(datestr(cellfun(@timeofday,sunsetData(idx0)),'hh:MM'));
    uniqueMonthsForDisplay = datetime(floor(unique(monthsForWeather)),round(mod(unique(monthsForWeather),1).*100),1);

    indicesForDisplay = find(uniqueMonthsForDisplay>=groupTrafficTableMonths.monthStartDateTimes(1)); % Display same range as traffic

    % Display total counts for all months
    totalCountMonthly = sum(totalMonthlyTraffic);
    totalCountMonthly2 = sum(totalAdjustedMonthlyTraffic2);

    % Plotting the total traffic counts per month
    figure('Position', [408 126 1132 921]);

    blueColor = [0 0.4471 0.7412 0.3];

    yyaxis right
    yyaxis right
    %hm00 = plot(groupTrafficTableMonths.monthStartDateTimes, meanTemperatureByMonth,'LineWidth',plotLineWidth,'Color',temperatureLineColor,'LineStyle','-','DisplayName','Temperature (°C)');
    hm00 = plot(uniqueMonthsForDisplay(indicesForDisplay), meanTemperatureByMonth(indicesForDisplay),'LineWidth',plotLineWidth,'Color',temperatureLineColor,'LineStyle','-','DisplayName','Temperature (°C)');
    hold on
    %hm01 = bubblechart(groupTrafficTableMonths.monthStartDateTimes, meanTemperatureByMonth,totalPrecipitationByMonth,'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerFaceAlpha',0.3,'DisplayName',...
    hm01 = bubblechart(uniqueMonthsForDisplay(indicesForDisplay), meanTemperatureByMonth(indicesForDisplay),totalPrecipitationByMonth(indicesForDisplay),'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerFaceAlpha',0.3,'DisplayName',...
        'Marker Size Scaled by Total Monthly Precipitation');
    hold off
    ylabel('Mean Temperature for Month (°C)', 'FontSize', labelFontSize+2,'Color',uptimeLabelColor,'FontWeight','bold');
    ylim([-10 max(ylim).*1.1]);

    % Set up datatips for additional info
    hm01.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Location:',repmat({locationString}, 1, length(totalMonthlyTraffic)).');
    hm01.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Date:',groupTrafficTableMonths.monthStartDateTimes);
    hm01.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Raw Monthly Counts:',num2sepstr(totalMonthlyTraffic));
    hm01.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Min Monthly Counts:',num2sepstr(minMonthlyTraffic));
    hm01.DataTipTemplate.DataTipRows(5) = dataTipTextRow('Max Monthly Counts:',num2sepstr(maxMonthlyTraffic));
    hm01.DataTipTemplate.DataTipRows(6) = dataTipTextRow('Adjusted Counts:',num2sepstr(totalAdjustedMonthlyTraffic2,'%.0f'));
    hm01.DataTipTemplate.DataTipRows(7) = dataTipTextRow('Mean Temperature (°C):',round(meanTemperatureByMonth,1));
    hm01.DataTipTemplate.DataTipRows(8) = dataTipTextRow('Min Temperature (°C):',minTemperatureByMonth);
    hm01.DataTipTemplate.DataTipRows(9) = dataTipTextRow('Max Temperature (°C):',maxTemperatureByMonth);
    hm01.DataTipTemplate.DataTipRows(10) = dataTipTextRow('Total Precipitation (mm):',totalPrecipitationByMonth);
    hm01.DataTipTemplate.DataTipRows(11) = dataTipTextRow('Mean Windspeed (km/h):',meanWindspeedByMonth);
    hm01.DataTipTemplate.DataTipRows(12) = dataTipTextRow('Mean Feels Like Temp (°C):',meanFeelslikeByMonth);
    hm01.DataTipTemplate.DataTipRows(13) = dataTipTextRow('Total Snow (mm):',totalSnowByMonth);
    hm01.DataTipTemplate.DataTipRows(14) = dataTipTextRow('Sunrise:',sunriseByMonth);
    hm01.DataTipTemplate.DataTipRows(15) = dataTipTextRow('Sunset:',sunsetByMonth);
    hm01.DataTipTemplate.DataTipRows(16) = dataTipTextRow('Uptime (%):',round(uptimeByMonth,1).');

    yyaxis left

    hold on;

    % Plot raw traffic counts per month
    hm1 = plot(groupTrafficTableMonths.monthStartDateTimes, totalMonthlyTraffic,'-','LineWidth',plotLineWidth,'DisplayName',['Raw Counts for Month ( Total = ' num2sepstr(totalCountMonthly,'%.0f') ' )'],'Color','k');
    plot(groupTrafficTableMonths.monthStartDateTimes, totalMonthlyTraffic,'o','MarkerSize',plotLineWidth.*2,'Color','k','MarkerFaceColor','k');

    hm2 = plot(groupTrafficTableMonths.monthStartDateTimes, totalMonthlyTraffic2,':','LineWidth',plotLineWidth./3,'DisplayName',['Raw Counts for Month (up to ' datestr(truncationCutoffTime) ')'],'Color','k');
    plot(groupTrafficTableMonths.monthStartDateTimes, totalMonthlyTraffic2,'o','MarkerSize',plotLineWidth.*2,'Color','k','MarkerFaceColor','k');

    hm3 = plot(groupTrafficTableMonths.monthStartDateTimes, totalAdjustedMonthlyTraffic2,'-','LineWidth',plotLineWidth,'DisplayName','Estimated Counts for Month (adjusted for Uptime and Sunset)','Color',adjustedCountColor);
    scatter(groupTrafficTableMonths.monthStartDateTimes, totalAdjustedMonthlyTraffic2,plotLineWidth.*40,'o',...%'MarkerSize',plotLineWidth.*2,...
        'DisplayName','Counts for Monty (Adjusted for Uptime and Sunset)','Color',adjustedCountColor,...
        'MarkerFaceColor',adjustedCountColor(1:3),'MarkerFaceAlpha',0.3);

    % Add plot labels and formatting
    ylabel(['Total ' modeDisplayString ' for Month'], 'FontSize', labelFontSize+2, 'FontWeight', 'bold');
    xlabel('Month', 'FontSize', labelFontSize);

    ax2 = gca;
    ax2.FontSize = axisFontSize;
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = blueColor;

    title(['Monthly ' modeDisplayString ' on Terrebonne (' locationString ')'], 'FontSize', titleFontSize);
    subtitle(['Total estimated counts: ' num2sepstr(totalCountMonthly2,'%.0f') ' (Complete Months Only)']);

    xtickangle(45);
    grid on;

    % Legend and display settings
    legend([hm1 hm2 hm3 hm00 hm01], 'Location', 'southwest','Color',0.8.*[1 1 1]);
    bubblelegend('Weekly Precipitation (mm)','Color',0.8.*[1 1 1],'Location','southeast');
    set(gca, 'Color', 0.8 * [1 1 1]);
    ax = gca;
    ax.FontSize = axisFontSize;

    ylim([0 max(ylim).*1.1]);
    xticks(groupTrafficTableMonths.monthStartDateTimes);

    ytick_positions = yticks;
    ytick_labels = arrayfun(@(v) num2sepstr(v,'%.0f'), ytick_positions, 'UniformOutput', false);
    yticklabels(ytick_labels);

    hold off;

end

%% Plot total traffic counts per day over the full date range

plotSizeFactor = 0.5;

% Create daily dates for grouping (removing time component)
inputTable.DayOnly = dateshift(inputTable.Date, 'start', 'day');

% Calculate corrected daily totals accounting for uptime and daylight
dailyTableAdjusted = inputTable(timeofday(inputTable.Date)<=truncationCutoffTime,:);
groupAdjustedTrafficTableDays = groupsummary(dailyTableAdjusted, 'DayOnly', 'sum', 'AdjustedCountsUptimeDaylight');
totalAdjustedDailyTraffic = groupAdjustedTrafficTableDays.('sum_AdjustedCountsUptimeDaylight');
uniqueDays = groupAdjustedTrafficTableDays.DayOnly;

% Calculate daily totals from the input table, properly summing hourly counts
dailyTable = inputTable(ismember(inputTable.DayOnly,uniqueDays),:);
groupTrafficTableDays = groupsummary(dailyTable, 'DayOnly', 'sum', modeString);
%groupTrafficTableDays = groupTrafficTableDays(groupTrafficTableDays.GroupCount==24,:);
totalDailyTraffic = groupTrafficTableDays.(['sum_' modeString]);

totalAdjustedDailyTraffic = max(totalDailyTraffic, totalAdjustedDailyTraffic);

% Get weather data for each unique day
uniqueDays = groupTrafficTableDays.DayOnly; % Excludes partial first day
dailyNoonTimes = uniqueDays + hours(12);

% Retrieve weather data
disp('Getting detailed weather data for daily plot...')
[precipitationData, temperatureData, sunriseData, sunsetData, sunhoursData, snowData, windspeedData, feelslikeData] = getWeatherstackData('Montreal', dailyNoonTimes);
disp('Done')

% Calculate total counts
totalCountDaily = sum(totalDailyTraffic);
totalCountDaily2 = round(sum(totalAdjustedDailyTraffic), 0);
minCountDaily = min(totalDailyTraffic);
maxCountDaily = max(totalDailyTraffic);
minAdjustedCountDaily = round(min(totalAdjustedDailyTraffic),0);
maxAdjustedCountDaily = round(max(totalAdjustedDailyTraffic),0);

% Create the figure
figure('Position', [408 126 1132 921]);

blueColor = [0 0.4471 0.7412 0.3];
temperatureLineColor = blueColor;
uptimeLabelColor = [1 0 0 0.5] + 0.5.*[0 1 1 0];
adjustedCountColor = [0 0 0 0.3];

% Plot temperature data on right y-axis
yyaxis right
%hd0 = plot(uniqueDays, temperatureData, 'LineWidth', plotLineWidth.*plotSizeFactor, ...
hd0 = plot(uniqueDays, feelslikeData, 'LineWidth', plotLineWidth.*plotSizeFactor, ...
    'Color', temperatureLineColor, 'LineStyle', '-', 'DisplayName', '''Feels Like'' Temperature(°C)');
hold on
%hd1 = bubblechart(uniqueDays, temperatureData, precipitationData, ...
hd1 = bubblechart(uniqueDays, feelslikeData, precipitationData, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.3, ...
    'DisplayName', 'Marker Size Scaled by Daily Precipitation');
hold off

ylabel('Mean Daily Temperature (°C)', 'FontSize', labelFontSize+2, 'Color', uptimeLabelColor, 'FontWeight', 'bold');

% Set up datatips for the bubble chart
hd1.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Location:', repmat({locationString}, 1, length(totalDailyTraffic)).');
hd1.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Date:', uniqueDays);
hd1.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Raw Daily Counts:', num2sepstr(totalDailyTraffic));
hd1.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Adjusted Counts:', num2sepstr(totalAdjustedDailyTraffic, '%.0f'));
hd1.DataTipTemplate.DataTipRows(5) = dataTipTextRow('Temperature (°C):', round(temperatureData, 1));
hd1.DataTipTemplate.DataTipRows(6) = dataTipTextRow('Precipitation (mm):', precipitationData);
hd1.DataTipTemplate.DataTipRows(7) = dataTipTextRow('Windspeed (km/h):', windspeedData);
hd1.DataTipTemplate.DataTipRows(8) = dataTipTextRow('Feels Like Temp (°C):', feelslikeData);
hd1.DataTipTemplate.DataTipRows(9) = dataTipTextRow('Sunrise:', cellfun(@(x) datestr(x, 'HH:MM'), sunriseData, 'UniformOutput', false));
hd1.DataTipTemplate.DataTipRows(10) = dataTipTextRow('Sunset:', cellfun(@(x) datestr(x, 'HH:MM'), sunsetData, 'UniformOutput', false));
hd1.DataTipTemplate.DataTipRows(11) = dataTipTextRow('Sun Hours:', sunhoursData);

% Plot traffic counts on left y-axis
yyaxis left

% Plot raw counts
hd2 = plot(uniqueDays, totalDailyTraffic, '-', 'LineWidth', plotLineWidth.*plotSizeFactor, ...
    'DisplayName', ['Raw Counts per Day ( Min = ' num2sepstr(minCountDaily, '%.0f') ' ; Max = ' num2sepstr(maxCountDaily,'%.0f') ' )'], ...
    'Color', 'k');
hold on

% Plot adjusted counts [NOT ROBUST FOR DAILY COUNTS]
% hd3 = plot(uniqueDays, totalAdjustedDailyTraffic, '-', 'LineWidth', plotLineWidth.*plotSizeFactor, ...
%      'DisplayName', 'Estimated Counts per Day (adjusted for Uptime and Sunset)', ...
%      'Color', adjustedCountColor);
hd3 = plot(uniqueDays, totalAdjustedDailyTraffic, '-', 'LineWidth', plotLineWidth.*plotSizeFactor, ...
     'DisplayName', ['Adjusted Counts per Day ( Min = ' num2sepstr(minAdjustedCountDaily, '%.0f') ' ; Max = ' num2sepstr(maxAdjustedCountDaily,'%.0f') ' )'], ...
     'Color', adjustedCountColor);

%hd4 = plot([uniqueDays(1) uniqueDays(end)],[5 5],'k:','LineWidth',1.5,'DisplayName','Five Counts Per Day Reference Line');

% Add scatter points for both raw and adjusted data
% scatter(uniqueDays, totalAdjustedDailyTraffic, plotLineWidth*40.*plotSizeFactor, 'o', ...
%     'DisplayName', 'Adjusted Daily Counts', 'Color', adjustedCountColor, ...
%     'MarkerFaceColor', adjustedCountColor(1:3), 'MarkerFaceAlpha', 0.3);

hold off

% Format the plot
ylabel(['Total ' modeDisplayString ' per Day'], 'FontSize', labelFontSize+2, 'FontWeight', 'bold');
xlabel('Date', 'FontSize', labelFontSize);

ax = gca;
ax.FontSize = axisFontSize;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = blueColor;

title(['Daily ' modeDisplayString ' on Terrebonne (' locationString ')'], 'FontSize', titleFontSize);
subtitle(['Total estimated counts: ' num2sepstr(totalCountDaily2, '%.0f') ' (' ...
    datestr(min(uniqueDays), 'yyyy-mm-dd') ' to ' ...
    datestr(max(uniqueDays), 'yyyy-mm-dd') ')']);

% Add legend and set display properties
lh1 = legend([hd2 hd3 hd0 hd1], 'Location', 'southwest', 'Color', 0.8.*[1 1 1],'BackgroundAlpha',legendBackgroundAlpha);
%legend([hd2 hd0 hd1], 'Location', 'southwest', 'Color', 0.8.*[1 1 1]);
lh2 = bubblelegend('Daily Precipitation (mm)', 'Color', 0.8.*[1 1 1], 'Location', 'southeast');

% Set legend positions - modify as needed
set(lh1,'Position',[0.5177 0.8160 0.3768 0.1015]);
set(lh2,'Position',[0.7452 0.6487 0.1484 0.1546]);


set(gca, 'Color', 0.8 * [1 1 1]);
grid on;
xtickangle(45);

% Format y-axis tick labels with separators
ytick_positions = yticks;
ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
yticklabels(ytick_labels);

% Adjust y-axis limits to add some padding
ylim([0 max(ylim).*1.1]);
xlim([uniqueDays(1) uniqueDays(end)])

% Save values for aggregate plot
saveName = [ strrep(locationString,' ','_') '_Daily_Counts_Saved.mat' ];
save(saveName, 'uniqueDays', 'totalDailyTraffic','totalAdjustedDailyTraffic');

saveName2 = [ strrep(locationString,' ','_') '_Daily_Weather_Saved.mat' ];
save(saveName2, 'uniqueDays', 'temperatureData', 'precipitationData','windspeedData','feelslikeData');

% Create a table combining all the data
combinedData = table(uniqueDays, totalDailyTraffic, totalAdjustedDailyTraffic, ...
                     temperatureData, precipitationData, windspeedData, feelslikeData, ...
                     'VariableNames', {'Date', 'TotalTraffic', 'AdjustedTraffic', ...
                                       'Temperature', 'Precipitation', 'WindSpeed', 'FeelsLike'});

% Create CSV filename based on the location
csvFileName = [strrep(locationString,' ','_') '_Count_and_Weather_Data.csv'];

% Write the table to a CSV file
writetable(combinedData, csvFileName);

% Display confirmation message
fprintf('Combined data successfully saved to %s\n', csvFileName);

return;

%% Plot total traffic counts per day over the full date range AT TWO LOCATIONS

clc

plotLineWidth = 10.0;
axisFontSize = 16.0;
labelFontSize = 20.0;
titleFontSize = 24.0;
legendFontSize = 16.0;
axisBackgroundColor = 0.8.*[1 1 1];

plotSizeFactor = 0.35;

% % Create daily dates for grouping (removing time component)
% inputTable.DayOnly = dateshift(inputTable.Date, 'start', 'day');
% 
% % Calculate daily totals from the input table, properly summing hourly counts
% dailyTable = inputTable;
% groupTrafficTableDays = groupsummary(dailyTable, 'DayOnly', 'sum', modeString);
% totalDailyTraffic = groupTrafficTableDays.(['sum_' modeString]);
% 
% % Calculate corrected daily totals accounting for uptime and daylight
% groupAdjustedTrafficTableDays = groupsummary(dailyTable, 'DayOnly', 'sum', 'AdjustedCountsUptimeDaylight');
% totalAdjustedDailyTraffic = groupAdjustedTrafficTableDays.('sum_AdjustedCountsUptimeDaylight');
% 
% % Ensure adjusted values are at least as large as raw counts
% totalAdjustedDailyTraffic = max(totalDailyTraffic, totalAdjustedDailyTraffic);
% 
% % Get weather data for each unique day
% uniqueDays = groupTrafficTableDays.DayOnly;
% noonTimes = uniqueDays + hours(12);
% 
% % Retrieve weather data
% disp('Getting detailed weather data for daily plot...')
% %[precipitationData, temperatureData, sunriseData, sunsetData, sunhoursData] = getWeatherstackData('Montreal', noonTimes);
% disp('Done')

% % Calculate total counts
% totalCountDaily = sum(totalDailyTraffic);
% totalCountDaily2 = round(sum(totalAdjustedDailyTraffic), 0);
% minCountDaily = min(totalDailyTraffic);
% maxCountDaily = max(totalDailyTraffic);

% Create the figure
figure('Position', [408 126 1132 921]);

blueColor = [0 0.4471 0.7412 0.3];
temperatureLineColor = blueColor;
uptimeLabelColor = [1 0 0 0.5] + 0.5.*[0 1 1 0];
adjustedCountColor = [0 0 0 0.3];

% Plot temperature data on right y-axis

load 'Western_Segment_Daily_Weather_Saved.mat'

yyaxis right
hd0 = plot(uniqueDays, temperatureData, 'LineWidth', plotLineWidth.*plotSizeFactor, ...
    'Color', temperatureLineColor, 'LineStyle', '-', 'DisplayName', 'Temperature (°C)');
hold on
hd1 = bubblechart(uniqueDays, temperatureData, precipitationData, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.3, ...
    'DisplayName', 'Marker Size Scaled by Daily Precipitation');
hold off

ylabel('Mean Daily Temperature (°C)', 'FontSize', labelFontSize+2, 'Color', uptimeLabelColor, 'FontWeight', 'bold');

% Set up datatips for the bubble chart
% hd1.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Location:', repmat({locationString}, 1, length(totalDailyTraffic)).');
% hd1.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Date:', uniqueDays);
% hd1.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Raw Daily Counts:', num2sepstr(totalDailyTraffic));
% hd1.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Adjusted Counts:', num2sepstr(totalAdjustedDailyTraffic, '%.0f'));
% hd1.DataTipTemplate.DataTipRows(5) = dataTipTextRow('Temperature (°C):', round(temperatureData, 1));
% hd1.DataTipTemplate.DataTipRows(6) = dataTipTextRow('Precipitation (mm):', precipitationData);
% hd1.DataTipTemplate.DataTipRows(7) = dataTipTextRow('Sunrise:', cellfun(@(x) datestr(x, 'HH:MM'), sunriseData, 'UniformOutput', false));
% hd1.DataTipTemplate.DataTipRows(8) = dataTipTextRow('Sunset:', cellfun(@(x) datestr(x, 'HH:MM'), sunsetData, 'UniformOutput', false));
% hd1.DataTipTemplate.DataTipRows(9) = dataTipTextRow('Sun Hours:', sunhoursData);

% Plot traffic counts on left y-axis
yyaxis left

% Plot raw counts

useRawCounts = true;

clear uniqueDays totalDailyTraffic totalAdjustedDailyTraffic;
load('Eastern_Segment_Daily_Counts_Saved.mat');

if useRawCounts
    totalDailyTrafficForDisplay = totalDailyTraffic;
    countString = 'Raw';
else
    totalDailyTrafficForDisplay = totalAdjustedDailyTraffic;
    countString = 'Adjusted';
end

% Calculate total counts
totalCountDaily = sum(totalDailyTrafficForDisplay);
minCountDaily = min(totalDailyTrafficForDisplay);
maxCountDaily = max(totalDailyTrafficForDisplay);

hd2a = plot(uniqueDays, totalDailyTrafficForDisplay, '-', 'LineWidth', plotLineWidth.*plotSizeFactor, ...
    'DisplayName', ['Counts per Day East Segment ( Min = ' num2sepstr(minCountDaily, '%.0f') ' ; Max = ' num2sepstr(maxCountDaily,'%.0f') ' ; Total = ' num2sepstr(totalCountDaily,'%.0f') ' )'], ...
    'Color', 'b');
hold on

clear uniqueDays totalDailyTraffic totalAdjustedDailyTraffic;
load('Western_Segment_Daily_Counts_Saved.mat');

% Calculate total counts

if useRawCounts
    totalDailyTrafficForDisplay = totalDailyTraffic;
    countString = 'Raw';
else
    totalDailyTrafficForDisplay = totalAdjustedDailyTraffic;
    countString = 'Adjusted';
end
totalCountDaily2 = sum(totalDailyTrafficForDisplay);
minCountDaily2 = min(totalDailyTrafficForDisplay);
maxCountDaily2 = max(totalDailyTrafficForDisplay);

hd2b = plot(uniqueDays, totalDailyTrafficForDisplay, '-', 'LineWidth', plotLineWidth.*plotSizeFactor, ...
     'DisplayName', ...
     ['Counts per Day West Segment ( Min = ' num2sepstr(minCountDaily2, '%.0f') ' ; Max = ' num2sepstr(maxCountDaily2,'%.0f') ' ; Total = ' num2sepstr(totalCountDaily2,'%.0f') ' )'], ...
     'Color', 'k');

% Plot adjusted counts [NOT ROBUST FOR DAILY COUNTS]
% hd3 = plot(uniqueDays, totalAdjustedDailyTraffic, '-', 'LineWidth', plotLineWidth.*plotSizeFactor, ...
%     'DisplayName', 'Estimated Counts per Day (adjusted for Uptime and Sunset)', ...
%     'Color', adjustedCountColor);

%hd3 = plot([uniqueDays(1) uniqueDays(end)],[5 5],'k:','LineWidth',1.5,'DisplayName','Five Counts Per Day Reference Line');

% Add scatter points for both raw and adjusted data
% scatter(uniqueDays, totalAdjustedDailyTraffic, plotLineWidth*40.*plotSizeFactor, 'o', ...
%     'DisplayName', 'Adjusted Daily Counts', 'Color', adjustedCountColor, ...
%     'MarkerFaceColor', adjustedCountColor(1:3), 'MarkerFaceAlpha', 0.3);

hold off

% Format the plot
ylabel(['Total ' modeDisplayString ' per Day'], 'FontSize', labelFontSize+2, 'FontWeight', 'bold');
xlabel('Date', 'FontSize', labelFontSize);

ax = gca;
ax.FontSize = axisFontSize;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = blueColor;

locationString2 = 'Eastern and Western Segments';
title(['Daily ' modeDisplayString ' on Terrebonne (' locationString2 ')'], 'FontSize', titleFontSize);
subtitle(['Total ' countString ' counts: ' num2sepstr(totalCountDaily + totalCountDaily2, '%.0f') ' (' ...
    datestr(min(uniqueDays), 'yyyy-mm-dd') ' to ' ...
    datestr(max(uniqueDays), 'yyyy-mm-dd') ')']);

% Add legend and set display properties
lh1 = legend([hd2a hd2b hd0 hd1], 'Location', 'southwest', 'Color', 0.8.*[1 1 1],'BackgroundAlpha',legendBackgroundAlpha);
%legend([hd2 hd0 hd1], 'Location', 'southwest', 'Color', 0.8.*[1 1 1]);
lh2 = bubblelegend('Daily Precipitation (mm)', 'Color', 0.8.*[1 1 1], 'Location', 'southeast');

% Set legend positions - modify as needed
set(lh1,'Position',[0.4271    0.8138    0.4678    0.1015]);
set(lh2,'Position',[0.7452    0.6477    0.1484    0.1546]);


set(gca, 'Color', 0.8 * [1 1 1]);
grid on;
xtickangle(45);

% Format y-axis tick labels with separators
ytick_positions = yticks;
ytick_labels = arrayfun(@(v) num2sepstr(v, '%.0f'), ytick_positions, 'UniformOutput', false);
yticklabels(ytick_labels);

% Adjust y-axis limits to add some padding
ylim([0 max(ylim).*1.3]);
xlim([uniqueDays(2) uniqueDays(end)])

% Save values for aggregate plot
saveName = [ strrep(locationString,' ','_') '_Daily_Counts_Saved.mat' ];
save(saveName, 'uniqueDays', 'totalDailyTrafficForDisplay');

% Aggregate the values into a single table



%% FUNCTION to compute Monday for a given year/week combination
function mondayDate = weekNumToMonday(year, weekNum)
    % Input validation
    if length(year) ~= length(weekNum)
        error('Input vectors year and weekNum must be the same length');
    end
    
    % Preallocate output array
    mondayDate = NaT(size(year));
    
    % Process each year/week combination
    for i = 1:length(year)
        % Create January 4th at noon
        jan4 = datetime(year(i), 1, 4, 12, 0, 0);
        
        % Get the day of week (1 = Sunday, 2 = Monday, etc.)
        dayNum = weekday(jan4);
        
        % Calculate days to subtract to get to Monday of week 1
        if dayNum == 1  % Sunday
            daysToSubtract = 6;
        else
            daysToSubtract = dayNum - 2;
        end
        
        % Get to week 1 Monday and add required weeks
        mon1 = jan4 - days(daysToSubtract);
        mondayDate(i) = mon1 + days(7*(weekNum(i)-1));
        
        % Explicitly set the time to noon
        mondayDate(i).Hour = 12;
        mondayDate(i).Minute = 0;
        mondayDate(i).Second = 0;
    end
end

%% FUNCTION to convert num to string with thousands separator

function out = num2sepstr(numin, format, sep)
% NUM2SEPSTR Convert to string with separation at thousands.
%
% out = NUM2SEPSTR(numin,[format],[sep]) formats numin to a string
%   according to the specified format ('%f' by default) and adds the
%   sepcified thousands seperators (commas by default).
%
% For non-scalar numin, num2sepstr outpts a cell array of the same shape
%   as numin where num2sepstr is called on each value in numin.
%
% String length from format, when specified, is applied before commas are
%   added; Instead of...
%
%   >> num2sepstr(1e6,'% 20.2f') % length = 22
%   ans =
%       '          1,000,000.00'
%
% ...try...
%
%   >> sprintf('% 20s',num2sepstr(1e6,'%.2f')) % length = 20
%   ans =
%       '        1,000,000.00'
%
% See also SPRINTF, NUM2STR
%
% Created by:
%   Robert Perrotta
% Thanks to MathWorks community members Stephen Cobeldick and Andreas J.
% for suggesting a faster, cleaner implementation using regexp and
% regexprep.
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


%% FUNCTION to fetch weather data from Weatherstack
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

    %disp(['Getting weather data for date ' string(currentDate)]);

    % Convert the current date to string format
    dateStr = datestr(currentDate, 'yyyy-mm-dd');

    % Construct the API request URL
    requestURL = sprintf('%s?access_key=%s&query=%s&historical_date=%s&hourly=1&interval=24&units=m', baseURL, accessKey, location, dateStr);

    % Submit the request
    response = webread(requestURL);

    % Extract precipitation data
    %if isfield(response, 'historical') && isfield(response.historical, dateStr)
    precipitation = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.precip;
    windspeed = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.wind_speed;
    feelslike = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.feelslike;
    temperature = response.historical.(['x' strrep(dateStr,'-','_')]).hourly.temperature;
    avgtemp = response.historical.(['x' strrep(dateStr,'-','_')]).avgtemp;
    totalsnow = response.historical.(['x' strrep(dateStr,'-','_')]).totalsnow;
    sunhours = response.historical.(['x' strrep(dateStr,'-','_')]).sunhour;
    sunriseTimeStr = response.historical.(['x' strrep(dateStr,'-','_')]).astro.sunrise;
    sunsetTimeStr = response.historical.(['x' strrep(dateStr,'-','_')]).astro.sunset;
    %else
    %    precipitation = NaN; % Handle missing data gracefully
    %end

    precipitationValues(ix) = precipitation;
    windspeedValues(ix) = windspeed;
    feelslikeValues(ix) = feelslike;
    %temperatureValues(ix) = temperature;
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

function aggregateCurve = estimateMissingTraffic(canonicalTrafficModel, avgTraffic, daylightIndices, darkIndices)
    % estimateMissingTraffic
    % Estimates traffic counts for missing times using a canonical model.
    %
    % INPUTS:
    % canonicalTrafficModel: A vector containing the canonical traffic model
    %                        at hourly intervals over 24 hours.
    % avgTraffic: A vector containing the average hourly traffic counts.
    % daylightIndices: Indices to consider for fitting.
    % darkIndices: Indices to exclude from fitting (potentially missing data).
    %
    % OUTPUT:
    % aggregateCurve: A vector containing the estimated traffic counts for all
    %                 hours, including interpolations for darkIndices based on
    %                 the canonical traffic model.

    % Ensure inputs are column vectors for consistency
    canonicalTrafficModel = canonicalTrafficModel(:);
    avgTraffic = avgTraffic(:);

    % Extract data for fitting based on daylightIndices
    canonicalDaylight = canonicalTrafficModel(daylightIndices);
    avgTrafficDaylight = avgTraffic(daylightIndices);

    % Fit the canonical model to the observed daylight traffic using least squares
    scaleFactor = (canonicalDaylight' * avgTrafficDaylight) / ...
                  (canonicalDaylight' * canonicalDaylight);

    % Apply the scale factor to the canonical model to estimate traffic
    scaledCanonicalTraffic = scaleFactor * canonicalTrafficModel;

    % Replace the dark indices with estimates from the scaled canonical model
    aggregateCurve = avgTraffic; % Start with the original traffic data
    aggregateCurve(darkIndices) = scaledCanonicalTraffic(darkIndices); % Fill missing values
end

% Function to calculate min/max daily counts per time period (weekly or monthly)
function [minDailyCounts, maxDailyCounts] = calculatePeriodDailyMinMax(inputTable, modeString, periodStartDateTimes, periodColumnName)
    % inputTable: timetable with hourly data
    % modeString: column name containing the count data
    % periodStartDateTimes: vector of period start datetimes (weeks or months)
    % periodColumnName: name of the column containing period start dates (default: 'weekStartDateTimes')
    %
    % Returns:
    % minDailyCounts: vector of minimum daily counts for each period
    % maxDailyCounts: vector of maximum daily counts for each period
    
    % Set default period column name if not provided
    if nargin < 4 || isempty(periodColumnName)
        periodColumnName = 'weekStartDateTimes';
    end
    
    % Filter the input table by the provided period start dates
    if nargin >= 3 && ~isempty(periodStartDateTimes)
        filteredTable = inputTable(ismember(inputTable.(periodColumnName), periodStartDateTimes), :);
    else
        filteredTable = inputTable;
        % If no specific periods provided, use all unique periods in the data
        periodStartDateTimes = unique(inputTable.(periodColumnName));
    end
    
    % Ensure we have data after filtering
    if isempty(filteredTable)
        warning('No data found for the specified %s', periodColumnName);
        minDailyCounts = NaN(size(periodStartDateTimes));
        maxDailyCounts = NaN(size(periodStartDateTimes));
        return;
    end
    
    % Prepare output arrays matching the input periodStartDateTimes
    minDailyCounts = NaN(length(periodStartDateTimes), 1);
    maxDailyCounts = NaN(length(periodStartDateTimes), 1);
    
    % Process each period
    for i = 1:length(periodStartDateTimes)
        periodStart = periodStartDateTimes(i);
        
        % Get data for this period
        periodData = filteredTable(filteredTable.(periodColumnName) == periodStart, :);
        
        if isempty(periodData)
            % No data for this period, leave as NaN
            continue;
        end
        
        % Get day dates for this period's data
        dayDates = dateshift(periodData.Properties.RowTimes, 'start', 'day');
        
        % Get unique days
        uniqueDays = unique(dayDates);
        
        % Calculate daily totals for each unique day
        dailyTotals = zeros(length(uniqueDays), 1);
        
        for j = 1:length(uniqueDays)
            dayData = periodData(dayDates == uniqueDays(j), :);
            dailyTotals(j) = sum(dayData.(modeString));
        end
        
        % Calculate min and max
        if ~isempty(dailyTotals)
            minDailyCounts(i) = min(dailyTotals);
            maxDailyCounts(i) = max(dailyTotals);
        end
    end
end

% Convenience functions for weekly and monthly periods
function [minDailyCounts, maxDailyCounts] = calculateWeeklyDailyMinMax(inputTable, modeString, weekStartDateTimes)
    [minDailyCounts, maxDailyCounts] = calculatePeriodDailyMinMax(inputTable, modeString, weekStartDateTimes, 'weekStartDateTimes');
end

function [minDailyCounts, maxDailyCounts] = calculateMonthlyDailyMinMax(inputTable, modeString, monthStartDateTimes)
    [minDailyCounts, maxDailyCounts] = calculatePeriodDailyMinMax(inputTable, modeString, monthStartDateTimes, 'monthStartDateTimes');
end

% Backward compatibility function
function maxDailyCounts = calculateWeeklyMaxDailyCounts(inputTable, modeString, weekStartDateTimes)
    [~, maxDailyCounts] = calculateWeeklyDailyMinMax(inputTable, modeString, weekStartDateTimes);
end