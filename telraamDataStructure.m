%% Diagnostic Script for Data Structure
% Run this after runTelraamAnalysis.m to understand the data structure

fprintf('\n========== DATA STRUCTURE DIAGNOSTIC ==========\n\n');

%% Check workspace variables
fprintf('--- Workspace Variables ---\n');
if exist('locationData', 'var')
    fprintf('✓ locationData exists\n');
    locationNames = fieldnames(locationData);
    fprintf('  Locations: %s\n', strjoin(locationNames, ', '));
else
    fprintf('✗ locationData NOT FOUND\n');
end

if exist('weatherData', 'var')
    fprintf('✓ weatherData exists\n');
else
    fprintf('✗ weatherData NOT FOUND\n');
end

if exist('analysis', 'var')
    fprintf('✓ analysis exists\n');
    fprintf('  Mode: %s\n', analysis.modeString);
else
    fprintf('✗ analysis NOT FOUND\n');
end

fprintf('\n');

%% Examine locationData structure
if exist('locationData', 'var')
    fprintf('--- Location Data Structure ---\n');
    locationNames = fieldnames(locationData);
    
    for i = 1:length(locationNames)
        locationName = locationNames{i};
        fprintf('\nLocation: %s\n', locationName);
        
        % Check what fields exist
        locFields = fieldnames(locationData.(locationName));
        fprintf('  Fields: %s\n', strjoin(locFields, ', '));
        
        % Check the data field
        if isfield(locationData.(locationName), 'data')
            data = locationData.(locationName).data;
            fprintf('  Data class: %s\n', class(data));
            
            if istimetable(data)
                fprintf('  Timetable with %d rows\n', height(data));
                fprintf('  Variables: %s\n', strjoin(data.Properties.VariableNames, ', '));
                fprintf('  Time range: %s to %s\n', ...
                    datestr(min(data.Properties.RowTimes)), ...
                    datestr(max(data.Properties.RowTimes)));
                
                % Check if mode exists
                if exist('analysis', 'var') && isfield(analysis, 'modeString')
                    if any(strcmp(data.Properties.VariableNames, analysis.modeString))
                        fprintf('  ✓ Mode "%s" found in data\n', analysis.modeString);
                        modeData = data.(analysis.modeString);
                        fprintf('    Min: %.0f, Max: %.0f, Mean: %.1f\n', ...
                            min(modeData), max(modeData), mean(modeData));
                    else
                        fprintf('  ✗ Mode "%s" NOT found in data\n', analysis.modeString);
                    end
                end
                
            elseif istable(data)
                fprintf('  Table with %d rows\n', height(data));
                fprintf('  Variables: %s\n', strjoin(data.Properties.VariableNames, ', '));
                
                % Check for date column
                if any(strcmp(data.Properties.VariableNames, 'Date and Time (Local)'))
                    dateCol = data.('Date and Time (Local)');
                    fprintf('  Date range: %s to %s\n', ...
                        datestr(min(dateCol)), datestr(max(dateCol)));
                end
                
            else
                fprintf('  Data is neither table nor timetable\n');
            end
        end
    end
end

%% Examine weatherData structure
if exist('weatherData', 'var')
    fprintf('\n--- Weather Data Structure ---\n');
    weatherFields = fieldnames(weatherData);
    fprintf('Fields: %s\n', strjoin(weatherFields, ', '));
    
    if isfield(weatherData, 'dates')
        fprintf('Number of dates: %d\n', length(weatherData.dates));
        fprintf('Date range: %s to %s\n', ...
            datestr(min(weatherData.dates)), datestr(max(weatherData.dates)));
    end
    
    if isfield(weatherData, 'temperature')
        fprintf('Temperature: min=%.1f, max=%.1f, mean=%.1f\n', ...
            min(weatherData.temperature), max(weatherData.temperature), ...
            mean(weatherData.temperature));
    end
    
    if isfield(weatherData, 'precipitation')
        fprintf('Precipitation: min=%.1f, max=%.1f, mean=%.1f\n', ...
            min(weatherData.precipitation), max(weatherData.precipitation), ...
            mean(weatherData.precipitation));
    end
end

fprintf('\n========== DIAGNOSTIC COMPLETE ==========\n\n');