%% Inspect Column Structure of Telraam Data Files
% This script loads both 60-minute and 15-minute resolution data files
% and displays their column structure for comparison

clear all; clc;

fprintf('================================================\n');
fprintf('Telraam Data Column Structure Inspector\n');
fprintf('================================================\n\n');

%% Define file pairs to check
filePairs = {
    % {60-min file, 15-min file, description}
    {'raw-data-2024East60.xlsx', 'raw-data-2024East15.xlsx', '2024 Eastern'};
    {'raw-data-2025East60.xlsx', 'raw-data-2025East15.xlsx', '2025 Eastern'};
    {'raw-data-2024West60.xlsx', 'raw-data-2024West15.xlsx', '2024 Western'};
    {'raw-data-2025West60.xlsx', 'raw-data-2025West15.xlsx', '2025 Western'};
};

%% Process each file pair
for pairIdx = 1:length(filePairs)
    file60 = filePairs{pairIdx}{1};
    file15 = filePairs{pairIdx}{2};
    description = filePairs{pairIdx}{3};
    
    fprintf('========================================\n');
    fprintf('%s Data Comparison\n', description);
    fprintf('========================================\n\n');
    
    % Load 60-minute data
    if exist(file60, 'file')
        fprintf('Loading 60-minute data: %s\n', file60);
        data60 = readtable(file60);
        cols60 = data60.Properties.VariableNames;
        fprintf('  -> Loaded successfully: %d rows, %d columns\n\n', height(data60), width(data60));
    else
        fprintf('60-minute file not found: %s\n\n', file60);
        cols60 = {};
    end
    
    % Load 15-minute data
    if exist(file15, 'file')
        fprintf('Loading 15-minute data: %s\n', file15);
        data15 = readtable(file15);
        cols15 = data15.Properties.VariableNames;
        fprintf('  -> Loaded successfully: %d rows, %d columns\n\n', height(data15), width(data15));
    else
        fprintf('15-minute file not found: %s\n\n', file15);
        cols15 = {};
    end
    
    % Compare column structures if both files exist
    if ~isempty(cols60) && ~isempty(cols15)
        fprintf('Column Comparison:\n');
        fprintf('------------------\n');
        
        % Check if columns are identical
        if isequal(cols60, cols15)
            fprintf('✓ Column structures are IDENTICAL\n\n');
        else
            fprintf('✗ Column structures DIFFER\n\n');
            
            % Find differences
            only60 = setdiff(cols60, cols15);
            only15 = setdiff(cols15, cols60);
            
            if ~isempty(only60)
                fprintf('Columns only in 60-min data:\n');
                for i = 1:length(only60)
                    fprintf('  - %s\n', only60{i});
                end
                fprintf('\n');
            end
            
            if ~isempty(only15)
                fprintf('Columns only in 15-min data:\n');
                for i = 1:length(only15)
                    fprintf('  - %s\n', only15{i});
                end
                fprintf('\n');
            end
        end
        
        % Display all columns for first pair only (for reference)
        if pairIdx == 1
            fprintf('Full Column List (from %s):\n', description);
            fprintf('--------------------------------------------\n');
            for i = 1:length(cols60)
                fprintf('[%2d] %s', i, cols60{i});
                
                % Check for mode-related columns
                if contains(lower(cols60{i}), 'bike')
                    fprintf(' <- BIKE');
                elseif contains(lower(cols60{i}), 'car')
                    fprintf(' <- CAR');
                elseif contains(lower(cols60{i}), 'pedestrian')
                    fprintf(' <- PEDESTRIAN');
                elseif contains(lower(cols60{i}), 'date') || contains(lower(cols60{i}), 'time')
                    fprintf(' <- DATE/TIME');
                end
                fprintf('\n');
            end
            fprintf('\n');
            
            % Show first few rows of data to understand structure
            if height(data60) > 0
                fprintf('Sample data (first 3 rows of 60-min):\n');
                fprintf('--------------------------------------\n');
                disp(data60(1:min(3,height(data60)), :));
            end
            
            if height(data15) > 0
                fprintf('Sample data (first 3 rows of 15-min):\n');
                fprintf('--------------------------------------\n');
                disp(data15(1:min(3,height(data15)), :));
            end
        end
    end
    
    fprintf('\n');
end

%% Summary
fprintf('================================================\n');
fprintf('Summary\n');
fprintf('================================================\n\n');

% Try to identify the exact column names for key fields
if exist('data60', 'var') && ~isempty(data60)
    fprintf('Detected Column Names (from last loaded file):\n\n');
    
    % Find date/time column
    dateColumns = cols60(contains(lower(cols60), 'date') | contains(lower(cols60), 'time'));
    if ~isempty(dateColumns)
        fprintf('Date/Time columns:\n');
        for i = 1:length(dateColumns)
            fprintf('  "%s"\n', dateColumns{i});
        end
        fprintf('\n');
    end
    
    % Find bike columns
    bikeColumns = cols60(contains(lower(cols60), 'bike'));
    if ~isempty(bikeColumns)
        fprintf('Bike-related columns:\n');
        for i = 1:length(bikeColumns)
            fprintf('  "%s"\n', bikeColumns{i});
        end
        fprintf('\n');
    end
    
    % Find car columns  
    carColumns = cols60(contains(lower(cols60), 'car'));
    if ~isempty(carColumns)
        fprintf('Car-related columns:\n');
        for i = 1:length(carColumns)
            fprintf('  "%s"\n', carColumns{i});
        end
        fprintf('\n');
    end
    
    % Find total columns
    totalColumns = cols60(contains(lower(cols60), 'total'));
    if ~isempty(totalColumns)
        fprintf('Total count columns:\n');
        for i = 1:length(totalColumns)
            fprintf('  "%s"\n', totalColumns{i});
        end
        fprintf('\n');
    end
    
    % Check data type of first column (usually datetime)
    if ~isempty(dateColumns)
        firstDateCol = dateColumns{1};
        fprintf('Data type of "%s": %s\n', firstDateCol, class(data60.(firstDateCol)));
        if isdatetime(data60.(firstDateCol))
            fprintf('  -> First value: %s\n', datestr(data60.(firstDateCol)(1)));
            fprintf('  -> Last value: %s\n', datestr(data60.(firstDateCol)(end)));
        end
    end
end

fprintf('\n================================================\n');
fprintf('Inspection Complete\n');
fprintf('================================================\n');