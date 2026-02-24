% UNformat listserv files for table - CFI site owner/contact categories
% (PIs/collaborators/other)
% Rosie Howard
% 6 January 2026


%% Settings

addpath(genpath('/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux'));

filepath = '/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux';
inputFile  = fullfile(filepath, 'CFI_siteOwnerCategories_contactList_fromSitesSpreadsheetComparison.txt');
% inputFile  = fullfile(filepath, 'CanFlux contacts(CanFlux listserv).csv');
outputFile = fullfile(filepath, 'CanFlux_PI_contact_categories.csv');

% Desired column order (example: reverse column order)
% Change this vector as needed, e.g. [3 1 2 4]
% columnOrder = 'reverse';  % options: 'reverse' or numeric vector

%% Read CSV
T = readtable(inputFile, 'TextType', 'string','Delimiter',' ');

%% Reorder columns
% if ischar(columnOrder) || isstring(columnOrder)
%     if strcmp(columnOrder, 'reverse')
%         T = T(:, end:-1:1);
%     else
%         error('Unknown columnOrder option.');
%     end
% else
%     T = T(:, columnOrder);
% end

%% Convert all columns to strings
T = varfun(@string, T);

%% Merge columns with single spaces
mergedColumn = join(T.Variables, " ", 2);

%% Create output table
Tout = table(mergedColumn, 'VariableNames', {'MergedText'});

%% Write to CSV
writetable(Tout, outputFile);

disp('Done: columns reordered and merged with spaces.');


% PIs listserv
%% parameters

filepath = '/Users/rosie/Documents/Micromet/CANFLUX_Database/ListServ';
inputFile  = fullfile(filepath, 'CanFluxPIs.xlsx');
outputFile = fullfile(filepath, 'CanFluxPIs_unique.csv');

% Read Excel file
raw = readcell(inputFile);

% Flatten to one column
oneCol = raw(:);

% Remove empty cells
oneCol = oneCol(~cellfun(@isempty, oneCol));


% Convert to string for consistency
oneCol = string(oneCol);

% Strip leading whitespace only
oneCol = strip(oneCol, 'left');

% Initialize output columns
col1 = strings(size(oneCol));
col2 = strings(size(oneCol));

% Identify entries containing "@"
% hasAt = contains(oneCol, "@");
delimiter = "(";
output1 = [];
output2 = [];
for i = 1:length(oneCol)
    [newStr, match] = split(oneCol(i), delimiter);

    [a,b] = size(newStr);
    name = newStr(1);
    if a == 2
        email = newStr(2);
        email = erase(email,')');
    else
        email = 'no email yet';
    end

    output1 = [output1; name];
    output2 = [output2; email];

end
output = [output1, output2];

[~, unique_idx] = unique(output(:, 1), 'stable'); % 'stable' preserves original order

% Keep only the rows corresponding to the unique indices
output_unique = output(unique_idx, :);

% Create table
T = table(output_unique(:,1), output_unique(:,2), ...
    'VariableNames', {'Name','Email'});
writetable(T,outputFile);