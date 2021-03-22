function exportGroupedData(title, columnData, groupTitles, filename,...
    sheetNum)
%exportGroupedData writes data to an Excel file that is formatted for
%easy importing to GraphPad
%   Input arguments:
%       title: column heading of the data to be exported
%       columnData: single column of data, **SORTED** by group
%       replicates: number of repeat measures per group
%       filename: full path to the file to be exported to
%       sheetNum: sheet that the data should be written to
%   Output arguments:
%       none

numColumns = length(columnData);
numGroups = length(groupTitles);
% Check that the input data will fit in the specified output size.
if rem(numColumns, numGroups) > 0
    error(['columnData with length ' num2str(numColumns)...
        ' is not evenly divisible by group count '...
        num2str(numGroups) '!']);
end
replicates = numColumns / numGroups;

% Reshape the input data.
outputArray = reshape(columnData, [replicates, numGroups])';

% Export the table.
outputHeaders = cell(1, replicates);
for i = 1:replicates
    outputHeaders{i} = strcat(title, num2str(i));
end
writetable(cell2table(outputArray, 'VariableNames', outputHeaders,...
    'RowNames', groupTitles), filename, 'Sheet', sheetNum,...
    'WriteRowNames', true);
end

