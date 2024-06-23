clc; 
clear all; 
close all;

% Load the Excel file

filename = 'newdataset.xlsx';
sheetName = 'selected_data';
selected_data = readtable(filename, 'Sheet', sheetName);

% % Ensure the 'time' column is numeric
% if iscell(data.time)
%     data.time = cellfun(@str2double, data.time);
% end
% 
% % Find the most frequent 'time' value and filter data
% [most_frequent_time, ~] = mode(data.time);
% selected_data = data(data.time == most_frequent_time, :);
% 
% % Remove duplicate IDs
% [~, unique_idx] = unique(selected_data.id);
% selected_data = selected_data(unique_idx, :);

% XML preprocessing to create a mapping of lane id to length
xml_data = xmlread('Rute.net.xml');

lane_length_map = containers.Map;
edges = xml_data.getElementsByTagName('edge');

for i = 0:edges.getLength-1
    edge = edges.item(i);
    lanes = edge.getElementsByTagName('lane');
    for j = 0:lanes.getLength-1
        lane = lanes.item(j);
        lane_id = strtrim(char(lane.getAttribute('id'))); % Trim whitespace
        lane_length = str2double(lane.getAttribute('length'));
        lane_length_map(lane_id) = lane_length; % Add to map
    end
end

% Add the lane length to the selected data based on lane ID
selected_data.length = zeros(height(selected_data), 1); % Initialize length column
for i = 1:height(selected_data)
    lane_id = strtrim(selected_data.lane{i}); % Trim whitespace
    if isKey(lane_length_map, lane_id)
        selected_data.length(i) = lane_length_map(lane_id);
    else
        warning('Lane ID %s not found in XML', lane_id);
    end
end



% Export the selected_data to an Excel file
writetable(selected_data, 'selected_data.xlsx');

