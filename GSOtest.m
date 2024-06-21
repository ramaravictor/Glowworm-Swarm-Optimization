clc; 
clear all; 
close all;

% Load the Excel file
filename = 'dataset.xlsx';
data = readtable(filename);

% Ensure the 'time' column is numeric
if iscell(data.time)
    data.time = cellfun(@str2double, data.time);
end

% Find the most frequent 'time' value and filter data
[most_frequent_time, ~] = mode(data.time);
selected_data = data(data.time == most_frequent_time, :);

% Remove duplicate IDs
[~, unique_idx] = unique(selected_data.id);
selected_data = selected_data(unique_idx, :);

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

% Ensure the 'meanspeed' column is numeric
if iscell(selected_data.meanspeed)
    selected_data.meanspeed = cellfun(@str2double, selected_data.meanspeed);
end

% Ensure 'x' and 'y' columns are numeric
if iscell(selected_data.x)
    x = cellfun(@str2double, selected_data.x);
else
    x = selected_data.x;
end

if iscell(selected_data.y)
    y = cellfun(@str2double, selected_data.y);
else 
    y = selected_data.y;
end

% Initialize variables
n = height(selected_data);  % Number of glowworms set to the number of rows in data_at_time
s = 0.8;                    % Step size
L0 = 5;                     % Initial luciferin
r0 = 30;                    % Initial decision range
rho = 0.4;                  % Luciferin decay constant
gamma = 0.6;                % Luciferin enhancement constant
B = 0.08;                   % Decision range update constant
rs = 30;                    % Maximum decision range
nt = 5;                     % Threshold for decision range update
maxIter = 150;              % Maximum number of iterations
showplot = true;

L = selected_data.length;
AVt = selected_data.meanspeed;

Fitness = L ./ AVt;

Agent = [x, y];             % Initializes the Search Agents (SA)

luciferin = zeros(n, maxIter);  % Initialize luciferin levels for all glowworms and time steps
luciferin(:, 1) = L0;           % Set initial luciferin levels to L0
rd = r0 * ones(n, 1);

figure; % Open figure for plotting
hold on;
xlabel('X');
ylabel('Y');
grid on;

% Plot initial positions
plot(x, y, 'o');
hold off;

% Main simulation loop
for t = 1:maxIter
    % Update Luciferin
    if t > 1
        for i = 1:n
            luciferin(i,t) = (1 - rho) * luciferin(i,t-1) + gamma * Fitness(i);
        end
    end
    
end

disp(luciferin);  % Display the luciferin levels for verification
