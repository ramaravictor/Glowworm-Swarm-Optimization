% clc; clear all; close all;

% Load the Excel file
filename = 'simulasisumo.xlsx';
data = readtable(filename);

% Ensure the 'time' column is numeric
if iscell(data.time)
    data.time = cellfun(@str2double, data.time);
end

% Find the most frequent 'time' value and filter data
[most_frequent_time, ~] = mode(data.time);
selected_data = data(data.time == most_frequent_time, :);

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
s = 0.8;                   % Step size
L0 = 5;                    % Initial luciferin
r0 = 30;                   % Initial decision range
rho = 0.4;                 % Luciferin decay constant
gamma = 0.6;               % Luciferin enhancement constant
B = 0.08;                  % Decision range update constant
rs = 30;                   % Maximum decision range
nt = 5;                    % Threshold for decision range update
maxIter = 150;             % Maximum number of iterations
showplot = true;

L = 777.91;                % Length of the road
data_speed = selected_data.speed;
Agent = [x, y];            % Initializes the Search Agents (SA)

luciferin = L0 * ones(n, 1);
rd = r0 * ones(n, 1);


figure; % Open figure for plotting
hold on;
xlabel('X');
ylabel('Y');
grid on;

% Plot initial positions
plot(x, y, 'o');
hold off;
