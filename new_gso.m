clc; clear all; close all;

% Load the Excel file
filename = 'Hsimulasi.xlsx';
data = readtable(filename);

% Filter data for time 86
data_at_time = data(data.time == 86, :);

% Select columns 'id', 'x', 'y', and 'speed'
selected_data = data_at_time(:, {'id', 'x', 'y', 'speed'});
selected_data.id = cellfun(@(x) str2double(x(3:end)), cellstr(selected_data.id));

% Initialize variables
n = 75;           % Number of glowworms
s = 0.8;          % Step size
L0 = 5;           % Initial luciferin
r0 = 30;          % Initial decision range
rho = 0.4;        % Luciferin decay constant
gamma = 0.6;      % Luciferin enhancement constant
B = 0.08;         % Decision range update constant
rs = 30;          % Maximum decision range
nt = 5;           % Threshold for decision range update
maxIter = 150;    % Maximum number of iterations
L = 777.91;       % Length of the road
showplot = true;

% Extract positions
x = selected_data.x;
y = selected_data.y;

% Calculate the fitness for each vehicle
epsilon = 0.001;  % A small constant to prevent division by zero
selected_data.fitness = L ./ max(selected_data.speed, epsilon);

% Initialize positions and luciferin
positions = [x, y]; 
luciferin = L0 * ones(n, 1);
rd = r0 * ones(n, 1);

% Initialize trails for plotting
colors = lines(n); % Different colors for each glowworm
trails = cell(n, 1);
for i = 1:n
    trails{i} = positions(i, :);
end

% Plot initial positions
figure; % Open figure for plotting
hold on;
title('Glowworm Swarm Optimization');
xlabel('X Coordinate');
ylabel('Y Coordinate');
grid on;
plot(positions(:,1), positions(:,2), 'o', 'LineWidth', 2); % Initial positions

% Main simulation loop
for t = 1:maxIter
    old_positions = positions; % Save old positions for plotting trails

    % Update Luciferin
    for i = 1:n
        if ~isinf(selected_data.fitness(i)) % Skip if fitness is inf
            luciferin(i) = (1 - rho) * luciferin(i) + gamma * selected_data.fitness(i);
        end
    end

    % Movement of Glowworm
    for i = 1:n
        distances = sqrt(sum((positions - positions(i,:)).^2, 2));
        Ni = find((distances < rd(i)) & (luciferin(i) < luciferin) & (distances > 0));
        
        if ~isempty(Ni)
            weights = luciferin(Ni) - luciferin(i);
            if any(weights > 0)
                probabilities = weights / sum(weights);
                j = Ni(randsample(length(Ni), 1, true, probabilities));
                direction = (positions(j,:) - positions(i,:)) / norm(positions(j,:) - positions(i,:));
                positions(i,:) = positions(i,:) + s * direction;
            end
        end
    end

    % Update Decision Range
    for i = 1:n
        distances = sqrt(sum((positions - positions(i,:)).^2, 2));
        Ni = find((distances < rd(i)) & (luciferin(i) < luciferin) & (distances > 0));
        n_i = length(Ni);
        rd(i) = min(rs, max(0, rd(i) + B * (nt - n_i)));
    end

    % Draw trails
    for i = 1:n
        plot([old_positions(i,1), positions(i,1)], [old_positions(i,2), positions(i,2)], 'Color', colors(i,:), 'LineWidth', 1); % Thinner lines for paths
    end

    if showplot
        drawnow;
    end
end

% Mark final positions with 'o' for only the glowworms that moved
initial_positions = [x, y];
moved = sqrt(sum((positions - initial_positions).^2, 2)) > 0; % Check if each glowworm moved
scatter(positions(moved, 1), positions(moved, 2), 100, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y'); % Mark final positions

% Draw the complete trails
for i = 1:n
    plot(trails{i}(:, 1), trails{i}(:, 2), 'Color', colors(i,:), 'LineWidth', 0.5);
end

hold off;

% Display final positions and luciferin levels
disp('Final positions:');
disp(positions);
disp('Final luciferin levels:');
disp(luciferin);
