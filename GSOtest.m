clc; clear all; close all;

% Load the Excel file
filename = 'Hsimulasi.xlsx';
data = readtable(filename);

% Filter data for time 
data_at_time = data(data.time == 86, :);

% Select columns 'id', 'x', 'y', and 'speed'
selected_data = data_at_time(:, {'id', 'x', 'y', 'speed'});

% Convert 'id'
selected_data.id = cellfun(@(x) str2double(x(3:end)), cellstr(selected_data.id));

% length of the road
L = 777.91;

% Calculate the fitness for each vehicle
epsilon = 0.001;  % A small constant to prevent division by zero
min_fitness = 1; % Minimal fitness level for movement
selected_data.fitness = L ./ max(selected_data.speed, epsilon);
selected_data.fitness = max(selected_data.fitness, min_fitness); % Apply minimum fitness

% Introduce small variability to the fitness to avoid equal values
selected_data.fitness = selected_data.fitness + 0.01 * randn(size(selected_data.fitness));

% Initialize options
options = struct('population', 75, 'L0', 5, 'r0', 30, 'rho', 0.4, 'y', 0.6, 'B', 0.08, 's', 0.8, 'rs', 30, 'nt', 5, 'maxIter', 150, 'showplot', true);
L0 = options.L0;
r0 = options.r0;
rho = options.rho;
y = options.y;
B = options.B;
rs = options.rs;
nt = options.nt;
maxIter = options.maxIter;

m = 2; % 2D space for (x, y)
n = options.population;
s = options.s;

% Initialize positions and luciferin
positions = [selected_data.x, selected_data.y]; % ubah 2 dimensi

luciferin = options.L0 * ones(n, 1);
rd = options.r0 * ones(n, 1);
colors = lines(n); % Different colors for each glowworm

% Initialize a cell array to store the trails of each glowworm
trails = cell(n, 1);
for i = 1:n
    trails{i} = positions(i, :);
end

figure; % Open figure for plotting
hold on;
title('Glowworm Swarm Optimization');
xlabel('X Coordinate');
ylabel('Y Coordinate');
grid on;

% Plot initial positions   
plot(positions(:,1), positions(:,2), 'x'); % Initial positions with larger 'x'
hold on;
plot (selected_data.x, selected_data.y, 'o')
hold on;
% Main simulation loop
for t = 1:maxIter
    old_positions = positions;
    old_luciferin = luciferin; % Save the previous luciferin values

    % Update Luciferin
    for i = 1:n
        luciferin(i) = (1 - rho) * old_luciferin(i) + y * selected_data.fitness(i);
    end

    % Movement of Glowworm
    for i = 1:n
        % Calculate Euclidean distance between glowworm i and all other glowworms
        distances = sqrt(sum((positions - positions(i,:)).^2, 2));
        
        % Define neighbors Ni(t) based on the given formula
        Ni = find((distances < rd(i)) & (luciferin(i) < luciferin) & (distances > 0));

        if ~isempty(Ni)
            % Calculate weights and probabilities for neighbors
            weights = luciferin(Ni) - luciferin(i);
            if any(weights > 0)
                probabilities = weights / sum(weights);
                % Select a neighbor j based on the calculated probabilities
                j = Ni(randsample(length(Ni), 1, true, probabilities));
                % Update positions in 2D using the provided formula
                direction = (positions(j,:) - positions(i,:));
                norm_direction = norm(direction);
                if norm_direction > 0
                    positions(i,1) = positions(i,1) + s * (positions(j,1) - positions(i,1)) / norm_direction;
                    positions(i,2) = positions(i,2) + s * (positions(j,2) - positions(i,2)) / norm_direction;
                end
            end
        end

        % Update the trail for the current glowworm
        trails{i} = [trails{i}; positions(i, :)];
    end


     % Update Decision Range
    for i = 1:n
        % Calculate Ni(t) again for updating the decision range
        distances = sqrt(sum((positions - positions(i,:)).^2, 2));
        Ni = find((distances < rd(i)) & (luciferin(i) < luciferin) & (distances > 0));
        n_i = length(Ni);  % Number of neighbors
        
        % Update the decision range based on the formula
        rd(i) = min(options.rs, max(0, rd(i) + options.B * (options.nt - n_i))); % Update decision range
    end


    % Draw trails
    for i = 1:n
        plot([old_positions(i,1), positions(i,1)], [old_positions(i,2), positions(i,2)], 'Color', colors(i,:), 'LineWidth', 1); % Thinner lines for paths
    end

    if options.showplot
        drawnow;
    end
end

% Mark final positions with 'o' for only the glowworms that moved
moved = sqrt(sum((positions - initial_positions).^2, 2)) > 0; % Check if each glowworm moved
scatter(positions(moved, 1), positions(moved, 2), 70, 'o', 'LineWidth', 1, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm'); % Mark final positions with yellow filled circles
hold off;

% Display final positions and luciferin levels
disp('Final positions:');   
disp(positions);
disp('Final luciferin levels:');
disp(luciferin);


