clc; 
clear; 
close all;

%% ------- Load Data -------------------------------------------------------

filename = 'selected_data.xlsx';
data = readtable(filename);

% Ensure the 'meanspeed' column is numeric
if iscell(data.meanspeed)
    data.meanspeed = cellfun(@str2double, data.meanspeed);
end

% Ensure 'x' and 'y' columns are numeric
if iscell(data.x)
    x = cellfun(@str2double, data.x);
else
    x = data.x;
end

if iscell(data.y)
    y = cellfun(@str2double, data.y);
else 
    y = data.y;
end

%% ------- Initialize variables --------------------------------------------

n = height(data);           % Number of glowworms
s = 0.8;                    % Step size
L0 = 5;                     % Initial luciferin
r0 = 50;                    % Initial decision range
rho = 0.4;                  % Luciferin decay constant
gamma = 0.6;                % Luciferin enhancement constant
B = 0.08;                   % Decision range update constant
rs = 30;                    % Maximum decision range
nt = 5;                     % Threshold for decision range update
maxIter = 150;             
showplot = true;

Agent = [x, y];                    
luciferin = zeros(n, maxIter);      
luciferin(:, 1) = L0;               
decision_range = r0 * ones(n, 1);

%% ------- Calculate the fitness value -------------------------------------

L = data.length;
AVt = data.meanspeed;

Fitness = L ./ AVt;


%% ------- Start of plot ---------------------------------------------------

figure; 
hold on;
xlabel('X');
ylabel('Y');
grid on;

% Plot initial positions
h = plot(x, y, 'o'); % Handle to the plot
pause(0.2)
hold on;


%% ------- Iterasi ---------------------------------------------------------

for t = 1:maxIter
    % Update Luciferin
    if t > 1
        for i = 1:n
            luciferin(i,t) = (1 - rho) * luciferin(i,t-1) + gamma * Fitness(i);
        end
    end
    
    % Moving the Glow-worms
    for ii = 1:n
        curAgent = Agent(ii,:);
        curLuciferin = luciferin(ii, t);
        distance = EuclidDistance(Agent, repmat(curAgent, n, 1));
        
        Ni = find((distance < decision_range(ii)) & (luciferin(:, t) > curLuciferin));
        
        if isempty(Ni)  % If no glow-worm exists within its local range
            Agent(ii,:) = curAgent;
        else
            localRangeLuciferin = luciferin(Ni, t);
            localRangeAgent = Agent(Ni, :);

            probs = (localRangeLuciferin - curLuciferin) / sum(localRangeLuciferin - curLuciferin);
          
            selectedJ = localRangeAgent(SelectByRoulette(probs), :);
            
            Agent(ii, :) = curAgent + s * (selectedJ - curAgent) / EuclidDistance(selectedJ, curAgent);
        end
        neighborSz = length(Ni);
        decision_range(ii) = min([rs, max([0, decision_range(ii) + B * (nt - neighborSz)])]);
    end
    
    % Update plot
    set(h, 'XData', Agent(:,1), 'YData', Agent(:,2));
    drawnow; % Force MATLAB to update the plot
    pause(0.1); % Add a pause to see the movement clearly
end

%% ------- Helper functions ------------------------------------------------

function ret = EuclidDistance(pos1, pos2)
    ret = sqrt(sum((pos1 - pos2) .^ 2, 2));
end

function idx = SelectByRoulette(probs)
    cumulativeProbs = cumsum(probs);
    r = rand();
    idx = find(cumulativeProbs >= r, 1);
end
