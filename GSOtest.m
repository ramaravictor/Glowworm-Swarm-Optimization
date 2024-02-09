
clc; clear all; close all;
tic

global n m A_init A Ell gamma ro step1 r_d r_s beta r_min n_t Ave_d bound agent_id

m = 2;                          % No. of dimensions

% Parameter initialization
% -----------------------------------------------
n = 75;                        % No. of agents
r_s = 3;                        % Sensor range
r_d = r_s * ones(n, 1);         % Local decision range
r_min = 0;                      % Threshold decision range
gamma = 0.6;                    % Luciferin enhancement constant
ro = 0.4;                       % Luciferin decay constant
step1 = 0.03;                   % Distance moved by each glowworm when a decision is taken
beta = 0.08;                    % Decision range gain
n_t = 5;                        % Desired no. of neighbors

% Initialization of variables
% ----------------------------------------------------
bound = 3;                      % Parameter specifying the workspace range
DeployAgents;                   % Deploy the glowworms randomly
Ell = 5 * ones(n, 1);           % Initialization of Luciferin levels
j = 1;                          % Iteration index
iter = 250;                     % No. of iterations
%Ave_d = zeros(iter, 1);         % Average distance
agent_id = randperm(n)';             % Agent IDs

% Main loop
% ----------------------------------------------
while (j <= iter)
    UpdateLuciferin;            % Update the luciferin levels at glowworms' current positions
    Act;                        % Select a direction and move
    for k = 1 : n               % store the state histories
        agent_x(k, j) = A(k, 1);
        agent_y(k, j) = A(k, 2);
    end
    j = j + 1;
end

toc                             % Display the total computation time

% Plots
% -------------------------------------------------

figure(1);                      % Plot of trajectories of glowworms from their
                                % initial locations to final locations
plot(A_init(:, 1), A_init(:, 2), 'x');
xlabel('X'); ylabel('Y');
hold on;
DefineAxis;
for k = 1 : n
    plot(agent_x(k, :, :), agent_y(k, :, :));
end
DefineAxis;
grid on;
hold on;
plot([-0.0093; 1.2857; -0.46], [1.5814; -0.0048; -0.6292], 'ok');

%figure(2);                      % Plot of final locations of glowworms
%plot(A(:, 1), A(:, 2), '.');
%DefineAxis;
%grid on;
%hold on;
%plot([-0.0093; 1.2857; -0.46], [1.5814; -0.0048; -0.6292], 'ok');

% Calculate distance and time for each agent
distances = zeros(n, 1);
times = zeros(n, 1);
for i = 1:n
    distances(i) = norm(A_init(i, :) - A(i, :));
    times(i) = iter * step1;        % Assuming each iteration takes step1 time
end

% Combine data into a matrix
Result = [agent_id, A_init, distances, times];


% Function 1: DeployAgents.m
% ---------------------------------------------------
function DeployAgents
    global n m A_init A bound
    B = -bound * ones(n, m);
    A_init = B + 2 * bound * rand(n, m);
    A = A_init;
end

% Function 2: UpdateLuciferin.m
% -------------------------------------------------
function UpdateLuciferin
    global n A J Ell gamma ro
    for i = 1 : n
        x = A(i, 1); y = A(i, 2);
        % The Matlab 'Peaks' function is used here. Please replace it with
        % the multimodal function for which peaks are sought
        J(i, :) = 3 * (1 - x)^2 * exp(-(x^2) - (y + 1)^2) - 10 * (x / 5 - x^3 - y^5) * exp(-x^2 - y^2) - 1 / 3 * exp(-(x + 1)^2 - y^2);
        
        Ell(i, :) = (1 - ro) * Ell(i, :) + gamma * J(i, :);
        
    end
end

% Function 3: Act.m
% ---------------------------------------------------
function Act
    global n r_s r_d N N_a beta n_t
    N(:, :) = zeros(n, n);
    N_a(:, :) = zeros(n, 1);
    for i = 1 : n
        FindNeighbors(i);
        FindProbabilities(i);
        Leader(i) = SelectAgent(i);
    end
    for i = 1 : n
        Move(i, Leader(i));
        r_d(i) = max(0, min(r_s, r_d(i) + beta * (n_t - N_a(i))));
    end
end

% Function 4: FindNeighbors.m
% ---------------------------------------------------
function FindNeighbors(i)
    global n m A N r_d N_a Ell
    n_sum = 0;
    for j = 1 : n
        if (j ~= i)
            square_sum = 0;
            for k = 1 : m
                square_sum = square_sum + (A(i, k) - A(j, k))^2;
            end
            d = sqrt(square_sum);
            if (d <= r_d(i)) & (Ell(i) < Ell(j))
                N(i, j) = 1;
                n_sum = n_sum + 1;
            end
        end
    end
    N_a(i) = n_sum;
end

% Function 5: FindProbabilities.m
% ----------------------------------------------
function FindProbabilities(i)
    global n N Ell pb
    Ell_sum = 0;
    for j = 1 : n
        Ell_sum = Ell_sum + N(i, j) * (Ell(j) - Ell(i));
    end
    if (Ell_sum == 0)
        pb(i, :) = zeros(1, n);
    else
        for j = 1 : n
            pb(i, j) = (N(i, j) * (Ell(j) - Ell(i))) / Ell_sum;
        end
    end
end

% Function 6: SelectAgent.m
% ----------------------------------------------
function j = SelectAgent(i)
    global n pb
    bound_lower = 0;
    bound_upper = 0;
    toss = rand;
    j = 0;
    for k = 1 : n
        bound_lower = bound_upper;
        bound_upper = bound_upper + pb(i, k);
        if (toss > bound_lower) & (toss < bound_upper)
            j = k;
            break;
        end
    end
end

% Function 7: Move.m
% ----------------------------------------------
function Move(i, j)
    global A m step1 Ell bound
    if (j ~= 0) & (Ell(i) < Ell(j))
        temp(i, :) = A(i, :) + step1 * Path(i, j);
        flag = 0;
        for k = 1 : m
            if (temp(i, k) < -bound) | (temp(i, k) > bound)
                flag = 1;
                break;
            end
        end
        if (flag == 0)
            A(i, :) = temp(i, :);
        end
    end
end

% Function 8: Path.m
% ----------------------------------------------
function Del = Path(i, j)
    global A m
    square_sum = 0;
    for k = 1 : m
        square_sum = square_sum + (A(i, k) - A(j, k))^2;
    end
    hyp = sqrt(square_sum);
    for k = 1 : m
        Del(:, k) = (A(j, k) - A(i, k)) / hyp;
    end
end

% Function 9: DefineAxis.m
% ----------------------------------------------
function DefineAxis
    global bound
    axis([-bound bound -bound bound]);
    grid on;
end
