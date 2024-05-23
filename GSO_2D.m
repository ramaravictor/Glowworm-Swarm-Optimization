clc; clear all; close all;


% Multimodal test functions
% Peaks function
J1 = @(X_init) 3*(1 - X_init(:,1)).^2 .* exp(-(X_init(:,1).^2) - (X_init(:,2) + 1).^2) - 10*(X_init(:,1)/5 - X_init(:,1).^3 - X_init(:,2).^5) .* exp(-X_init(:,1).^2 - X_init(:,2).^2) - (1/3) .* exp(-(X_init(:,1) + 1).^2 - X_init(:,2).^2);
% Fungsi Rastrigin 2D
J2 = @(X_init) 20 + X_init(:,1).^2 - 10*cos(2*pi*X_init(:,1)) + X_init(:,2).^2 - 10*cos(2*pi*X_init(:,2)); 
% Circles function
J3 = @(X_init) (X_init(:,1).^2 + X_init(:,2).^2).^0.25 .* (sin(2 * 50 * (X_init(:,1).^2 + X_init(:,2).^2).^0.1).^2 + 1.0);
% Staircase function
J4 = @(x_init) 25 - x_init(1) - x_init(2);
% Equal-peaks-A function
J6 = @(x_init) sum(cos(x_init).^2);
% Equal-peaks-B function
J9 = @(x_init) cos(x_init(1)).^2 + sin(x_init(2)).^2;


% Memanggil fungsi GSO(objf,solSize)
[optX, optFcn, noIter, X_init, agent_x, agent_y] = GSO(J9, 2);

% Fungsi GSO
function [optX,optFcn,noIter, X_init, agent_x, agent_y] = GSO(objf,solSize)


% Menginisialisasi opsi
options.population = 80;
options.L0 = 5;
options.r0 = 30;
options.rho = 0.4;
options.y = 0.6;
options.B = 0.08;
options.s = 0.5;
options.rs = 30;
options.nt = 5;
options.maxIter = 250;
options.showplot = true;


% Inisialisasi variabel
noPop = options.population;
L0 = options.L0;
r0 = options.r0;
rho = options.rho;
y = options.y;
B = options.B;
s = options.s;
rs = options.rs;
nt = options.nt;
maxIter = options.maxIter;

% Mengambil populasi awal dari data Excel
data = xlsread('Hsimulasi.xlsx', 'Sheet2', 'C2:D81');
% Mengambil indeks acak
indeks_acak = randperm(size(data, 1), noPop);
X_init = data(indeks_acak, :);

% X_init = data;


% Fungsi objektif untuk optimasi
objfcn = @(X_init) ConvertToMin(X_init, objf);

% Inisialisasi luciferin dan decision_range
luciferin = L0 * ones(noPop, 1);
decision_range = r0 * ones(noPop, 1);
iter = 1;

% Inisialisasi matriks untuk menyimpan posisi koordinat cacing per iterasi
agent_x = zeros(noPop, maxIter);
agent_y = zeros(noPop, maxIter);

%% ---------- Start of plot --------------------------------------------
if solSize < 3 && options.showplot     % Determines if a plot should be made
    fig = figure;
end
if solSize == 1 && options.showplot
    oneDimPlotStarter(objfcn,data)
elseif solSize == 2
    twoDimPlotStarter(objfcn,data)
end
if solSize < 3 && options.showplot
    pause(0.2)
    hold on
end
if solSize == 1
    hd = oneDimPlotUpdate(objfcn,X_init);
elseif solSize == 2
    hd = twoDimPlotUpdate(objfcn,X_init);
end
%% ------- End of plot -----------------------------------------------------
%% Iterasi
while iter <= maxIter

%      if solSize < 3 && options.showplot           % Updates the plot
%          pause(0.2);
%          delete(hd)
%      end


    % Mengupdate luciferin
    luciferin = (1 - rho) * luciferin + y * getFcn(objfcn, X_init);

    [~, ~] = max(luciferin);

    % Movement-phase Glow-worms
    for ii = 1:noPop
        curX = X_init(ii, :);
        curLuciferin = luciferin(ii);
        
        distFromI = EuclidDistance(X_init, repmat(curX, noPop, 1));

        Ni = find((distFromI < decision_range(ii)) & (luciferin > curLuciferin));

        if isempty(Ni)                          % If no glow-worm exists within its local range
            X_init(ii, :) = curX;
        else
            localRangeL = luciferin(Ni);
            localRangeX = X_init(Ni, :);
            
            probs = (localRangeL - curLuciferin) / sum(localRangeL - curLuciferin);
            selectedPos = SelectByRoulete(probs);
            selectedX = localRangeX(selectedPos, :);
            X_init(ii, :) = curX + s * (selectedX - curX) / EuclidDistance(selectedX, curX);
        end
        neighborSz = length(Ni);
        
        decision_range(ii) = min([rs, max([0, decision_range(ii) + B * (nt - neighborSz)])]); % Pencarian tetangga
    end

    iter = iter + 1;

    % Updating the plot
    if solSize == 1 && options.showplot
        hd = oneDimPlotUpdate(objfcn,X_init);
    elseif solSize == 2
        hd = twoDimPlotUpdate(objfcn,X_init);
    end

    % Menyimpan posisi agent saat ini
    agent_x(:, iter) = X_init(:, 1);
    agent_y(:, iter) = X_init(:, 2);
end

if solSize < 3 && options.showplot
    hold off
end

% optX = Xs(bestPos,:);
optX = X_init;
optFcn = getFcn(objfcn, optX);
noIter = iter - 1;

% Fungsi bantu
    function Fx = getFcn(objfcn, X_init)
        n = size(X_init, 1);
        Fx = ones(n, 1);
        for k = 1:n
            Fx(k) = objfcn(X_init(k, :));
        end
    end

    function ret = EuclidDistance(pos1, pos2)
        ret = sqrt(sum((pos1 - pos2).^2, 2));
    end

    function ret = SelectByRoulete(allProb)
        cumProb = cumsum(allProb);
        rn = rand;
        hd = find(cumProb >= rn);
        ret = hd(1);
    end

    function twoDimPlotStarter(~, X_init)
        xrange = X_init(:, 1);
        yrange = X_init(:, 2);
        plot(xrange, yrange, 'x');          % Tambahkan tanda 'x' pada posisi awal cacing
        grid on;
        hold on;
    end

    function hd = twoDimPlotUpdate(~, X_init)
        xrange = X_init(:, 1);
        yrange = X_init(:, 2);
        hd = plot(xrange, yrange, 'ro', 'markersize', 1, 'markerfacecolor', 'm');
    end

    function minObjFcn = ConvertToMin(X_init, objfcn)
        fcn = objfcn(X_init);
        if fcn >= 0
            minObjFcn = 1 / (1 + fcn);
        else
            minObjFcn = 1 + abs(fcn);
        end
    end
end
