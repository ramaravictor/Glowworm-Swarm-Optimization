clc;
clear;
close all;

%% ------- Muat Data -------------------------------------------------------
data = readtable('Hsimulasi.xlsx');  % Ganti dengan nama file jika berbeda
selectedData = data(5:317, :);  % Ambil 102 data teratas

%% ------- Inisialisasi variabel --------------------------------------------
s = 1;                          % Ukuran langkah
L0 = 5;                         % Luciferin awal
r0 = 100;                       % Jangkauan keputusan awal
rho = 0.4;                      % Konstanta peluruhan luciferin
gamma = 0.6;                    % Konstanta peningkatan luciferin
beta = 0.08;                    % Konstanta pembaruan jangkauan keputusan
rs = 100;                       % Jangkauan keputusan maksimum
nt = 10;                        % Ambang batas pembaruan jangkauan keputusan
maxIter = 200;                  % Jumlah iterasi maksimum
showplot = true;                % Flag untuk menampilkan plot

%% ------- Muat Data Lane --------------------------------------------------
xml_data = xmlread('Rute.net.xml');
lanes_xml = xml_data.getElementsByTagName('lane');

for i = 1:height(selectedData)
    current_lane = selectedData.lane{i}; % Lane dari Excel
    
    for j = 0:lanes_xml.getLength-1
        lane = lanes_xml.item(j); % Ambil elemen lane ke-j dari XML
        lane_id = char(lane.getAttribute('id')); % Dapatkan ID lane dari XML
        
        if strcmp(current_lane, lane_id) % Padankan lane di Excel dan XML
            length_attr = lane.getAttribute('length'); % Dapatkan panjang lane dari XML
            selectedData.length(i) = str2double(length_attr); % Simpan panjang ke dalam array
            break;
        end
    end
end

%% ------- Hitung Rata-rata Kecepatan --------------------------------------
selectedData.avgSpeed = zeros(height(selectedData), 1);
uniqueIds = unique(selectedData.id);

for i = 1:length(uniqueIds)
    currentId = uniqueIds{i};
    
    for currentTime = 0:20
        previousData = selectedData(strcmp(selectedData.id, currentId) & selectedData.time <= currentTime, :);
        
        if ~isempty(previousData)
            avgSpeed = mean(previousData.speed);
        else
            avgSpeed = 0;
        end
        
        selectedData.avgSpeed(strcmp(selectedData.id, currentId) & selectedData.time == currentTime) = avgSpeed;
    end
end

%% ------- Hitung Fitness Value --------------------------------------------
selectedData.fitnessValue = zeros(height(selectedData), 1);

for i = 1:height(selectedData)
    lane_length = selectedData.length(i);
    Speed = selectedData.speed(i);
    
    if Speed == 0
        fitnessValue = 0;
    else
        fitnessValue = lane_length / Speed;
    end
    
    selectedData.fitnessValue(i) = fitnessValue;
end

%% ------- Jalankan Loop Utama Berdasarkan Time ----------------------------
uniqueTimes = unique(selectedData.time);  % Dapatkan nilai time yang unik

% Inisialisasi array untuk menyimpan waktu pemilihan Cluster Head
cluster_head_time = zeros(length(uniqueTimes), 1);

% Jalankan Loop Utama Berdasarkan Time
for timeIdx = 1:length(uniqueTimes)
    currentTime = uniqueTimes(timeIdx);
    
    % Mulai pengukuran waktu
    tic;

    % Pilih data untuk waktu saat ini
    currentAgents = selectedData(selectedData.time == currentTime, :);
    n = height(currentAgents);  % Jumlah agen pada waktu ini
    
    % Ekstrak data posisi awal agen dari tabel
    x = currentAgents.x;
    y = currentAgents.y;
    Agent = [x, y];  % Posisi awal agen
    
    % Inisialisasi luciferin dan jangkauan keputusan
    luciferin = zeros(n, 1);  
    luciferin(:) = L0;  % Set luciferin awal
    decision_range = r0 * ones(n, 1);  % Set jangkauan keputusan awal
    
    % Inisialisasi kepercayaan dengan nilai awal 1 untuk semua glowworm
    directTrust = ones(n, 1);  
    reputationTrust = ones(n, 1);  % Untuk trust berbasis reputasi

    % Menyimpan riwayat posisi cacing
    trail = zeros(n, 2, maxIter+1);  % Memastikan maxIter didefinisikan
    trail(:, :, 1) = Agent;
    
    for t = 1:maxIter  % Looping untuk maxIter iterasi per time
        %% Langkah 1: Pembaruan luciferin
        for i = 1:n
            luciferin(i) = (1 - rho) * luciferin(i) + gamma * currentAgents.fitnessValue(i);
        end
        
        %% Langkah 2: Tentukan tetangga
        neighbors = cell(n, 1);
        for i = 1:n
            neighbors{i} = [];
            for j = 1:n
                if i ~= j
                    % Hitung jarak euclidean antara agen i dan j
                    distance = sqrt((Agent(i, 1) - Agent(j, 1))^2 + (Agent(i, 2) - Agent(j, 2))^2);
                    
                    if distance < decision_range(i) && luciferin(j) > luciferin(i)
                        neighbors{i} = [neighbors{i}; j];
                    end
                end
            end
        end
        
        %% Langkah Tambahan: Candidate Centroid Set Construction
        % Agen dengan luciferin lebih tinggi dibandingkan dengan tetangganya
        % akan dipertimbangkan sebagai kandidat pusat cluster (centroid).
        candidate_centroids = [];
        for i = 1:n
            if isempty(neighbors{i}) || all(luciferin(i) >= luciferin(neighbors{i}))
                candidate_centroids = [candidate_centroids; i]; % Simpan indeks agen sebagai kandidat centroid
            end
        end
        
        %% Langkah 3: Hitung probabilitas dan perbarui posisi
        for i = 1:n
            if ~isempty(neighbors{i})
                prob = zeros(numel(neighbors{i}), 1);
                sum_luciferin_diff = sum((directTrust(neighbors{i}) + reputationTrust(neighbors{i})) .* (luciferin(neighbors{i}) - luciferin(i)));
                for k = 1:numel(neighbors{i})
                    j = neighbors{i}(k);
                    prob(k) = (directTrust(j) + reputationTrust(j)) * (luciferin(j) - luciferin(i)) / sum_luciferin_diff;
                end
                % Pilih tetangga berdasarkan probabilitas
                selected_neighbor = neighbors{i}(randsample(1:numel(neighbors{i}), 1, true, prob));
                
                % Perbarui posisi x dan y secara terpisah
                distance = sqrt((Agent(i, 1) - Agent(selected_neighbor, 1))^2 + (Agent(i, 2) - Agent(selected_neighbor, 2))^2);
                Agent(i, 1) = Agent(i, 1) + s * (Agent(selected_neighbor, 1) - Agent(i, 1)) / distance;
                Agent(i, 2) = Agent(i, 2) + s * (Agent(selected_neighbor, 2) - Agent(i, 2)) / distance;
            end
        end

        %% Langkah Tambahan: Pembaruan Trust
        for i = 1:n
            if ~isempty(neighbors{i})
                % Pembaruan direct trust berdasarkan hasil interaksi
                if currentAgents.fitnessValue(i) > 0
                    directTrust(i) = min(1, directTrust(i) + 0.1);  % Tingkatkan direct trust
                else
                    directTrust(i) = max(0, directTrust(i) - 0.1);  % Kurangi direct trust
                end

                % Reputation Trust diperbarui berdasarkan neighbor
                reputationTrust(i) = mean(directTrust(neighbors{i}));  % Reputation didapat dari rata-rata kepercayaan tetangga
            end
        end
        
        %% Langkah 4: Pembaruan jangkauan keputusan
        for i = 1:n
            decision_range(i) = min(rs, max(0, decision_range(i) + beta * (nt - numel(neighbors{i})) * (directTrust(i) + reputationTrust(i))));
        end
        
        % Simpan jejak posisi
        trail(:, :, t+1) = Agent;
        
        %% Plot posisi agen dengan jejaknya
        if showplot
            figure(5); % Plot dalam figure 5
            clf;
            hold on;
            % Plot jejak pergerakan
            for i = 1:n
                plot(squeeze(trail(i, 1, 1:t+1)), squeeze(trail(i, 2, 1:t+1)), '-');
            end
            % Plot posisi akhir agen pada iterasi saat ini
            plot(Agent(:, 1), Agent(:, 2), 'ro');
            
            % Plot centroid candidate
            plot(Agent(candidate_centroids, 1), Agent(candidate_centroids, 2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
            
            % Pastikan sumbu x dan y memiliki skala yang tetap
            xlim([min(selectedData.x), max(selectedData.x)]);
            ylim([min(selectedData.y), max(selectedData.y)]);
            
            axis equal;  % Menjaga skala x dan y tetap proporsional
            title(['Time ', num2str(currentTime), ' - Iterasi ke-', num2str(t)]);
            hold off;
            drawnow;
        end

    end

    % Simpan waktu pemilihan Cluster Head untuk waktu ini
    cluster_head_time(timeIdx) = toc;
end

% Plot hasil waktu pemilihan Cluster Head
figure;
plot(uniqueTimes, cluster_head_time, '-o');
title('Waktu Pemilihan Cluster Head untuk Setiap Data Time');
xlabel('Time (detik)');
ylabel('Waktu Pemilihan (detik)');
grid on;

