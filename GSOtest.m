clc;
clear;
close all;


%% ------- Muat Data -------------------------------------------------------
filename = 'Hsimulasi.xlsx';
sheet = 'Sheet2';
data = readtable(filename, 'Sheet', sheet);

t = 1;
maxIterations = height(data);

% Inisialisasi tabel untuk menyimpan hasil
result = table('Size', [317, 8], ...
    'VariableTypes', {'double', 'double', 'string', 'double', 'double', 'string', 'double', 'string'}, ...
    'VariableNames', {'t', 'd', 'id', 'x', 'y', 'lane', 'speed', 'type'});

% Inisialisasi matriks untuk menyimpan jarak antar titik
jarakAntarTitik = zeros(maxIterations, maxIterations);

while t + 1 <= maxIterations
    % Increment t
    t = t + 1;

    % Kalkulasi nilai d hanya untuk titik tertentu
    d = sqrt((data.x(t) - data.x(t- 1)).^2 + (data.y(t) - data.y(t- 1)).^2);

    % Menyimpan nilai t, d, id, x, dan y ke dalam result
    result.t(t) = data.time(t);
    result.d(t) = d;
    result.id{t} = data.id{t};
    result.x(t) = data.x(t);
    result.y(t) = data.y(t);
    result.speed(t) = data.speed(t);
    result.type(t) = data.type(t);
    result.lane(t) = data.lane(t);
    result.RREPSN = zeros(height(result), 1);
    result(result.t == 0, :) = [];

    % Menyimpan jarak antar titik ke dalam matriks
    jarakAntarTitik(t-1, t) = d;
    jarakAntarTitik(t, t-1) = d;
end

% Inisialisasi variabel baru untuk menyimpan data
group = table('Size', [100, 1], ...
    'VariableTypes', {'cell'}, ...
    'VariableNames', {'Result'});

% Menggabungkan data t dan id menjadi data baru 'sequence' di tabel result
result.sequence = strcat(string(result.id), '_', string(result.t));

% Inisialisasi struktur untuk menyimpan jumlah kemunculan setiap ID pada setiap iterasi
id_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
id_count = containers.Map('KeyType', 'char', 'ValueType', 'double');

for t = 1:max(result.t)
    % Mendapatkan ID yang muncul pada iterasi saat ini
    ids_current = unique(result.id(result.t == t));
    
    % Loop melalui setiap ID yang muncul pada iterasi saat ini
    for id_idx = 1:numel(ids_current)
        id = ids_current{id_idx};
        % Jika ID tidak ada dalam struktur id_count, tambahkan dan atur nilai awalnya menjadi 0
        if ~isKey(id_count, id)
            id_count(id) = 0;
        end
        % Mendapatkan jumlah kemunculan ID pada iterasi sebelumnya
        count_prev = id_count(id);
        
        % Mendapatkan indeks ID pada iterasi saat ini
        idx_current = find(strcmp(result.id, id) & result.t == t);
        
        % Memperbarui sequence untuk ID pada iterasi saat ini dengan indeks unik yang tepat
        for i = 1:numel(idx_current)
            % Mengubah tipe data SSN menjadi integer dan memulai pengurutan dari time 1
            result.SSN(idx_current(i)) = count_prev + i;
        end
        
        % Mengupdate jumlah kemunculan ID
        id_count(id) = count_prev + numel(idx_current);
    end
end

% Inisialisasi variabel
numNodes = height(unique(result));
validDValues = zeros(numNodes, numNodes);

% Tentukan jumlah baris yang ingin digunakan
jumlah_baris = 313;

% Ambil sejumlah baris tertentu dari tabel result
data_terbatas = result(1:jumlah_baris, :);

% Mengambil jumlah unik dari kolom 'id' dalam tabel 'data_terbatas' untuk mendapatkan jumlah node
numNodes = numel(unique(data_terbatas.sequence));

% Menginisialisasi matriks validDValues dengan jarak antar node
for i = 1:numNodes
    for j = 1:numNodes
        % Perhitungan jarak antar node i dan j
        validDValues(i, j) = sqrt((data_terbatas.x(i) - data_terbatas.x(j))^2 + (data_terbatas.y(i) - data_terbatas.y(j))^2);
    end
end

%% ----------------- Inisialisasi AODV ------------------------------------------
status = '!';
dist = inf(1, numNodes);
next = zeros(1, numNodes);

% Inisialisasi status, dist, dan next
for i = 1:numNodes
    if i == 1
        status(i) = '!';
        dist(i) = 0;
        next(i) = 0;
    else
        status(i) = '?';
        % Gunakan hasil perhitungan jarak dari tabel result
        dist(i) = data_terbatas.d(i);
        next(i) = 1;
    end
end

% Inisialisasi variabel lainnya
flag = 0;
temp = 0;

% Set goalNode
goalNode = 20; % Sesuaikan dengan node tujuan

% Initialize variables to store ping information
pingResults = cell(numNodes, numNodes);
rrepsn = zeros(max(result.t), numel(unique(result.id)));
threshold_lower = -14;
threshold_upper = 14;
result.Status = repmat("Connected", height(result), 1);

%% ----------------- Main loop -----------------------------------
while flag ~= 1 && temp < numNodes
    temp = temp + 1; % Tambahkan iterasi

    % Pilih node dengan dist terkecil dan status '?'
    [minDist, vert] = min(dist(status == '?'));

    % Perbarui status
    status(vert) = '!';

    % Perbarui dist dan next untuk node tetangga
    for i = 1:numNodes
        if status(i) == '?' && dist(i) > dist(vert) + validDValues(vert, i)
            dist(i) = dist(vert) + validDValues(vert, i);
            next(i) = vert;

            % Simulasi RREQ hanya jika tidak dalam keadaan Timeout
            if validDValues(vert, i) < 300
                pingResults{vert, i} = 'Ping: Reply 100%';

                % Log RREQ
                disp(['Node ' num2str(vert) ' sends RREQ message to node ' num2str(i)]);
                % Simulasikan penerimaan RREQ dan kirimkan RREP
                % Update RREPSN values
                rrepsn(i) = rrepsn(i) + 1; % Tingkatkan nilai rrepsn untuk node yang membalas
                % Update tableSSN with RREPSN values
                result.RREPSN(i) = rrepsn(i); % Update nilai RREPSN untuk node yang membalas
                disp(['Node ' num2str(i) ' sends RREP message (RREPSN=' num2str(rrepsn(i)) ') to node ' num2str(vert)]);

                if i == goalNode
                    flag = 1;
                    break;
                end
                
            else
                % Jika mencapai Timeout, set hasil ping menjadi Timeout
                pingResults{vert, i} = 'Timeout';
                pingResults{vert, i} = ['Node ' num2str(vert) ' timeout to Node ' num2str(i)];
                disp(pingResults{vert, i});
                result.Status(vert) = "Timeout";

            end

            % Tambahkan kondisi untuk keluar dari loop jika goalNode tercapai
            if i == goalNode
                flag = 1;
                break;
            end
        end
    end
    result.Difference = result.RREPSN - result.SSN;
    result.Status(result.Difference < threshold_lower | result.Difference > threshold_upper) = "Disconnected";
    result.Status(contains(result.Status, 'Timeout')) = "Timeout";

    if all(status == '!')
        flag = 1;
        break;
    end
end


%% ----- Check for nodes that initiated RREQ but did not receive RREP (Timeout)---------
for i = 1:numNodes
    initiatedRREQ = find(~cellfun('isempty', pingResults(i, :)));
    
    % Initialize receivedRREP as an empty array
    receivedRREP = [];
    
    % Loop through each node's ping result at time t
    for j = 1:numNodes
        % Check if there's a ping result for the current node at time t
        if ~isempty(pingResults{i, j})
            % Check if the ping result indicates a successful reply at time t
            if contains(pingResults{i, j}, 'Ping: Reply 100%')
                receivedRREP = [receivedRREP j];
            end
        end
    end
end

% Inisialisasi variabel untuk menyimpan rute
i = goalNode; % Ganti dengan goalNode
count = 1;
route(count) = goalNode;

% Bangun rute dari node terakhir ke node pertama
while next(i) ~= 0 % Ganti dengan node awal
    count = count + 1;
    route(count) = next(i);
    i = next(i);
end


%% ----------------------------------------------------
% Membuat loop untuk mengecek setiap nilai t
for t = 1:max(result.t)
    % Mendapatkan ID yang muncul pada iterasi saat ini
    ids_current = unique(result.id(result.t == t));
    
    % Loop melalui setiap ID yang muncul pada iterasi saat ini
    for id_idx = 1:numel(ids_current)
        id = ids_current{id_idx};
        % Jika ID tidak ada dalam struktur id_counts, tambahkan dan atur nilai awalnya menjadi 0
        if ~isKey(id_counts, id)
            id_counts(id) = 0;
        end
        % Mendapatkan jumlah kemunculan ID pada iterasi sebelumnya
        count_prev = id_counts(id);
        
        % Mendapatkan indeks ID pada iterasi saat ini
        idx_current = find(strcmp(result.id, id) & result.t == t);
        
        % Memperbarui sequence untuk ID pada iterasi saat ini dengan indeks unik yang tepat
        for i = 1:numel(idx_current)
            result.sequence{idx_current(i)} = [id, '_', num2str(count_prev + i)];
        end
        
        % Mengupdate jumlah kemunculan ID
        id_counts(id) = count_prev + numel(idx_current);
    end
end

%% ----------------------------------------------------
for t = 1:50
    % Mengambil data dengan nilai 't' sesuai iterasi
    resultTime = result(result.t == t, :);

    % Perhitungan nilai d
    if t > 1
        d = sqrt((data.x(t) - data.x(t-1)).^2 + (data.y(t) - data.y(t-1)).^2);
    else
        d = 0; 
    end
    
    % Jika data tidak mencapai 80 baris, tambahkan baris dengan nilai 0
    if size(resultTime, 1) < 317
        rowsTotal = 317 - size(resultTime, 1);
        rowsZero = array2table(zeros(rowsTotal, width(resultTime)), 'VariableNames', resultTime.Properties.VariableNames);
        resultTime = [resultTime; rowsZero];
    end

    % Simpan resultTime ke dalam group
    group.Result{t} = resultTime;

    % Hapus variabel yang tidak ingin ditampilkan di workspace
    clear nonZeroDIdx rowsTotal rowsZero;
end


%% ----------------------------------------------------

for t = 1:50
    % Mengambil tabel dari dalam cell array
    resultTableTime = group.Result{t};

    % Menambahkan kolom warna ke dalam tabel hanya jika d > 0
    resultTableTime.color = cell(height(resultTableTime), 1);

    % Temukan indeks baris dengan nilai d terkecil dan terbesar
    minD = find(resultTableTime.d == min(resultTableTime.d(resultTableTime.d > 0)), 1, 'first');
    maxD = find(resultTableTime.d >= 300);

    % Berikan warna hijau untuk nilai d terkecil jika d > 0
    if ~isempty(minD)
        resultTableTime.color{minD} = 'green';
    end
    
    % Berikan warna merah untuk nilai d terbesar jika d > 0
    if ~isempty(maxD)
        resultTableTime.color(maxD)= {'magenta'};
        % Ubah status menjadi 'Timeout'
        resultTableTime.Status(maxD) = {'Timeout'};
    end
    
    % Isi nilai biru hanya untuk baris dengan nilai d sama dengan 0
    zeroDIdx = resultTableTime.d == 0;
    
    % Hapus node biru dengan nilai d = 0 dari hasil plot
    resultTableTime(zeroDIdx, :) = [];
    
    % Isi nilai biru untuk baris dengan nilai d tidak sama dengan 0 dan tidak memiliki warna
    nonZeroDIdx = find(resultTableTime.d > 0 & cellfun('isempty', resultTableTime.color));
    resultTableTime.color(nonZeroDIdx) = {'blue'};
    
    % Menyimpan indeks baris dengan nilai d terkecil sebagai Head Cluster (warna hijau)
    headClusterIdx = find(strcmp(resultTableTime.color, 'green'));
    if ~isempty(headClusterIdx)
        resultTableTime.color{headClusterIdx} = 'Head Cluster';
    end

    % Menghasilkan nilai pt dalam rentang [200, 300] berdasarkan t
    pt = 50 + (t - 1) * 10; % Pertambahan 10 setiap iterasi t
    
    % Pastikan pt tidak melebihi 300
    if pt > 100
        pt = 100;
    end
    
    % Membuat kolom pt untuk setiap baris
    resultTableTime.pt = repmat(pt, height(resultTableTime), 1);
    
    % Mengatur semua nilai dalam rt menjadi 40
    rt = repmat(20, height(resultTableTime), 1);

    resultTableTime.rt = rt;


    % Menyimpan tabel yang telah dimodifikasi ke dalam cell array
    group.Result{t} = resultTableTime;

    % Hapus variabel yang tidak ingin ditampilkan di workspace
    clear nonZeroDIdx zeroDIdx;
    clear headClusterIdx maxD minD;
end

% Inisialisasi variabel baru untuk warna pada result
result.color = cell(height(result), 1);


%% ----------------------------------------------------

for t = 1:50
    % Mengambil tabel dari dalam cell array
    resultTableTimeSerangan = group.Result{t};

    % Menambahkan kolom warna ke dalam tabel hanya jika d > 0
    resultTableTimeSerangan.color = cell(height(resultTableTimeSerangan), 1);

    % Temukan indeks baris dengan nilai d terkecil dan terbesar
    minD = find(resultTableTimeSerangan.d == min(resultTableTimeSerangan.d(resultTableTimeSerangan.d > 0)), 1, 'first');
    maxD = find(resultTableTimeSerangan.d >= 300);

    % Berikan warna hijau untuk nilai d terkecil jika d > 0
    if ~isempty(minD)
        resultTableTimeSerangan.color{minD} = 'green';
    end
    
    % Berikan warna merah untuk nilai d terbesar jika d > 0
    if ~isempty(maxD)
        resultTableTimeSerangan.color(maxD)= {'magenta'};
        % Ubah status menjadi 'Timeout'
        resultTableTimeSerangan.Status(maxD) = {'Timeout'};
    end
    
    % Isi nilai biru hanya untuk baris dengan nilai d sama dengan 0
    zeroDIdx = resultTableTimeSerangan.d == 0;
    
    % Hapus node biru dengan nilai d = 0 dari hasil plot
    resultTableTimeSerangan(zeroDIdx, :) = [];
    
    % Isi nilai biru untuk baris dengan nilai d tidak sama dengan 0 dan tidak memiliki warna
    nonZeroDIdx = find(resultTableTimeSerangan.d > 0 & cellfun('isempty', resultTableTimeSerangan.color));
    resultTableTimeSerangan.color(nonZeroDIdx) = {'blue'};
    
    
    % Menyimpan indeks baris dengan nilai d terkecil sebagai Head Cluster (warna hijau)
    % headClusterIdx = find(strcmp(resultTableTimeSerangan.color, 'green'));
    % if ~isempty(headClusterIdx)
    %     resultTableTimeSerangan.color{headClusterIdx} = 'Head Cluster';
    % end
    

    % Menyimpan indeks baris dengan status "Disconnected" dan mengubah warna menjadi merah
    disconnectedIdx = find(strcmp(resultTableTimeSerangan.Status, 'Disconnected'));
    if ~isempty(disconnectedIdx)
        resultTableTimeSerangan.color(disconnectedIdx) = {'red'};
    end

    % Check if the status is "Disconnected" and update SSN accordingly
    disconnectedIdx = find(result.Status == "Disconnected");
    for idx = 1:numel(disconnectedIdx)
        % Generate a random RREPSN for disconnected nodes
        result.RREPSN(disconnectedIdx(idx)) = randi([0, 1000000000]); % Assuming the range for RREPSN
    end

    % Inisialisasi matriks koneksi
    resultTableTimeSerangan.koneksi = zeros(size(resultTableTimeSerangan, 1), size(resultTableTimeSerangan, 1));
    
    % Mendapatkan indeks node yang belum terkoneksi
    unconnectedNodesIdx = find(sum(resultTableTimeSerangan.koneksi, 2) == 0);
    
    % Urutkan node yang belum terkoneksi berdasarkan nilai d dari terkecil hingga terbesar
    [~, sortedIdx] = sort(resultTableTimeSerangan.d(unconnectedNodesIdx));
    sortedUnconnectedNodesIdx = unconnectedNodesIdx(sortedIdx);
    
    % Membuat koneksi ulang berdasarkan node yang tidak terkoneksi yang sudah diurutkan
    for i = 1:length(sortedUnconnectedNodesIdx)
        currentNode = sortedUnconnectedNodesIdx(i);
        for j = (i+1):length(sortedUnconnectedNodesIdx)
            nextNode = sortedUnconnectedNodesIdx(j);
            if resultTableTimeSerangan.d(nextNode) < 300 % Jika jarak antara node saat ini dengan node berikutnya kurang dari 300
                if sum(resultTableTimeSerangan.koneksi(currentNode, :)) < 2 && sum(resultTableTimeSerangan.koneksi(nextNode, :)) < 2 % Pastikan setiap node memiliki maksimal 2 koneksi
                    resultTableTimeSerangan.koneksi(currentNode, nextNode) = 1;
                    resultTableTimeSerangan.koneksi(nextNode, currentNode) = 1;
                    break; % Hanya satu koneksi yang perlu ditambahkan
                end
            end
        end
    end
    
    % Nonaktifkan koneksi ke dan dari node-node merah
    redNodesIdx = find(strcmp(resultTableTimeSerangan.color, 'red'));
    if ~isempty(redNodesIdx)
        for i = 1:length(redNodesIdx)
            redNode = redNodesIdx(i);
            resultTableTimeSerangan.koneksi(redNode, :) = 0; % Nonaktifkan koneksi ke node lain
            resultTableTimeSerangan.koneksi(:, redNode) = 0; % Nonaktifkan koneksi dari node lain
        end
    end
    
    % Membuat koneksi ulang berdasarkan node yang tidak terkoneksi dan bukan berwarna merah
    for i = 1:size(resultTableTimeSerangan.koneksi, 1)
        if sum(resultTableTimeSerangan.koneksi(i, :)) == 0 && ~strcmp(resultTableTimeSerangan.color(i), 'red') % Jika node belum terkoneksi dengan siapa pun dan bukan berwarna merah
            for j = 1:size(resultTableTimeSerangan.koneksi, 2)
                if i ~= j && sum(resultTableTimeSerangan.koneksi(j, :)) < 2 && resultTableTimeSerangan.d(i) < 300 && resultTableTimeSerangan.d(j) < 300
                    if sum(resultTableTimeSerangan.koneksi(i, :)) < 2 && sum(resultTableTimeSerangan.koneksi(j, :)) < 2 % Pastikan setiap node memiliki maksimal 2 koneksi
                        resultTableTimeSerangan.koneksi(i, j) = 1;
                        resultTableTimeSerangan.koneksi(j, i) = 1;
                        break; % Hanya satu koneksi yang perlu ditambahkan
                    end
                end
            end
        end
    end
    
    % Menghubungkan node "Disconnected" ke node "Connected"
    disconnectedNodesIdx = find(strcmp(resultTableTimeSerangan.Status, 'Disconnected'));
    connectedNodesIdx = find(~strcmp(resultTableTimeSerangan.Status, 'Disconnected') & ~strcmp(resultTableTimeSerangan.color, 'red'));
    
    for i = 1:length(disconnectedNodesIdx)
        disconnectedNode = disconnectedNodesIdx(i);
        for j = 1:length(connectedNodesIdx)
            connectedNode = connectedNodesIdx(j);
            if resultTableTimeSerangan.d(disconnectedNode) < 300 && resultTableTimeSerangan.d(connectedNode) < 300
                if sum(resultTableTimeSerangan.koneksi(disconnectedNode, :)) < 2 && sum(resultTableTimeSerangan.koneksi(connectedNode, :)) < 2 % Pastikan setiap node memiliki maksimal 2 koneksi
                    resultTableTimeSerangan.koneksi(disconnectedNode, connectedNode) = 1;
                    resultTableTimeSerangan.koneksi(connectedNode, disconnectedNode) = 1;
                    break; % Hanya satu koneksi yang perlu ditambahkan
                end
            end
        end
    end

    % Menghasilkan nilai pt dalam rentang [200, 300] berdasarkan t
    pt = 50 + (t - 1) * 10; % Pertambahan 10 setiap iterasi t
    
    % Pastikan pt tidak melebihi 300
    if pt > 100
        pt = 100;
    end
    
    % Membuat kolom pt untuk setiap baris
    resultTableTimeSerangan.pt = repmat(pt, height(resultTableTimeSerangan), 1);
  
    % Mengatur semua nilai dalam rt menjadi 40
    rt = repmat(15, height(resultTableTimeSerangan), 1);

    resultTableTimeSerangan.rt = rt;

    % Set pt dan rt menjadi 0 untuk node yang memiliki warna merah
    if ~isempty(redNodesIdx)
        for i = 1:length(redNodesIdx)
            redNode = redNodesIdx(i);
            resultTableTimeSerangan.pt(redNode, :) = -0.5;
            resultTableTimeSerangan.rt(redNode, :) = -0.5;
        end
    end
    
    % Menyimpan tabel yang telah dimodifikasi ke dalam cell array
    group.ResultTime{t} = resultTableTimeSerangan;

    % Hapus variabel yang tidak ingin ditampilkan di workspace
    clear nonZeroDIdx zeroDIdx;
    clear headClusterIdx maxD minD;
end

%% ------------ P L O T T I N G ----------------------------------------------------------
warna = {'blue', 'red', 'green', 'black', 'cyan', 'magenta', 'yellow', 'white'};

% Inisialisasi delay dan throughput
delay1 = zeros(1, 100);
throughput1 = zeros(1, 100);

% Inisialisasi delay dan throughput
delay2 = zeros(1, 100);
throughput2 = zeros(1, 100);

% Membuat plot untuk setiap nilai t dari 1 hingga 20
for t_idx = 1:20
    
    % Mengambil tabel dari dalam cell array untuk plot kedua
    resultTableTimeSerangan = group.Result{t};

    % Hitung jumlah node merah
    redNodeCount = sum(strcmp(resultTableTimeSerangan.color, 'red'));

    % Membersihkan figur pertama sebelum memplot iterasi berikutnya
    figure(1);
    clf;
    axis([-50 350 -40 120]);
    title(['Simulasi 1 Tanpa Serangan - Iterasi ', num2str(t_idx)]);
    xlabel('Data x');
    ylabel('Data y');
    grid on;
    hold on;

    % Membersihkan figur kedua sebelum memplot iterasi berikutnya
    figure(2);
    clf;
    axis([-50 350 -40 120]);
    title(['Simulasi 2 Serangan - Iterasi ', num2str(t_idx), ' - Malicious Nodes: ', num2str(redNodeCount)]);
    xlabel('Data x');
    ylabel('Data y');
    grid on;
    hold on;

    % Membersihkan figur delay sebelum memplot iterasi berikutnya
    figure(3);
    axis('auto');
    xlabel('Jumlah Kendaraan (s)');
    ylabel('Delay (ms)');
    grid on;
    hold on;

    % Membersihkan figur throughput sebelum memplot iterasi berikutnya
    figure(4);
    axis('auto');
    xlabel('Jumlah Kendaraan (s)');
    ylabel('Throughput (kbps)');
    grid on;
    hold on;

    % Mengambil tabel dari dalam cell array untuk plot pertama
    resultTableTime = group.Result{t_idx};

    % Urutkan berdasarkan nilai d
    [~, idxSorted] = sort(resultTableTime.d);
    resultTableTime = resultTableTime(idxSorted, :);

    for i = 1:size(resultTableTime, 1)
        if strcmp(resultTableTime.color{i}, 'Head Cluster')
            figure(1);
            % scatter(resultTableTime.x(i), resultTableTime.y(i), 100, 'green', 'X', 'LineWidth', 1.5); % Symbol X for Head Cluster
        elseif strcmp(resultTableTime.color{i}, 'blue')
            figure(1);
            scatter(resultTableTime.x(i), resultTableTime.y(i), 64, 'blue', 'o', 'filled'); % Blue dots
        elseif strcmp(resultTableTime.color{i}, 'magenta')
            figure(1);
            scatter(resultTableTime.x(i), resultTableTime.y(i), 64, 'magenta', 'o', 'filled'); % Magenta dots
        end

        % Plot garis antar node
        if i < size(resultTableTime, 1)
            figure(1);
            plot([resultTableTime.x(i), resultTableTime.x(i+1)], [resultTableTime.y(i), resultTableTime.y(i+1)], 'b--', 'LineWidth', 1);
        end
    end

    % Menambahkan legenda untuk subplot pertama
    figure(1);
    hold on; 
    % h1 = scatter(NaN, NaN, 100, 'green', 'X', 'LineWidth', 1.5); 
    h2 = scatter(NaN, NaN, 64, 'blue', 'o', 'filled'); 
    h3 = scatter(NaN, NaN, 64, 'magenta', 'o', 'filled');
    h4 = scatter(NaN, NaN, 64, 'red', 'o', 'filled');
    leg1 = legend([h2], 'Node Kendaraan', 'Location', 'northeast');
    set(leg1, 'Box', 'on');
    hold off;

    % Plot data pada subplot pertama
    figure(2);
    xlabel('Data x');
    ylabel('Data y');
    clf;
    grid on;
    hold on;

    % Mengambil tabel dari dalam cell array untuk plot kedua
    resultTableTimeSerangan = group.ResultTime{t_idx};

    % Menentukan newHeadCluster berdasarkan nilai d terkecil yang tidak 'red'
    minD = min(resultTableTimeSerangan.d(~strcmp(resultTableTimeSerangan.color, 'red')));
    newHeadClusterIndex = find(resultTableTimeSerangan.d == minD, 1);

    % Iterate over all nodes to plot them based on their properties
    for i = 1:size(resultTableTimeSerangan, 1)
        if i == newHeadClusterIndex
            figure(2);
            scatter(resultTableTimeSerangan.x(i), resultTableTimeSerangan.y(i), 100, 'g', 'X', 'LineWidth', 1.5);
        elseif strcmp(resultTableTimeSerangan.color{i}, 'red') || strcmp(resultTableTimeSerangan.color{i}, 'Malicious')
            figure(2);
            scatter(resultTableTimeSerangan.x(i), resultTableTimeSerangan.y(i), 64, 'r', 'filled'); 
        elseif strcmp(resultTableTimeSerangan.color{i}, 'magenta')
            figure(2);
            scatter(resultTableTimeSerangan.x(i), resultTableTimeSerangan.y(i), 64, 'm', 'filled');
        else
            figure(2);
            scatter(resultTableTimeSerangan.x(i), resultTableTimeSerangan.y(i), 64, 'b', 'filled');
        end

        % Plot connections between nodes
        for j = i + 1:size(resultTableTimeSerangan.koneksi, 2)
            if resultTableTimeSerangan.koneksi(i, j) == 1
                if strcmp(resultTableTimeSerangan.color{i}, 'red') || strcmp(resultTableTimeSerangan.color{j}, 'red')
                    % Draw connections involving red nodes with a different style
                    plot([resultTableTimeSerangan.x(i), resultTableTimeSerangan.x(j)], [resultTableTimeSerangan.y(i), resultTableTimeSerangan.y(j)], 'r--', 'LineWidth', 1.5);
                else
                    % Draw connections between non-red nodes
                    plot([resultTableTimeSerangan.x(i), resultTableTimeSerangan.x(j)], [resultTableTimeSerangan.y(i), resultTableTimeSerangan.y(j)], 'b--', 'LineWidth', 1);
                end
            end
        end
    end

    figure(2);
    hold on; 
    % h1 = scatter(NaN, NaN, 100, 'green', 'X', 'LineWidth', 1.5); 
    h2 = scatter(NaN, NaN, 64, 'blue', 'o', 'filled');
    h3 = scatter(NaN, NaN, 64, 'red', 'o', 'filled');
    leg2 = legend([h2, h3], 'Node Kendaraan', 'Malicious', 'Location', 'northeast');
    set(leg2, 'Box', 'on');
    hold off;

    % Perhitungan delay dan throughput pada detik t_idx untuk group.Result
    total_pt_1 = sum(group.Result{t_idx}.pt);
    total_rt_1 = sum(group.Result{t_idx}.rt);
    Delay1 = total_pt_1 / total_rt_1;
    
    % Perhitungan throughput pada detik t_idx untuk group.Result
    paket_diterima_1 = group.Result{t_idx}.rt; % paket data yang diterima dalam kb
    waktu_pengiriman_1 = group.Result{t_idx}.pt; % waktu pengiriman dalam detik
    Throughput1 = paket_diterima_1 ./ max(waktu_pengiriman_1, 1);
    
    % Perhitungan delay pada detik t_idx untuk group.ResultTime
    total_pt_2 = sum(group.ResultTime{t_idx}.pt);
    total_rt_2 = sum(group.ResultTime{t_idx}.rt);
    Delay2 = total_pt_2 / total_rt_2;
    
    % Perhitungan throughput pada detik t_idx untuk group.ResultTime
    paket_diterima_2 = group.ResultTime{t_idx}.rt; % paket data yang diterima dalam kb
    waktu_pengiriman_2 = group.ResultTime{t_idx}.pt; % waktu pengiriman dalam detik
    Throughput2 = paket_diterima_2 ./ max(waktu_pengiriman_2, 1);
    
    % Menyimpan hasil perhitungan delay dan throughput
    delay1(t_idx) = mean(Delay1);
    throughput1(t_idx) = mean(Throughput1); % Menggunakan mean untuk mendapatkan nilai rata-rata jika ada beberapa elemen
    delay2(t_idx) = mean(Delay2);
    throughput2(t_idx) = mean(Throughput2); % Menggunakan mean untuk mendapatkan nilai rata-rata jika ada beberapa elemen

    % Plot delay
    figure(3);
    plot(1:t_idx, delay1(1:t_idx), 'g.-'); % Plot delay dari figure 1
    hold on;
    plot(1:t_idx, delay2(1:t_idx), 'r.-');
    h_delay = legend('Normal', 'Under Attack', 'Location', 'northeast');
    set(h_delay, 'Box', 'on');
    hold off;

    % Plot throughput
    figure(4);
    plot(1:t_idx, throughput1(1:t_idx), 'g.-'); % Plot throughput dari figure 1
    hold on;
    plot(1:t_idx, throughput2(1:t_idx), 'r.-');
    h_throughput = legend('Normal', 'Under Attack', 'Location', 'northeast');
    set(h_throughput, 'Box', 'on');
    hold off;

    pause(0.0);

end

hold off;

% Tampilkan hasil rute
disp('AODV Route:');
disp(route);





%% --------------------------Trust Evaluation-----------------------------------------


numNodes = size(result, 1); 

trust_values = ones(numNodes, 1); 
trust_opinion_values = zeros(numNodes, 1); 
total_trust = zeros(numNodes, 1); 
trust_threshold = 1; 

% Loop through each node to evaluate trust
for i = 1:numNodes 
    % Update the status for the current node from the data
    status = result.Status(i); % Get the status of the current node (element-wise)
    
    % Trust evaluation based on behavior
    if strcmp(status, "Connected")
        trust_values(i) = trust_values(i) + 0.25; % Increase trust
    else
        trust_values(i) = trust_values(i) - 0.25; % Decrease trust
    end
    % Clamp trust values between 0 and max (e.g., 1 or any chosen upper limit)
    trust_values(i) = max(0, trust_values(i)); 
    
    % Now compute T_(V_BOi), sum of opinions from neighboring nodes
    % Simplified: assume some direct and indirect opinions
    direct_opinions = 0.9 * trust_values(i); % Example of direct opinion D_T_VOi
    indirect_opinions = 0.1 * trust_values(i); % Example of indirect opinion I_T_VOi
    trust_opinion_values(i) = direct_opinions + indirect_opinions;
    
    % Combine trust values as per formula T_(V_AB) = T_(V_B) + sum(T_(V_BOi))
    total_trust(i) = 0.8 * trust_values(i) + 0.2 * trust_opinion_values(i);
    
    % If trust is less than threshold, broadcast message about low trust
    if total_trust(i) <= trust_threshold
        fprintf('Node %d has low trust: %.2f\n', i, total_trust(i));
        % Broadcast logic can be added here
    end
end

% Menyimpan total_trust ke dalam kolom baru di result
result.total_trust = total_trust;




%% ------- Muat Data Untuk GSO -------------------------------------------------------
selectedData = result(1:313, :); 

%% ------- Inisialisasi variabel -----------------------------------------------------
s = 1;                          % step size
L0 = 5;                         % Luciferin awal
r0 = 300;                       % Jangkauan keputusan awal
rho = 0.4;                      % Konstanta peluruhan luciferin
gamma = 0.6;                    % Konstanta peningkatan luciferin
beta = 0.08;                    % Konstanta pembaruan jangkauan keputusan
rs = 300;                       % Jangkauan keputusan maksimum
nt = 5;                         % Ambang batas pembaruan jangkauan keputusan
maxIter = 150;                  % Jumlah iterasi maksimum
showplot = true;                % Flag untuk menampilkan plot
trust_threshold = 1;            % Ambang batas kepercayaan

%% ------- Muat Data Lane ------------------------------------------------------------
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

%% ------- Hitung Rata-rata Kecepatan ------------------------------------------------
selectedData.avgSpeed = zeros(height(selectedData), 1);
uniqueIds = unique(selectedData.id);

for i = 1:length(uniqueIds)
    currentId = uniqueIds{i};
    
    for currentTime = 0:20
        previousData = selectedData(strcmp(selectedData.id, currentId) & selectedData.t <= currentTime, :);
        
        if ~isempty(previousData)
            avgSpeed = mean(previousData.speed);
        else
            avgSpeed = 0;
        end
        
        selectedData.avgSpeed(strcmp(selectedData.id, currentId) & selectedData.t == currentTime) = avgSpeed;
    end
end

%% ------- Hitung Fitness Value ------------------------------------------------------
selectedData.fitnessValue = zeros(height(selectedData), 1);

for i = 1:height(selectedData)
    lane_length = selectedData.length(i);  
    speed = selectedData.speed(i);         
    trust = selectedData.total_trust(i);   
    
    if speed == 0
        fitnessValue = 0;
    else
        fitnessValue = (lane_length / speed) + trust;
    end
    
    selectedData.fitnessValue(i) = fitnessValue;
end

%% --------- LOOP UTAMA GSO ----------------------------------------------------------
uniqueTimes = unique(selectedData.t(selectedData.t <= 20));  % Filter data sampai time = 20

% Inisialisasi variabel untuk menyimpan evaluasi metrik
clusterHeadSelectionTime = zeros(length(uniqueTimes), 1);
totalPacketsGSO = zeros(length(uniqueTimes), 1);
deliveredPacketsGSO = zeros(length(uniqueTimes), 1);
delayGSO = zeros(length(uniqueTimes), 1);
energyConsumptionGSO = zeros(length(uniqueTimes), 1);
clusteringEfficiency = zeros(length(uniqueTimes), 1);

% Tambahan variabel untuk durasi
clusterHeadDuration = zeros(height(selectedData), 1);
clusterMemberDuration = zeros(height(selectedData), 1);

for timeIdx = 1:length(uniqueTimes)
    currentTime = uniqueTimes(timeIdx);
    
    % Pilih data untuk waktu saat ini
    currentAgents = selectedData(selectedData.t == currentTime, :);
    n = height(currentAgents); 
    
    Agent = [currentAgents.x, currentAgents.y];  % Posisi awal agen
    
    % Inisialisasi luciferin dan jangkauan keputusan
    luciferin = L0 * ones(n, 1);
    decision_range = r0 * ones(n, 1);

    % Simpan jumlah paket yang dihasilkan dan berhasil dikirim
    totalPacketsGSO(timeIdx) = n;
    deliveredPackets = 0;
    totalDelay = 0;
    
    % Menyimpan riwayat posisi cacing
    trail = zeros(n, 2, maxIter+1);
    trail(:, :, 1) = Agent;
    
    % Indeks untuk agen yang bergerak dan tidak bergerak
    moving_agents = false(n, 1);
    
    % Inisialisasi konsumsi energi
    totalEnergyConsumption = 0;
    
    % Mulai pengukuran waktu pemilihan Cluster Head
    tic;
    
    for t = 1:maxIter
        energyConsumedIteration = 0;
        
        % Langkah 1: Pembaruan luciferin
        for i = 1:n
            luciferin(i) = (1 - rho) * luciferin(i) + gamma * currentAgents.fitnessValue(i);
        end

        % Langkah 2: Tentukan tetangga
        neighbors = cell(n, 1);
        for i = 1:n
            neighbors{i} = [];
            for j = 1:n
                if i ~= j
                    distance = sqrt((Agent(i, 1) - Agent(j, 1))^2 + (Agent(i, 2) - Agent(j, 2))^2);
                    if distance < decision_range(i) && luciferin(j) > luciferin(i)
                        neighbors{i} = [neighbors{i}; j];
                    end
                end
            end
        end

        % Langkah 3: Perbarui posisi
        for i = 1:n
            if currentAgents.total_trust(i) >= trust_threshold
                moving_agents(i) = true;
                if ~isempty(neighbors{i})
                    prob = zeros(numel(neighbors{i}), 1);
                    sum_luciferin_diff = sum(luciferin(neighbors{i}) - luciferin(i));
                    for k = 1:numel(neighbors{i})
                        j = neighbors{i}(k);
                        prob(k) = (luciferin(j) - luciferin(i)) / sum_luciferin_diff;
                    end
                    
                    selected_neighbor = neighbors{i}(randsample(1:numel(neighbors{i}), 1, true, prob));
                    distance = norm(Agent(i, :) - Agent(selected_neighbor, :));
                    
                    Agent(i, 1) = Agent(i, 1) + s * (Agent(selected_neighbor, 1) - Agent(i, 1)) / distance;
                    Agent(i, 2) = Agent(i, 2) + s * (Agent(selected_neighbor, 2) - Agent(i, 2)) / distance;
                    
                    energyConsumedIteration = energyConsumedIteration + distance * 0.5;
                    packetDelay = distance / currentAgents.avgSpeed(i);
                    totalDelay = totalDelay + packetDelay;
                end
            else
                moving_agents(i) = false;
            end
        end

        totalEnergyConsumption = totalEnergyConsumption + energyConsumedIteration;
        
        % Langkah 4: Pembaruan jangkauan keputusan
        for i = 1:n
            decision_range(i) = min(rs, max(0, decision_range(i) + beta * (nt - numel(neighbors{i}))));
        end
        
        trail(:, :, t+1) = Agent;
    end
    
    deliveredPacketsGSO(timeIdx) = deliveredPackets;
    delayGSO(timeIdx) = totalDelay / max(1, deliveredPackets);
    
    % Perhitungan evaluasi metrik
    clusteredVehicles = sum(moving_agents);
    clusteringEfficiency(timeIdx) = (clusteredVehicles / n) * 100;
    
    % Update durasi
    clusterHeadDuration(moving_agents) = clusterHeadDuration(moving_agents) + toc;
    clusterMemberDuration(moving_agents) = clusterMemberDuration(moving_agents) + toc;
    
    % Proses pemilihan Cluster Head
    [~, clusterHeadIdx] = max(luciferin);
    clusterHeadSelectionTime(timeIdx) = toc;
    energyConsumptionGSO(timeIdx) = totalEnergyConsumption;
    
    % Plot posisi agen
    if showplot
        figure(5);
        clf;
        hold on;
        for i = 1:n
            if moving_agents(i)
                plot(squeeze(trail(i, 1, 1:t+1)), squeeze(trail(i, 2, 1:t+1)), '-');
            end
        end
        plot(Agent(moving_agents, 1), Agent(moving_agents, 2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        plot(Agent(~moving_agents, 1), Agent(~moving_agents, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        legend('Trust Node', 'Malicious Node');
        title(['Time ', num2str(currentTime), ' - Iterasi ke-', num2str(t)]);
        hold off;
        drawnow;
    end
end

% Rata-rata durasi
averageClusterHeadDuration = mean(clusterHeadDuration(clusterHeadDuration > 0));
averageClusterMemberDuration = mean(clusterMemberDuration(clusterMemberDuration > 0));
disp(['Average Cluster Head Duration: ', num2str(averageClusterHeadDuration), ' seconds']);
disp(['Average Cluster Member Duration: ', num2str(averageClusterMemberDuration), ' seconds']);


% Hitung PDR untuk GSO
PDR_GSO = sum(deliveredPacketsGSO) / sum(totalPacketsGSO) * 100;  % Dalam persen
disp(['PDR GSO: ', num2str(PDR_GSO), '%']);

%% ---------- Plot Delay untuk GSO ---------------------------------------------
% figure;
% plot(uniqueTimes, delayGSO, '-o', 'LineWidth', 2, 'DisplayName', 'Delay GSO');
% xlabel('Time');
% ylabel('Delay (s)');
% title('Average Delay per Time for GSO');
% legend('Location', 'northwest');
% grid on;

%% ---------- Plot Konsumsi Energi GSO ---------------------------------------------
% figure;
% plot(uniqueTimes, energyConsumptionGSO, '-o', 'LineWidth', 2, 'DisplayName', 'GSO');
% xlabel('Time');
% ylabel('Total Energy Consumption');
% title('Energy Consumption of GSO up to time = 20');
% legend('Location', 'northwest');
% grid on;







%% ----------Particle Swarm Optimization (Corrected)----------------------------------

% Load Data for PSO from Result
selectedData = result(1:313, :);  % Using data from result(1:313, :)

unique_times = unique(selectedData.t);  % List of unique times
energyConsumptionPSO = zeros(length(unique_times), 1);  % Energy consumption for PSO
totalPacketsPSO = zeros(length(unique_times), 1);  % Total packets generated
deliveredPacketsPSO = zeros(length(unique_times), 1);  % Packets successfully delivered
execution_times = zeros(length(unique_times), 1);  % Execution times

% Additional variables for durations
clusterHeadDurationPSO = zeros(height(selectedData), 1);  % Duration as Cluster Head
clusterMemberDurationPSO = zeros(height(selectedData), 1);  % Duration as Cluster Member

% Dangerous nodes are identified based on the "Disconnected" status
dangerousNodes = strcmp(selectedData.Status, 'Disconnected');

% Initialize figure for plotting
figure_handle = figure;

% Main loop for PSO
for t_idx = 1:length(unique_times)
    current_time = unique_times(t_idx);
    current_data = selectedData(selectedData.t == current_time, :);
    X = current_data.x;
    Y = current_data.y;
    Type = current_data.type;  % Vehicle type (car/taxi)
    Speed = current_data.speed;
    Lane = current_data.lane;

    % Convert X and Y columns to numeric if necessary
    if iscell(X)
        X = str2double(X);
    end
    if iscell(Y)
        Y = str2double(Y);
    end
    if iscell(Speed)
        Speed = str2double(Speed);
    end
    if iscell(Lane)
        Lane = str2double(Lane);
    end

    % Convert 'Type' to numeric (car = 1, taxi = 2)
    type_numeric = zeros(size(Type));
    type_numeric(strcmp(Type, 'mobil')) = 1;
    type_numeric(strcmp(Type, 'taxi')) = 2;

    % Combine features into one matrix (X, Y, type_numeric, Speed)
    features = [X Y type_numeric Speed];

    nNode = length(X);  % Number of nodes (vehicles)
    nCH = 2;  % Number of Cluster Heads to be selected

    % Track packets and energy consumption
    totalPacketsPSO(t_idx) = nNode;
    deliveredPackets = 0;

    % PSO parameters
    VarSize = [1 nCH];
    VarMin = 1;
    VarMax = nNode;

    CostFunction = @(ch) ClusterCostPSO(ch, features, nNode, nCH);

    MaxIt = 100;  % Maximum iterations
    energyConsumed = 0;
    nPop = 50;  % Population size

    w = 1;  % Inertia weight
    wdamp = 0.99;  % Damping ratio
    c1 = 1.5;  % Personal learning coefficient
    c2 = 2.0;  % Global learning coefficient

    VelMax = 0.1 * (VarMax - VarMin);
    VelMin = -VelMax;

    % Initialize particles
    empty_particle.Position = [];
    empty_particle.Cost = [];
    empty_particle.Velocity = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    particle = repmat(empty_particle, nPop, 1);
    GlobalBest.Cost = inf;

    % Start the timer for Cluster Head selection
    selection_start_time = tic;  % Use tic to start timing

    % Initialize particle positions and costs
    for i = 1:nPop
        particle(i).Position = randi([VarMin VarMax], VarSize);
        particle(i).Velocity = zeros(VarSize);
        particle(i).Cost = CostFunction(particle(i).Position);
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end
    end

    BestCost = zeros(MaxIt, 1);

    % PSO main loop
    for it = 1:MaxIt
        for i = 1:nPop
            particle(i).Velocity = w * particle(i).Velocity ...
                + c1 * rand(VarSize) .* (particle(i).Best.Position - particle(i).Position) ...
                + c2 * rand(VarSize) .* (GlobalBest.Position - particle(i).Position);

            particle(i).Velocity = max(particle(i).Velocity, VelMin);
            particle(i).Velocity = min(particle(i).Velocity, VelMax);

            particle(i).Position = particle(i).Position + particle(i).Velocity;
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);
            particle(i).Position = round(particle(i).Position);
            particle(i).Cost = CostFunction(particle(i).Position);

            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end
            end

            % Simulate particle movement and energy consumption
            distanceMoved = rand * 10;
            energyConsumed = energyConsumed + distanceMoved * 0.3;

            % Simulate data delivery
            if ~dangerousNodes(i)
                if rand() > 0.1
                    deliveredPackets = deliveredPackets + 1;
                end
            else
                disp(['Node ', num2str(i), ' is dangerous and does not forward packets.']);
            end
        end

        BestCost(it) = GlobalBest.Cost;
    end

    % Calculate execution time using toc
    duration = toc(selection_start_time);  % Get elapsed time
    execution_times(t_idx) = duration;

    % Determine Cluster Head and Member roles
    currentCH = GlobalBest.Position;
    isClusterHead = false(nNode, 1);
    isClusterHead(currentCH) = true;

    % Update durations incrementally
    clusterHeadDurationPSO(isClusterHead) = clusterHeadDurationPSO(isClusterHead) + duration;
    clusterMemberDurationPSO(~isClusterHead) = clusterMemberDurationPSO(~isClusterHead) + duration;

    % Display Cluster Head selection info
    disp(['Time: ' num2str(current_time) ', Best Cluster Head Indices: ' num2str(currentCH)]);
    disp(['Time: ' num2str(current_time) ', Best Cost: ' num2str(GlobalBest.Cost)]);

    % Plot vehicle positions
    clf(figure_handle);
    hold on;
    plot(X, Y, 'bo');
    plot(X(currentCH), Y(currentCH), 'r*', 'MarkerSize', 10);
    title(['VANET Clustering with PSO (Time = ' num2str(current_time) ')']);
    xlabel('X Position');
    ylabel('Y Position');
    legend('Vehicles', 'Cluster Heads');
    grid on;
    pause(1);
end



% Calculate and display average durations
averageClusterHeadDurationPSO = mean(clusterHeadDurationPSO(clusterHeadDurationPSO > 0));
averageClusterMemberDurationPSO = mean(clusterMemberDurationPSO(clusterMemberDurationPSO > 0));
disp(['Average Cluster Head Duration PSO: ', num2str(averageClusterHeadDurationPSO), ' seconds']);
disp(['Average Cluster Member Duration PSO: ', num2str(averageClusterMemberDurationPSO), ' seconds']);


% Hitung PDR untuk PSO
PDR_PSO = sum(deliveredPacketsPSO) / sum(totalPacketsPSO) * 100;
disp(['PDR PSO: ', num2str(PDR_PSO), '%']);



% Plot comparison of Average Cluster Head Duration and Member Duration
figure;
hold on;
bar(1, averageClusterHeadDuration, 'b');
bar(2, averageClusterHeadDurationPSO, 'r');
% bar(3, averageClusterMemberDuration, 'r');
% bar(4, averageClusterMemberDurationPSO, 'b');
xlabel('Cluster Head');
ylabel('Duration (seconds)');
title('Comparison of Cluster Head Durations: GSO vs PSO');
legend('GSO', 'PSO');
grid on;
hold off;


%% ----------- Plot Perbandingan PDR ---------------------------------------

% figure;
% bar([PDR_GSO, PDR_PSO]);
% set(gca, 'XTickLabel', {'GSO', 'PSO'});
% ylabel('Packet Delivery Ratio (PDR) %');
% title('Comparison of PDR: GSO vs PSO');
% grid on;

%% ---------- Plot Perbandingan Konsumsi Energi -------------------------------------
% figure;
% % plot(uniqueTimes, energyConsumptionGSO, '-o', 'LineWidth', 2, 'DisplayName', 'GSO');
% % plot(1:maxIter, energyConsumptionGSO, '-o', 'LineWidth', 2, 'DisplayName', 'GSO');
% hold on;
% plot(1:length(unique_times), energyConsumptionPSO, '-s', 'LineWidth', 2, 'DisplayName', 'PSO');
% xlabel('Iteration');
% ylabel('Energy Consumption');
% % title('Comparison of Energy Consumption: GSO vs PSO');
% legend('Location', 'northwest');
% grid on;
% hold off;

%% ----------- Plot waktu pemilihan Cluster Head untuk GSO -------------------------
figure;
plot(uniqueTimes, clusterHeadSelectionTime, '-o', 'LineWidth', 2, 'DisplayName', 'GSO Selection Time');
hold on;

% Plot waktu eksekusi PSO untuk setiap waktu
plot(unique_times, execution_times, '-s', 'LineWidth', 2, 'DisplayName', 'PSO Execution Time');

% Tambahkan judul, label sumbu, dan legenda
title('Comparison of Cluster Head Selection Time: GSO vs PSO');
xlabel('Iteration');
ylabel('CH Selection Time');

% Ubah penempatan legenda ke pojok kiri atas
legend('Location', 'northwest');

grid on;
hold off;




%% Fungsi Cost untuk clustering
function cost = ClusterCostPSO(ch, features, nNode, nCH)
    % Menghitung total jarak dari setiap node ke Cluster Head terdekat
    totalDistance = 0;
    
    for i = 1:nNode
        minDist = inf;
        for j = 1:nCH
            % Hitung jarak Euclidean antar node dan CH
            dist = sqrt(sum((features(i,1:2) - features(ch(j),1:2)).^2)); % Jarak 2D (X, Y)
            if dist < minDist
                minDist = dist;
            end
        end
        % Tambah jarak terdekat ke total jarak
        totalDistance = totalDistance + minDist;
    end
    
    % Average distance as the cost
    cost = totalDistance / nNode;
end







