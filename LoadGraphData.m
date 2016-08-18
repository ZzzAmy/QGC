function [ n, e ] = LoadGraphData(file_path)

    global matG
    global cellData
    
    fid = fopen(file_path, 'r');
    cellData = textscan(fid, '%f%f', 'delimiter', '\t');
    fclose(fid);

    n = max(max(cellData{1, 1}), max(cellData{1, 2}));
    e = length(cellData{1, 1});
    
    % -------------- Construct MatG --------------- %
    
    X = cellData{1, 1};
    Y = cellData{1, 2};
    V = ones(e, 1);
    
    matG = sparse(X, Y, V, n, n);
    matG = (matG' + matG);
end

