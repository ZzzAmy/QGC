function [ n, e ] = LoadGraph_DBLP(graph_path, vetexName_path)

    global matG
    global cellVertexNames
    
    fileTags = fopen(vetexName_path);
    cellVertexNames = textscan(fileTags, '%d%s', 'delimiter', '\t');
    
    fid = fopen(graph_path, 'r');
    cellData = textscan(fid, '%f%f%f', 'delimiter', '\t');
    fclose(fid);

    n = max(max(cellData{1, 1}), max(cellData{1, 2}));
    e = length(cellData{1, 1});
    
    % -------------- Construct MatG --------------- %
    
    X = cellData{1, 1};
    Y = cellData{1, 2};
    V = cellData{1, 3};
    
    matG = sparse(X, Y, V, n, n);
    matG = (matG' + matG);
end

