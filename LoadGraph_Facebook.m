function [ n, e ] = LoadGraph_Facebook(graph_path)

    global matG
    global cellVertexNames
    
    fid = fopen(graph_path, 'r');
    cellData = textscan(fid, '%f%f', 'delimiter', '\t');
    fclose(fid);

    n = max(max(cellData{1, 1}), max(cellData{1, 2})) + 1;
    e = length(cellData{1, 1});
    
    % -------------- Construct matG --------------- %
    
    X = cellData{1, 1};
    Y = cellData{1, 2};
    
    matG = sparse(X+1, Y+1, ones(e,1), n, n);
    matG = 0.5 * (matG' + matG);
    
    
    % -------------- Construct cellVertexNames --------------- %
    cellVertexNames = cell(1,2);
    cellVertexNames{1,1} = zeros(n,1);
    cellVertexNames{1,2} = cell(n,1);
    for i = 1:n
        cellVertexNames{1,1}(i,1) = i;
        cellVertexNames{1,2}{i,1} = i-1;
    end
end

