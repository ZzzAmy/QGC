function [ n, e, mapUID2ID, mapID2UID] = LoadGraph_Twitter(graph_path)

    global matG
    global cellVertexNames
    
    fid = fopen(graph_path, 'r');
    cellData = textscan(fid, '%f%f', 'delimiter', '\t');
    fclose(fid);
    
    % -------------- Construct matG --------------- %
    
    X = cellData{1, 1};
    Y = cellData{1, 2};
    
    e = size(X, 1);
    
    Xid = zeros(e,1);
    Yid = zeros(e,1);
    
    id = 1;
    mapUID2ID = java.util.HashMap;
    mapID2UID = java.util.HashMap;
    for i = 1:size(X,1)
        if mod(i, 10000) == 0
            fprintf('i = %f \n', i);
        end
        uid1 = X(i, 1);
        uid2 = Y(i, 1);
        if (mapUID2ID.containsKey(uid1) == false)
            mapUID2ID.put(uid1, id);
            mapID2UID.put(id, uid1);
            id = id + 1;
        end
        if (mapUID2ID.containsKey(uid2) == false)
            mapUID2ID.put(uid2, id);
            mapID2UID.put(id, uid2);
            id = id + 1;
        end
        X_id(i,1) = mapUID2ID.get(uid1);
        Y_id(i,1) = mapUID2ID.get(uid2);
    end
    n = id - 1;
    
    matG = sparse(X_id, Y_id, ones(e,1), n, n);
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

