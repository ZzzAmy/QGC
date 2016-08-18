function [ matGT ] = LoadGroundTruth_Twitter( graph_path, N )

    global mapUID2ID
    
    fid = fopen(graph_path, 'r');
    
    cellData = {};
    
    matGT = [];
    
    while ~feof(fid)
        tline=fgetl(fid);
        t = str2num(tline);
        vecX = t(2:length(t));
        for i = 1:length(vecX)
            vecX(i) = mapUID2ID.get(vecX(i));
        end
        vecY = ones(length(t)-1, 1) * 1;
        vecVal = ones(length(t)-1, 1);
        
        vecLabel = sparse(vecX, vecY, vecVal, N, 1);
        
        matGT = [matGT vecLabel];
    end
end

