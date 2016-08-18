function [ matGT ] = LoadGroundTruth_FB( graph_path, N )

    fid = fopen(graph_path, 'r');
    
    cellData = {};
    
    matGT = [];
    
    while ~feof(fid)
        tline=fgetl(fid);
        t = str2num(tline(7:length(tline)));
        vecX = t(2:length(t)) + 1;
        vecY = ones(length(t)-1, 1) * 1;
        vecVal = ones(length(t)-1, 1);
        
        vecLabel = sparse(vecX, vecY, vecVal, N, 1);
        
        matGT = [matGT vecLabel];
    end
end

