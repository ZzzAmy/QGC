function [ accuracy, precision, recall, f1score, running_time ] = LoadSCDresults( query, SCD_path, data_flag)
    
    global matG
    global N
    
    hop_length = 1;
    
    %% Load the groubd truth data
    if strcmp(data_flag, 'facebook')
        % OSX file path
        [ X_grnd ] = LoadGroundTruth_FB(strcat('/Users/iankuoli/Dataset/sna_facebook/facebook/', int2str(query), '.circles'), N);
        % Windows file path
        %[ matGT ] = LoadGroundTruth_FB(strcat('C:/Dataset/sna_facebook/facebook/', int2str(query-1), '.circles'), N);
    elseif strcmp(data_flag, 'twitter')
        % OSX file path
        [ X_grnd ] = LoadGroundTruth_Twitter(strcat('/Users/iankuoli/Dataset/sna_twitter/twitter/', int2str(mapID2UID.get(queries(i))), '.circles'), N);
        % Windows file path
        %[ matGT ] = LoadGroundTruth_Twitter(strcat('C:/Dataset/sna_twitter/twitter/', int2str(queries_UID(i)), '.circles'), N);
    else
    end
    
    
    % Set the consider
    vecConsider = zeros(N, 1);
    vecConsider(query+1, 1) = 1;
    for p = 1:hop_length
        vecConsider = vecConsider + matG * vecConsider;
        vecConsider = vecConsider > 0;
    end
    vecConsider(query+1, 1) = 0;
    
    
    fid = fopen(SCD_path, 'r');
    cellData = textscan(fid, '%s%s', 'delimiter', '=');
    fclose(fid);
    
    str_cluster = cellData{1,2}{6};
    C = regexp(str_cluster,'],[','split');
    C{1,1} = C{1,1}(1,3:end);
    C{1,end} = C{1,end}(1,1:end-2);
    
    indx_x = [];
    indx_y = [];
    indx_v = [];
    
    for k = 1:length(C)
        cc = textscan(C{1,k},'%f', 'delimiter', ',');
        indx_x = [indx_x cc{1,1}'];
        indx_y = [indx_y (k*ones(1,length(cc{1,1})))];
    end
    indx_v = ones(1, length(indx_x));
    indx_x = indx_x + 1;
    
    X1 = sparse(indx_x, indx_y, indx_v, N, length(C));
    
    running_time = str2num(cellData{1,2}{9});
    

    [accuracy, precision, recall, f1score] = MeasureF1Score(X1, X_grnd, vecConsider, 'f1');
end

