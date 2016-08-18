function [ n ] = LoadData( file_path, t )

    global matG
    global matA
    global label
    global cellData
    
    fid = fopen(file_path, 'r');
    cellData = textscan(fid, '%f%f%f', 'delimiter', '\t');
    fclose(fid);

    n = length(cellData{1, 1});
    
    % -------------- Construct MatG --------------- %
    
    % calculate t-nearest neighbors
    mat_distance_l2 = zeros(n, t);
    mat_distance_idx = zeros(n, t);
    
    for i = 1:n
        distance = zeros(n, 1) + 10e6;
        for j = 1:n
            if i ~= j
                x_i = [cellData{1, 1}(i ,1), cellData{1, 2}(i ,1)];
                x_j = [cellData{1, 1}(j ,1), cellData{1, 2}(j ,1)];
                distance(j ,1) = (x_i - x_j) * (x_i - x_j)';
            end
        end
        
        for k = 1:t
            [val, indx] = min(distance);
            mat_distance_l2(i, k) = val;
            mat_distance_idx(i, k) = indx;
            distance(indx, 1) = 10e6;
        end
    end
    
    matG = sparse(n ,n);
    
    for i = 1:n
        sigma_i = mat_distance_l2(i, t);
        for k = 1:t
            j = mat_distance_idx(i, k);
            sigma_j = mat_distance_l2(j, t);
            val = exp(- mat_distance_l2(i, k) / (2 * sqrt(sigma_i * sigma_j)));
            
            matG(i, j) = val;
            matG(j, i) = val;
        end
    end
    
    matG(n, n) = 0;
    
%     vecX = zeros(2 * n * t);
%     vecY = zeros(2 * n * t);
%     vecVal = zeros(2 * n * t);
%     count = 1;
%     for i = 1:n
%         sigma_i = mat_distance_l2(i, t);
%         for k = 1:t
%             j = mat_distance_idx(i, k);
%             sigma_j = mat_distance_l2(j, t);
%             val = exp(- mat_distance_l2(i, k) / (2 * sigma_i * sigma_j));
%             
%             vecX(count) = i;
%             vecY(count) = j;
%             vecVal(count) = val;
%             
%             vecX(count + 1) = i;
%             vecY(count + 1) = j;
%             vecVal(count + 1) = val;
%             
%             count = count+2;
%         end
%     end
            

end

