function [ oPhi, K ] = VectorClustering(weight_eigvec, K, threshold)
    % Function: VectorClustering
    % This function is to cluster obtained eigenvectors in the eigenspace.
    % Notice that the similarity is defined as Cosine Similarity.
    % The scale of a vector represents its 'strength'.
    %--------------------------------------------------------------------------------------
    % weight_eigvec: weighted eigenvectors
    % K: the number of clustering. (K = 0 ===> auto-determine the number of clusters K)
    %
    
    p = 10;
    topP = 3;
    centroid_types = {'mean', 'power mean', 'top-P', 'ndist'};
    cent_type = centroid_types{3};
    
    H = size(weight_eigvec, 2);

    scale_weight_eigvec = zeros(size(weight_eigvec,1), 1);
    for i = 1:size(weight_eigvec,1)
        scale_weight_eigvec(i,1) = weight_eigvec(i,:) * weight_eigvec(i,:)';
    end

    %% Initialize the clustering tree leaves
    %oPhi = weight_eigvec .* repmat(sign(sum(weight_eigvec)),[N-q_size 1]);
    oPhi = weight_eigvec;
    matPartitionLabel = [];
    while nnz(oPhi) > 0
        %i = 1;
        %while nnz(oPhi(i,:)) == 0
        %    i = i + 1;
        %end
        i = find(sum(oPhi, 2), 1, 'first');

        matMatch = bsxfun(@times, weight_eigvec, weight_eigvec(i,:));
        matMatch = matMatch >= 0;
        vecMatch = matMatch * ones(H, 1);
        vecMatch = vecMatch >= H;
        matPartitionLabel = [matPartitionLabel vecMatch];

        vecZ = vecMatch == 0;
        oPhi = bsxfun(@times, oPhi, vecZ);
    end

    %% Initialize the similarity matrix
    G = size(matPartitionLabel,2);
    cent_norm = zeros(G, 1);
    centroids = zeros(G, H);
    variance = zeros(G*H, H);
    for i = 1:G
        [x, y] = find(matPartitionLabel(:,i));
        
        if strcmp(cent_type, 'mean')
            % mean of the all nodes to be centroid
            centroid4 = mean(weight_eigvec(x, :));
            centroid4 = centroid4 / norm(centroid4);
            centroids(i, :) = centroid / sqrt(centroid * centroid');
        elseif strcmp(cent_type, 'power mean')
            % quadratic mean of the all nodes to be centroid
            we_square = (weight_eigvec(x, :) .^ p) .* sign(weight_eigvec(x, :));
            centroid = (sum(we_square) / length(x));
            centroid = abs(centroid) .^ (1/p) .* sign(centroid);
            centroid = centroid / norm(centroid);
            centroids(i, :) = centroid;
        elseif strcmp(cent_type, 'top-P')
            % mean of the top-P nodes to be centroid
            
            if size(scale_weight_eigvec(x, :),1) > topP
                [v1, x1] = maxk(scale_weight_eigvec(x, :), topP);
            else
                [v1, x1] = maxk(scale_weight_eigvec(x, :), size(scale_weight_eigvec(x, :),1));
            end
            centroid2 = mean(weight_eigvec(x(x1,1), :));
            cent_norm(i,1) = norm(centroid2);
            centroid2 = centroid2 / norm(centroid2);
            centroids(i, :) = centroid2;
            
        elseif strcmp(cent_type, 'ndist')
            % mean of the all nodes to be centroid
            [vecMean, matVar] = MultiVarNormalDist(weight_eigvec(x, :));
            centroids(i, :) = vecMean;
            variance((H*(i-1)+1):(H*i), :) = matVar + 10e-5 * eye(H);
        else
        end
    end
    
    G = size(centroids, 1);
    
    matClusterSim = zeros(G);
    if strcmp(cent_type, 'ndist')
        % model-based HAC (by Gaussian prior)
        for i = 1:(G-1)
            for j = (i+1):G
                tmp = JensenDistOfNormal(centroids(i,:), variance((H*(i-1)+1):(H*i),:), centroids(j,:), variance((H*(j-1)+1):(H*j),:));
                matClusterSim(i,j) = tmp / (nnz(matPartitionLabel(:,i)) * nnz(matPartitionLabel(:,j)));
            end
        end
        matClusterSim = matClusterSim + matClusterSim';
        matClusterSim = matClusterSim + 1;
        matClusterSim = 1 ./ matClusterSim;
    else
        % non-model prior
        matClusterSim = centroids * centroids';
        scale_cluster = sqrt(diag(matClusterSim));
        matClusterSim = bsxfun(@rdivide, matClusterSim, scale_cluster);
        matClusterSim = bsxfun(@rdivide, matClusterSim, scale_cluster');
    end
    matClusterSim = matClusterSim - eye(G);
    X = length(matClusterSim);
    
    %% Hierarchical agglomerative clustering greedily
    %
    % Facebook data: matClusterSim >= -0.02 
    % Twitter data: matClusterSim >= 0.9
    %
    while (X > K && K > 0) || (nnz(matClusterSim >= threshold) > 0 && K == 0)
        %% Choose the most similar eigenvectors and merge them
        [ max_v, max_x, max_y ] = MaxOfMatrix( matClusterSim );
        vecMerge = (matPartitionLabel(:,max_x) + matPartitionLabel(:,max_y)) > 0;     
        if max_x > max_y
            matPartitionLabel(:,max_x) = [];
            matPartitionLabel(:,max_y) = [];
        elseif max_x < max_y
            matPartitionLabel(:,max_y) = [];
            matPartitionLabel(:,max_x) = [];
        else
            fprintf('max_x: %f , ', max_x);
            fprintf('max_y: %f , ', max_y);
            fprintf('max_v: %f , \n', max_v);
        end
        matPartitionLabel = [matPartitionLabel vecMerge];

        %% Update the similarity matrix
        G = size(matPartitionLabel,2);
        centroids = zeros(G, H);
        cent_norm = zeros(G, 1);
        for i = 1:G
            [x, y] = find(matPartitionLabel(:,i));
            
            if strcmp(cent_type, 'mean')
                % mean of the all nodes to be centroid
                centroid4 = mean(weight_eigvec(x, :));
                centroid4 = centroid4 / norm(centroid4);
                centroids(i, :) = centroid / sqrt(centroid * centroid');
            elseif strcmp(cent_type, 'power mean')
                % quadratic mean of the all nodes to be centroid
                we_square = (weight_eigvec(x, :) .^ p) .* sign(weight_eigvec(x, :));
                centroid = (sum(we_square) / length(x));
                centroid = abs(centroid) .^ (1/p) .* sign(centroid);
                centroid = centroid / norm(centroid);
                centroids(i, :) = centroid;
            elseif strcmp(cent_type, 'top-P')
                % mean of the top-P nodes to be centroid
                [v1, x1] = maxk(scale_weight_eigvec(x, :), topP);
                centroid2 = mean(weight_eigvec(x(x1,1), :));
                cent_norm(i, 1) = norm(centroid2);
                centroid2 = centroid2 / norm(centroid2);
                centroids(i, :) = centroid2;
            elseif strcmp(cent_type, 'ndist')
                % mean of the all nodes to be centroid
                [vecMean, matVar] = MultiVarNormalDist(weight_eigvec(x, :));
                centroids(i, :) = vecMean;
                variance((H*(i-1)+1):(H*i), :) = matVar + 10e-5 * eye(H);
            else
            end
        end
        matClusterSim = zeros(G);
        if strcmp(cent_type, 'ndist')
            % model-based HAC (by Gaussian prior)
            for i = 1:(G-1)
                for j = (i+1):G
                    tmp = JensenDistOfNormal(centroids(i,:), variance((H*(i-1)+1):(H*i),:), centroids(j,:), variance((H*(j-1)+1):(H*j),:));
                    matClusterSim(i,j) = tmp / (nnz(matPartitionLabel(:,i)) * nnz(matPartitionLabel(:,j)));
                end
            end
            matClusterSim = matClusterSim + matClusterSim';
            matClusterSim = matClusterSim + 1;
            matClusterSim = 1 ./ matClusterSim;
        else
            % non-model prior
            matClusterSim = centroids * centroids';
            scale_cluster = sqrt(diag(matClusterSim));
            matClusterSim = bsxfun(@rdivide, matClusterSim, scale_cluster);
            matClusterSim = bsxfun(@rdivide, matClusterSim, scale_cluster');
            matClusterSim = matClusterSim - 10*eye(G);
        end
        matClusterSim = matClusterSim - eye(G);
        X = length(matClusterSim);
    end

    %% Determine which cluster a data belongs to
    %
    % Kmeans
    %
    %[a, c] = kmeans(weight_eigvec, K, 'emptyaction', 'drop');
    
    %
    % HAC, and then inner product of data and centroid
    
    index = weight_eigvec * centroids';
    scale_centroid = diag(centroids * centroids');
    index = bsxfun(@rdivide, index, sqrt(scale_weight_eigvec));
    index = bsxfun(@rdivide, index, sqrt(scale_centroid)');
    [x, y] = max(index, [], 2);
    n = length(y);
    oPhi = sparse(1:n, y, ones(n,1), n, size(centroids, 1));
    sum_oPhi = sum(oPhi);
    for i = 1:length(sum_oPhi)
        if sum_oPhi(length(sum_oPhi)+1-i) == 0
            oPhi(:,length(sum_oPhi)+1-i) = [];
        end
    end
    K = size(oPhi, 2);
    %
    % HAC
    %
    %oPhi = matPartitionLabel;
    %K = size(matPartitionLabel, 2);
end

