function [weight_eigvec, vecQ, K] = QOGC_QGC(query, K, H, vecRel, MyLancType, threshold, eta)
    % Function: QOGC_QGC
    % This function performs Query-oriented Graph Clustering.
    %--------------------------------------------------------------------------------------
    % query: queryID
    % K: the number of clustering. (K = 0 ===> auto-determine the number of clusters K)
    % H: the top-H eigenvectors are obtained for clustering
    % vecRel: relevance vector
    % MyLancType: approach selection index
    % eta: is used only when MyLancType = 2, often set to a huge value
    %
    
    global matG
    global matG_SC
    global N    
    global maxitr
    
    %% Initial setting
    max_clustering_itr = 200;
    
    zeta = 0.05;
    q_size = 1;

    matR = spdiags(vecRel(:), 0, N, N);
    vecD = sum(matG);
    matD = spdiags(vecD(:), 0, N, N);
    
    vecRD = sum(matR * (matG * matR))';
    
    vecRD_sqrt =  1./ (vecRD .^ 0.5);
    matRD_sqrt = spdiags(vecRD_sqrt(:), 0, N, N);
    
    matG_SC = matR * matG * matR;
    matG_SC = matRD_sqrt * matG_SC * matRD_sqrt;       
    matG_SC = 0.5 * (matG_SC + matG_SC');
 
    %% Query setting
    if query > 0
        q_size = 1;
        matG_SC(query,:) = [];
        matG_SC(:,query) = [];
        rel = vecRel;
        rel(query,:) = [];
    else
        q_size = 0;
        rel = vecRel;
    end
    
    %% Run Graph Clustering
    clustering_time = 0;
    kmeansRun = true;
    while (kmeansRun && clustering_time < max_clustering_itr)
        clustering_time = clustering_time + 1;
        
        %% Eigenmap for dimensionality reduction
        if MyLancType == 0
            %
            % Minimize the objective function without considering the balance constraint
            %
            opts.issym = 1;
            opts.maxit = maxitr;
            [eigVec, eigVal] = eigs(matG_SC, H, 'la', opts);
            weight_eigvec = eigVec;
            fprintf('finish eigenmaps => QGC \n');   
        elseif MyLancType == 1
            %
            % Minimize the objective function by using Langrange multiplier
            %
            rel = rel / sqrt(rel' * rel);
            opts.issym = 1;
            opts.maxit = maxitr;
            [eigVec, eigVal] = eigs(@myGG, N - q_size, H, 'la', opts, rel);
            weight_eigvec = eigVec - rel * (rel' * eigVec);
            fprintf('finish eigenmaps => QGCB 1 \n');
        elseif MyLancType == 2
            %
            % Minimize the objective function by setting huge multiplier of the Lagrangian term.
            %
            rel = rel / sqrt(rel' * rel);
            opts.issym = 1;
            opts.maxit = maxitr;
            [eigVec, eigVal] = eigs(@myGG2, N - q_size, H, 'la', opts, rel, eta);
            weight_eigvec = eigVec;
            fprintf('finish eigenmaps => QGCB 2 \n');
        else
        end

%         vecDeig = diag(eigVal);
%         gradDeig = vecDeig(1:(length(vecDeig)-1)) - vecDeig(2:length(vecDeig));
%         mean_gradDeig = mean(gradDeig);
%         flag = 0;
%         for i = 1:length(gradDeig)
%             if gradDeig(i) > mean_gradDeig && flag == 1
%                 H = i + 1;
%                 break;
%             end
%             if gradDeig(i) > mean_gradDeig
%                 flag = flag + 1;
%             end
%         end
        
        %% Label assignment        
        weight_eigvec = weight_eigvec * (eigVal+10e-20);
%         weight_eigvec = weight_eigvec(:, 1:H);
        [oPhi, K] = VectorClustering(weight_eigvec, K, threshold);
        %[oPhi, K] = VectorClustering2(weight_eigvec, K, 0.05);
        sizePhi = sum(oPhi);

        if any(sizePhi(:)==0) == 0
            kmeansRun = false;
        end
    end
               
    %% Insert the query into the clustering result
    if query > 0
        vecQ = insertrows(oPhi, zeros(1, K), query-1);
    else
        vecQ = oPhi;
    end
end