function [ NCut_rel, vecRankScore, vecPerform ] = QOCut( query, weight_eigvec, vecQ, type, vecRel, rank_type, M, m_type)
    % Function: QOGC_QGC
    % This function performs Query-oriented Graph Clustering.
    %--------------------------------------------------------------------------------------
    % query: queryID
    % weight_eigvec: weighted eigenvectors
    % vecQ: label assignment
    % type: 'ratiocut' or 'ncut'
    % vecRel: relevance vector
    % rank_type: choose a method to rank data for each cluster
    %            'rel', 'eig', 'eigXrel', 'eigXrel2', 'eig+rel'
    % M: QONCut@M => only consider top-M items for each cluster
    % m_type: 'TopNperCluster' or 'TotalN'
    %
    
    global matG

    global N
    K = size(vecQ, 2);
    
    %% Label assigning
    %{
    vecLabel = zeros(N, K);
    for i = 1:N
        [val, index] = max(vecQ(i,:));
        vecLabel(i, index) = 1;
        
        if mod(i, 10000) == 0
            fprintf('%d\n', i);
        end
    end
    vecLabel(query, :) = 0;
    %}
    
    [vecVal, vecIndex] = max(vecQ, [], 2);
    vecLabel = sparse(1:N, vecIndex, ones(N, 1), N, K);

    %% Ranking score is according to the relevance
    vecPerform = bsxfun(@times, vecLabel, vecRel);

    if query > 0
        vecPerform(query, :) = zeros(1, K);
    end
    
    if nargin > 6
        vecLabel = zeros(N, K);
        vecRel2 = zeros(N, 1);
            
        if strcmp(m_type, 'TopNperCluster')
            % top-M for each cluster
            for i = 1:K
                [topM, loc] = maxk(vecPerform(:, i), M);
                for j = 1:M
                    indx = find(vecPerform(:, i) == topM(j), 1);
                    vecLabel(indx, i) = vecPerform(indx, i);
                    vecPerform(indx, i) = -100;
                    vecRel2(indx) = vecRel(indx);
                end
            end
        elseif strcmp(m_type, 'TotalN')
            % top-M of the whole graph
            [topM, loc] = maxk(vecRel, M);
            mask = sparse(loc, ones(M,1), ones(M,1), N, 1);
            vecRel2 = mask .* vecRel;
            vecLabel = bsxfun(@times, vecPerform, mask);
        else
        end
        
        vecPerform = vecLabel;
        
        % Only consider the top M vertices, left others with relevance 0s.
        vecRel = vecRel2;
    end 
        
    %% Calculate the normalization term
    if strcmp(type, 'ratiocut')
        % the idea is the same with RadioCut
        vecPerform_sum = ones(1, N) * vecPerform;
        
        % ..... waiting for completeness
    else
        % the idea is the same with NCut
        vecPerform_sum = ones(1, N) * (repmat(vecRel, [1, K]) .* (matG * vecPerform));

        % Normalization    
        vecPerform_norm = vecPerform ./ repmat(vecPerform_sum, [N, 1]);

        % Calculate Cut
        vecTemp = full(repmat(vecRel,[1 K]) .* (matG * vecPerform_norm) - 10e5 * vecPerform_norm);
        vecTemp(vecTemp<0) = 0;
        NCut_rel = sum(sum(vecTemp));
    end
    
    %% Ranking for each cluster
    eigIndex = insertrows(weight_eigvec, zeros(1,size(weight_eigvec, 2)), query-1);
    
    if strcmp(rank_type, 'rel')
        %
        % according to the relevance
        %
        vecRankScore = vecPerform;
    elseif strcmp(rank_type, 'eig')
        %
        % according to the eigenvector
        %
        vecRankScore = zeros(N, K);
        for i = 1:K
            [x, y, z] = find(vecQ(:,i));
            if length(x) > 1 
                mean_vec = mean(eigIndex(x,:));
            else
                mean_vec = eigIndex(x,:);
            end
            mean_vec = mean(eigIndex(x,:));
            tmp = eigIndex * mean_vec';
            vecRankScore(:,i) = vecQ(:,i) .* tmp;
        end
    elseif strcmp(rank_type, 'eigXrel')
        %
        % according to the (eigenvector .* releance)
        %
        vecRankScore = zeros(N, K);
        for i = 1:K
            [x, y, z] = find(vecQ(:,i));
            if size(x,1) > 1 
                mean_vec = mean(eigIndex(x,:));
            else
                mean_vec = eigIndex(x,:);
            end
            tmp = eigIndex * mean_vec';
            vecRankScore(:,i) = vecQ(:,i) .* tmp;
        end
        vecRankScore = bsxfun(@times, vecRankScore, vecRel);
    elseif strcmp(rank_type, 'eigXrel2')
        %
        % according to the (eigenvector .* releance .^ 2)
        %
        vecRankScore = zeros(N, K);
        for i = 1:K
            [x, y, z] = find(vecQ(:,i));
            if size(x,1) > 1 
                mean_vec = mean(eigIndex(x,:));
            elseif size(x,1) == 1 
                mean_vec = eigIndex(x,:);
            else
                mean_vec = zeros(1, size(eigIndex, 2));
            end
            tmp = eigIndex * mean_vec';
            vecRankScore(:,i) = vecQ(:,i) .* tmp;
        end
        vecRankScore = bsxfun(@times, vecRankScore, vecRel);
    elseif strcmp(rank_type, 'eig+rel')
        %
        % according to the (eigenvector + relevance)
        %
        vecRankScore = zeros(N, K);
        for i = 1:K
            [x, y, z] = find(vecQ(:,i));
            if length(x) > 1 
                mean_vec = mean(eigIndex(x,:));
            else
                mean_vec = eigIndex(x,:);
            end
            tmp = eigIndex * mean_vec';
            vecRankScore(:,i) = vecQ(:,i) .* tmp;
        end
        vecRankScore = bsxfun(@plus, 0.1 * vecRankScore, 0.005 * vecRel);
    else
    end
end

