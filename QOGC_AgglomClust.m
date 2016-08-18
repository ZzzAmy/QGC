%
% QOGC: Agglomerative Clustering
% 
function [matrixV] = QOGC_AgglomClust(query, K, Q, vecRel)
    global matG
    global N
    
    if Q > N
        Q = N;
    end
    
    data = zeros(Q, N);
    vecIndex = zeros(Q, 1);
    
    [topQVal, topQIdx] = maxk(vecRel, Q);
    data = matG(topQIdx, :);
    
%     topQ = maxk(vecRel, Q);
%     for i = 1:Q
%         indx = find(vecRel == topQ(i), 1);
%         data(i, :) = matA(indx, :);
%         vecIndex(i, 1) = indx;
%         vecRel(indx) = -100;
%     end
    
    %Z = linkage(full(matA), 'ward', '@Distance_nodes');
    %Z = linkage(full(matA), 'ward', 'euclidean');
    %matrixV = full(sparse(1:N, c, 1, N, K));
    
    Z = linkage(full(data), 'ward', 'cosine');
    %Z = linkage(data, 'ward', 'euclidean');
    c = cluster(Z, 'maxclust', K);
    
    matrixV = zeros(N, K);
    for i = 1:Q
        if c(i) == 0
            c(i) = 1;
        end
        matrixV(topQIdx(i, 1), c(i)) = 1;
    end
end