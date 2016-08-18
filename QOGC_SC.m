function [eigVec_ret, matrixV] = QOGC_SC(query, K, Q, vecRel)
    global matG
    global N
    
    if Q > N
        Q = N;
    end
    
    data = zeros(Q, N);
    vecIndex = zeros(Q, 1);
    
    [topQval, topQidx] = maxk(vecRel, Q);
    matGG = matG(topQidx, topQidx);
    
    vecD = sum(matGG);
    matD_sqrt = spdiags(1 ./ (vecD(:) .^ 0.5), 0, Q, Q);
    matG_AG = matD_sqrt * matGG * matD_sqrt;
    matG_AG = 0.5 * (matG_AG + matG_AG');
    
    
    kmeansRun = true;
    while (kmeansRun)
        opts.issym = 1;
        [eigVec, eigVal] = eigs(matG_AG, K, 'la', opts);

        c = kmeans(eigVec, K, 'emptyaction', 'drop');
        n = length(c);

        oPhi = sparse(1:n, c, ones(n,1), n, K);

        sizePhi = sum(oPhi);

        if any(sizePhi(:)==0) == 0
            kmeansRun = false;
        end
    end
    
    matrixV = zeros(N, K);
    eigVec_ret = zeros(N, K);
    for i = 1:Q
        if c(i) == 0
            c(i) = 1;
        end
        matrixV(topQidx(i, 1), c(i)) = 1;
        eigVec_ret(topQidx(i, 1), :) = eigVec(i,:);
    end
    
    eigVec_ret(query,:) = [];
end