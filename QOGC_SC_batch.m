function [ vecNcut_rel, vecBestPerform, NCut_BestRel ] = QOGC_SC_batch( query, K, Q, M )
    global N
    global vecRel

    [eigQ, vecQ] = QOGC_SC(query, K, M, vecRel);
    
    fprintf('Calculating Spectral Clustering ...');
    
    [NCut_BestRel, vecBestPerform, vecTemp] = QOCut(query, [], vecQ, 'ncut', vecRel, 'rel', M, 'TotalN');
        
    vecNcut_rel = [NCut_BestRel];
end

