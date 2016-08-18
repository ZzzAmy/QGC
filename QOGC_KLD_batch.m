function [vecNcut_rel, vecBestPerform, NCut_BestRel] = QOGC_KLD_batch(query, K, times)
    global N
    
    global vecRel

    vecNcut_rel = zeros(times, 1);

    vecBestPerform = zeros(N, times);
    NCut_BestRel = inf;
    
    %% Run Query-oriented Graph Clustering
    for indx = 1:times
        
        [vecQ] = QOGC_KLD(query, K, vecRel);

        [NCut_rel, vecPerform, vecTemp] = QOCut(query, vecQ, 'ncut', vecRel);

        if isnan(NCut_rel) == 1
            indx = indx - 1;
        elseif NCut_rel < 0.1;
            vecNcut_rel(indx, 1) = NCut_rel;
        else
            vecNcut_rel(indx, 1) = NCut_rel;
        end
        
        if NCut_BestRel > NCut_rel
            NCut_BestRel = NCut_rel;
            vecBestPerform = vecPerform;
        end
    end
end