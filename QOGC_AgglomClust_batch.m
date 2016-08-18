function [ vecNcut_rel, vecBestPerform, NCut_BestRel ] = QOGC_AgglomClust_batch( query, K, Q, M )
    global N
    global vecRel

    [vecQ] = QOGC_AgglomClust(query, K, Q, vecRel);
    
    fprintf('Calculating QOCut...');
    
    if nargin > 3
        [NCut_BestRel, vecBestPerform, vecTemp] = QOCut(query, [], vecQ, 'ncut', vecRel, 'rel', M, 'TotalN');
    else
        [NCut_BestRel, vecBestPerform, vecTemp] = QOCut(query, [], vecQ, 'ncut', vecRel, 'rel');
    end
        
    vecNcut_rel = [NCut_BestRel];
end

