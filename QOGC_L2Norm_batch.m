function [vecNcut_rel, vecBestPerform, NCut_BestRel] = QOGC_L2Norm_batch(query, K, valDampFactor, alg_type)
    global N
    
    global vecRel
    
    [vecQ] = QOGC_L2Norm_2(query, K, valDampFactor, alg_type, vecRel);

    [NCut_BestRel, vecBestPerform, vecTemp] = QOCut(query, vecQ, 'ncut', vecRel);
    
    vecNcut_rel = [NCut_BestRel];
end