function [vecNcut_rel, vecBestPerform, NCut_BestRel] = QOGC_Tree_batch(query, K, distType, valDampFactor)
    global N
    global vecRel

    fprintf('QOGC Tree Merging...');
    [vecQ] = QOGC_Tree(query, K, distType, valDampFactor);
    
    fprintf('Calculating QOCut...');
    [NCut_BestRel, vecBestPerform, vecTemp] = QOCut(query, vecQ, 'ncut', vecRel);
        
    vecNcut_rel = [NCut_BestRel];
end