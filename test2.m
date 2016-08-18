%% Record test results
    if TEST_TYPE == 1
        matRes_Tree = zeros(topN, K);
        matRes_L2 = zeros(topN, K);
        matRes_KLD = zeros(topN, K);
        matRes_SC = zeros(topN, K);

        for k = 1:K
            topN_Res_Tree = maxk(Tree_vecBestPerform(:, k), topN);
            topN_Res_L2 = maxk(L2_vecBestPerform(:, k), topN);
            topN_Res_KLD = maxk(KLD_vecBestPerform(:, k), topN);
            topN_Res_SC = maxk(SC_vecBestPerform(:, k), topN);

            for n = 1:topN
                matRes_Tree(n, k) = find(Tree_vecBestPerform(:, k) == topN_Res_Tree(n), 1);
                matRes_L2(n, k) = find(L2_vecBestPerform(:, k) == topN_Res_L2(n), 1);
                matRes_KLD(n, k) = find(KLD_vecBestPerform(:, k) == topN_Res_KLD(n), 1);
                matRes_SC(n, k) = find(SC_vecBestPerform(:, k) == topN_Res_SC(n), 1);
            end
        end

        cellResult{i, 1} = query;
        cellResult{i, 2} = Tree_NCut_BestRel;
        cellResult{i, 3} = L2_NCut_BestRel;
        cellResult{i, 4} = KLD_NCut_BestRel;
        cellResult{i, 5} = SC_NCut_BestRel;
        cellResult{i, 6} = matRes_Tree;
        cellResult{i, 7} = matRes_L2;
        cellResult{i, 8} = matRes_KLD;
        cellResult{i, 9} = matRes_SC;
    else
        cellRes_Tree = cell(topN, K);
        cellRes_L2 = cell(topN, K);
        cellRes_KLD = cell(topN, K);
        cellRes_SC = cell(topN, K);
        
        for k = 1:K
            topN_Res_Tree = maxk(Tree_vecBestPerform(:, k), topN);
            topN_Res_L2 = maxk(L2_vecBestPerform(:, k), topN);
            topN_Res_KLD = maxk(KLD_vecBestPerform(:, k), topN);
            topN_Res_SC = maxk(SC_vecBestPerform(:, k), topN);
            
            fprintf('\nK = %d : ', k);

            for n = 1:topN
                id_Tree = find(Tree_vecBestPerform(:, k) == topN_Res_Tree(n), 1);
                id_L2N = find(L2_vecBestPerform(:, k) == topN_Res_L2(n), 1);
                id_KLD = find(KLD_vecBestPerform(:, k) == topN_Res_KLD(n), 1);
                id_SC = find(SC_vecBestPerform(:, k) == topN_Res_SC(n), 1);
                
                if Tree_vecBestPerform(id_Tree, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_Tree, 1);
                    if ~isempty(id)
                        cellRes_Tree(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_Tree{n, k} = 'n / a';
                    end
                    
                end
                if L2_vecBestPerform(id_L2N, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_L2N, 1);
                    if ~isempty(id)
                        cellRes_L2(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_L2{n, k} = 'n / a';
                    end
                end
                if KLD_vecBestPerform(id_KLD, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_KLD, 1);
                    if ~isempty(id)
                        cellRes_KLD(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_KLD{n, k} = 'n / a';
                    end
                end
                if SC_vecBestPerform(id_SC, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_SC, 1);
                    if ~isempty(id)
                        cellRes_SC(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_SC{n, k} = 'n / a';
                    end
                end
                
                fprintf('%d , ', id_SC);
                
                Tree_vecBestPerform(id_Tree, k) = 0;
                L2_vecBestPerform(id_L2N, k) = 0;
                KLD_vecBestPerform(id_KLD, k) = 0;
                SC_vecBestPerform(id_SC, k) = 0;
            end
        end
        cellResult{i, 1} = cellVertexNames{1,2}(find(cellVertexNames{1,1} == query));
        cellResult{i, 2} = Tree_NCut_BestRel;
        cellResult{i, 3} = L2_NCut_BestRel;
        cellResult{i, 4} = KLD_NCut_BestRel;
        cellResult{i, 5} = SC_NCut_BestRel;
        cellResult{i, 6} = cellRes_Tree;
        cellResult{i, 7} = cellRes_L2;
        cellResult{i, 8} = cellRes_KLD;
        cellResult{i, 9} = cellRes_SC;
        
        
        cellR = cell(10, 1);
        for n = 1:20
            [val , indx] = max(vecRel);
            cellR{n, 1} = cellVertexNames{1,2}(find(cellVertexNames{1,1} == indx, 1));
            vecRel(indx) = 0;
        end
    end
    
    fprintf('\n query: %d => ', query);
    
    fprintf('Tree: %f , ', Tree_NCut_BestRel);
    fprintf('L2Norm: %f , ', L2_NCut_BestRel);
    fprintf('KLD: %f , ', KLD_NCut_BestRel);
    fprintf('SC: %f \n', SC_NCut_BestRel);
    
    %fprintf('L2Norm_Greedy: %f , ', Tree_NCut_BestRel);
    %fprintf('L2Norm_GibSpl: %f , ', L2_NCut_BestRel);
    %fprintf('L2Norm_GdTurn: %f \n', KLD_NCut_BestRel);
    Tree_perf = full(Tree_vecBestPerform);
    L2_perf = full(L2_vecBestPerform);
    KLD_perf = full(KLD_vecBestPerform);
    SC_perf = full(SC_vecBestPerform);