function [weight_eigvec, vecQ, vecPerform1, vecBestPerform, NCut_BestRel, vecPerform2, vecBestPerform2, NCut_BestRel2, vecPerform3, vecBestPerform3, NCut_BestRel3] = QOGC_QGC_batch(query, K, H, M)
    global N
    
    global vecRel
    
    %rank_type = 'eig';
    rank_type = 'eigXrel';
    %rank_type = 'eigXrel2';
    
    %% QGC without clustering balance constraint
    ok = false;
    time = 0;
    while ok == false && time < 5
        time = time + 1;
        
        [weight_eigvec, vecQ] = QOGC_QGC(query, K, H, vecRel, 0, -0.02, 100);

        if nargin > 3 && length(M) > 1
            NCut_BestRel = zeros(length(M), 1);
            vecBestPerform = zeros(N, K, length(M));
            vecPerform1 = zeros(N, K, length(M));
            for i = 1:length(M)
                [A_NCut_BestRel, A_vecBestPerform, A_vecPerform1] = QOCut(query, weight_eigvec, vecQ, 'ncut', vecRel, rank_type, M(i), 'TotalN');
                NCut_BestRel(i) = A_NCut_BestRel;
                if size(A_vecBestPerform, 2) < K
                    A_vecBestPerform = [A_vecBestPerform zeros(N,1)];
                    A_vecPerform1 = [A_vecPerform1 zeros(N,1)];
                end
                vecBestPerform(:,:,i) = A_vecBestPerform;
                vecPerform1(:,:,i) = A_vecPerform1;
            end
        elseif nargin > 3 && M > 0
            [NCut_BestRel, vecBestPerform, vecPerform1] = QOCut(query, weight_eigvec, vecQ, 'ncut', vecRel, rank_type, M, 'TotalN');
        else
            [NCut_BestRel, vecBestPerform, vecPerform1] = QOCut(query, weight_eigvec, vecQ, 'ncut', vecRel, rank_type);
        end   
        vecNcut_rel = [NCut_BestRel];
        
        if nnz(sum(vecQ)) == K
            ok = true;
        end
    end
    
    %% QGC with cluster balance constraint by Lagrange multiplier
    ok = false;
    time = 0;
    while ok == false && time < 5
        time = time + 1;
        
        [weight_eigvec, vecQ] = QOGC_QGC(query, K, H, vecRel, 1, -0.02, 100);

        if nargin > 3 && length(M) > 1
            NCut_BestRel2 = zeros(length(M));
            vecBestPerform2 = zeros(N, K, length(M));
            vecPerform2 = zeros(N, K, length(M));
            for i = 1:length(M)
                [A_NCut_BestRel2, A_vecBestPerform2, A_vecPerform2] = QOCut(query, weight_eigvec, vecQ, 'ncut', vecRel, rank_type, M(i), 'TotalN');
                NCut_BestRel2(i) = A_NCut_BestRel2;
                if size(A_vecBestPerform2, 2) < K
                    A_vecBestPerform2 = [A_vecBestPerform2 zeros(N,1)];
                    A_vecPerform2 = [A_vecPerform2 zeros(N,1)];
                end
                vecBestPerform2(:,:,i) = A_vecBestPerform2;
                vecPerform2(:,:,i) = A_vecPerform2;
            end
        elseif nargin > 3 && M > 0
            [NCut_BestRel2, vecBestPerform2, vecPerform2] = QOCut(query, weight_eigvec, vecQ, 'ncut', vecRel, rank_type, M, 'TotalN');
        else
            [NCut_BestRel2, vecBestPerform2, vecPerform2] = QOCut(query, weight_eigvec, vecQ, 'ncut', vecRel, rank_type);
        end
        vecNcut_rel2 = [NCut_BestRel2];
        
        if nnz(sum(vecQ)) == K
            ok = true;
        end
    end
    
    %% QGC with cluster balance by setting large constant
%     [weight_eigvec3, vecQ] = QOGC_QGC(query, K, H, vecRel, 2, 100);
%     
%     if nargin > 3 && length(M) > 1
%         NCut_BestRel3 = zeros(length(M));
%         vecBestPerform3 = zeros(N, K, length(M));
%         vecPerform3 = zeros(N, K, length(M));
%         for i = 1:length(M)
%             [A_NCut_BestRel3, A_vecBestPerform3, A_vecPerform3] = QOCut(query, weight_eigvec3, vecQ, 'ncut', vecRel, 'eigXrel', M(i), 'TotalN');
%             NCut_BestRel3(i) = A_NCut_BestRel3;
%             vecBestPerform3(:,:,i) = A_vecBestPerform3;
%             vecPerform3(:,:,i) = A_vecPerform3;
%         end
%     elseif nargin > 3 && M > 0
%         [NCut_BestRel3, vecBestPerform3, vecPerform3] = QOCut(query, weight_eigvec3, vecQ, 'ncut', vecRel, 'eigXrel', M, 'TotalN');
%     else
%         [NCut_BestRel3, vecBestPerform3, vecPerform3] = QOCut(query, weight_eigvec3, vecQ, 'ncut', vecRel, 'eigXrel');
%     end
%     vecNcut_rel3 = [NCut_BestRel3];
    
    vecPerform3 = vecPerform2;
    vecBestPerform3 = vecBestPerform2;
    NCut_BestRel3 = NCut_BestRel2;
    
end