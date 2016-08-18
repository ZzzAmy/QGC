function [ vecResult ] = QOGC_Tree(query, K, distType, valDampFactor)
%
% Function: QueryOrientedClustering (bound by top K clusters)
% Abstract:
%   In this function, we introduce two parameters: #numM and #numN.
%   #numM determines how many vertices are labelized at most. 
%   #numN determines how many distributions remains at most after merge funcitons.
%
    global matG
    global matA
    global N
    
    global matResultVal
    global matResultIndex
    
    iterMax = 5;
    
    %
    % while the top-k probabilty distributions would not be merged, the
    % last-2k distribuiotns are all been merged. Thus, we have two
    % explanaitons. 
    % First, the top-k distribuitons are different with one another enough.
    % Second, the last-2k distribuiotns are irrelevant to the query enough.
    %
    numM = 3 * K;
    numN = 2 * K;

    % Create the query of the multi-label.
    % As such, the query is represented by a set of vectors.
    matQ = sparse(N, 1);
    matQ(query, 1) = 1;
    matV = matQ;
    k = 0;

    %% Tree Search & Merge
    % matParentLink(i, j) = 1  ----->  v_j is a parent of v_i.
    matParentLink = sparse(N, N);
    %matLabelStep = sparse(N, 1);

    isConverge = false;
    isTopKMerged = true;
    
    numNewLabels = 0;

    while isConverge == false  
        k = k + 1;
        if k == iterMax || isTopKMerged == false;
            isConverge = true;
        end
        
        [vecChildX, vecChildY, vecChildVal] = find(ones(1, N) * matQ);
        queryNum = length(vecChildY);
        vecV_q = matV(:, 1);

        %% Label propagation
        for index = 1:queryNum
            % get the index of query vector
            q = vecChildY(index);

            % find all parent label-vectors index
            vecParent = matParentLink(q, 1:length(matQ(1,:)));
            vectemp = sign(matQ*vecParent');

            % calculate the matrix I_q
            matI_q = sparse(1:N, 1:N, 1);    

            [X, Y, Val] = find(vectemp);
            for i = 1:length(X)
                matI_q(X(i), X(i)) = 0;
            end

            % label propagation
            vecV = matV(:, q);
            vecQ = matQ(:, q);
            vecV = (1-valDampFactor) * matI_q * (matA * (matI_q * vecV)) + valDampFactor * vecQ;
            
%             if index == 1
%                 %vecV = (1-valPageRankAlpha) * matG * matrixD * matrixD * vecV + valPageRankAlpha * vecQ;
%                 vecV = (1-valDampFactor) * matI_q * (matA * (matI_q * vecV)) + valDampFactor * vecQ;
%             else
%                 vecV = (1-valDampFactor) * matI_q * (matA * (matI_q * vecV)) + valDampFactor * vecQ;
%             end

            matV(:, q) = vecV;
        end
%         matV = (1-valDampFactor) * matA * matV + valDampFactor * matQ;

        %% When k is even
        if mod(k, 2) == 0
            
            [V_X, V_Y, V_V] = find(ones(1, N) * matV);
            
            if length(V_Y) > K
            
                vecTopK = V_Y(2: K + 1);
                isTopKMerged = false;

                % --- Sibling merge ---

                [vecPaX, vecPaY, vecPaVal] = find(ones(1, N) * matParentLink);

                boolAllowMerge = true;

                while boolAllowMerge == true 

                    minChildID_i = -1;
                    minChildID_j = -1;
                    minDistance = 100;

                    % for each parent
                    IsMergeHappens = false;
                    parentID = -1;
                    for idx = 1:length(vecPaY)
                        if parentID == vecPaY(length(vecPaY) - idx + 1) && IsMergeHappens == false
                            continue;
                        end
                        parentID = vecPaY(length(vecPaY) - idx + 1);
                        IsMergeHappens = false;
                        [vecChildX, vecChildY, vecChildVal] = find(matParentLink(:, parentID));

                        vecHasMergedIndex = zeros(N, 1);

                        i = 0;
                        while i < length(vecChildX)
                            i = i + 1;
                            ChildID_i = vecChildX(i);

                            if vecHasMergedIndex(ChildID_i, 1) == 1
                                continue;
                            end

                            j = i;

                            while j < length(vecChildX)
                                j = j + 1;
                                ChildID_j = vecChildX(j);

                                if vecHasMergedIndex(ChildID_j, 1) == 1
                                    continue;
                                end

                                v_i = matV(:, ChildID_i);
                                v_j = matV(:, ChildID_j);

                                dist = Distance(v_i, v_j , distType);

                                if dist < minDistance
                                    minChildID_i = ChildID_i;
                                    minChildID_j = ChildID_j;
                                    minDistance = dist;
                                end
                            end
                        end     
                    end

                    [valChild_X, valChild_Y, valChild_Val] = find(matParentLink * ones(length(matParentLink(1,:)), 1));
                    valChildNum = length(valChild_X);

                    if valChildNum > numN                   
                        % i, j are similar enough

                        if ismember(minChildID_i, vecTopK) || ismember(minChildID_j, vecTopK)
                            isTopKMerged = true;
                        end

                        IsMergeHappens = true;
                        dimY = length(matV(1, :));

                        vecHasMergedIndex(minChildID_i, 1) = 1;
                        vecHasMergedIndex(minChildID_j) = 1;
                        vecChildX(length(vecChildX) + 1) = dimY + 1;

                        % Create a new vector v_k, and remove v_i and v_j.
                        dimY = length(matV(1, :));
                        matV(:, dimY + 1) = 0.5 * (matV(:, minChildID_i) + matV(:, minChildID_j));
                        matQ(:, dimY + 1) = 0.5 * (matQ(:, minChildID_i) + matQ(:, minChildID_j));
                        matV(:, minChildID_i) = zeros(N, 1);
                        matQ(:, minChildID_i) = zeros(N, 1);
                        matV(:, minChildID_j) = zeros(N, 1);
                        matQ(:, minChildID_j) = zeros(N, 1);

                        % Delete parent related to v_i and v_j by replcing to v_k
                        matParentLink(dimY + 1, :) = sign(matParentLink(minChildID_i, :) + matParentLink(minChildID_j, :));
                        matParentLink(:, dimY + 1) = sign(matParentLink(:, minChildID_i) + matParentLink(:, minChildID_j));
                        matParentLink(minChildID_i, :) = zeros(1, N);
                        matParentLink(minChildID_j, :) = zeros(1, N);
                        matParentLink(:, minChildID_i) = zeros(N, 1);
                        matParentLink(:, minChildID_j) = zeros(N, 1);

                        fprintf(strcat(int2str(minChildID_i), ' >SM< ', int2str(minChildID_j), ' => ', int2str(dimY + 1), '\n'));

                        % check the existence of siblings to determine to
                        % trigger Bottom-up merge.
                        [X, Y, V] = find(matParentLink(dimY + 1, :));

                        for idx = 1:length(Y)
                            parentID = vecPaY(idx);
                            [X, Y, Val] = find(matParentLink(:, parentID));

                            if length(X) ~= 1
                                continue;
                            end
                            ChildVecID = X(1);

                            % if a parent has only a child, remove its child.
                            matParentLink(ChildVecID, :) = sparse(1, N);
                            matParentLink(:, ChildVecID) = sparse(N, 1);
                            matV(:, ChildVecID) = sparse(N, 1);
                            matQ(:, ChildVecID) = sparse(N, 1);
                        end

                        boolAllowMerge = true;
                    else
                        boolAllowMerge = false;
                    end  
                end
            end
        end

        % the vertex that is visited for the first time
        [lenthX, lenthY] = size(matQ);
        vecAllQ = matQ * ones(lenthY, 1);
        
        [vecNewX, vecNewY, vecNewVal] = find(matV(:, 1) .* (sign(matV(:, 1)) - sign(vecAllQ)));
        %[vecNewX, vecNewY, vecNewVal] = find(matV(:, 1) .* (sign(matV(:, 1)) - sign(vecV_q)));
        
        [vecSortedNewVal, vecSortedNewIndex] = sort(vecNewVal,'descend');

        % the vertex that has been visited before
        [vecOldX, vecOldY, vecOldVal] = find(sign(vecV_q));

        %% New label generation
        for i = 1:length(vecSortedNewIndex)
            
            [x, y ,v] = find(ones(1, N) * matQ);
            
            numNewLabels = length(y) - 1;
            
            if numNewLabels >= numM
                break;
            end
            
%             numNewLabels = numNewLabels + 1;
            
            newLabelID = vecNewX(vecSortedNewIndex(i));

            % Initialize a new label
            vecNewLabel = sparse(N, 1);
            vecNewLabel(newLabelID, 1) = matV(newLabelID, 1);
            dimY = length(matV(1, :));
            matV(:, dimY + 1) = vecNewLabel;                                  
            matQ(:, dimY + 1) = vecNewLabel;

            % the neighbor vertices of newLabelID 
            [vecNearX, vecNearY, vecNearVal] = find(matG(:, newLabelID));

            % the parent vertices (visited neighbors) of newLabelID
            setParentNodes = intersect(vecNearX, vecOldX);

            vecParentNodes = sparse(1, N);
            for j = 1:length(setParentNodes)
                vecParentNodes(1, setParentNodes(j)) = 1;
            end

            % find all parent vector IDs from parent vertices
            [vecParentVecIDX, vecParentVecIDY, vecParentVecIDVal] = find(vecParentNodes * matQ);

            for j = 1:length(vecParentVecIDY)
                matParentLink(dimY + 1, vecParentVecIDY(j)) = 1;
            end
        end
        
        fprintf('The Length of matQ in itr');
        fprintf(strcat(int2str(k), ' : ', int2str(length(matQ(1,:))), '\n'));
    end
    
    [x, y, v] = find(sum(matV));
    for i = 1:K
        vecResult(:, i) = matV(:, y(i+1));
    end
    
    %vecResult = matV(:, 2:1+K);
    %{
    matSortedVal = sparse(length(matV(:, 1)), length(matV(1, :)));
    matSortedIndex = sparse(length(matV(:, 1)), length(matV(1, :)));

    for i = 2:length(matV(1, :))
        if sum(matV(:, i)) == 0
            continue;
        end
        [matSortedVal(:, i), matSortedIndex(:, i) ] = sort(matV(:, i), 'descend');
    end
    [SV_X, SV_Y, SV_V] = find(matSortedVal);
    
    setYindex = intersect(SV_Y, SV_Y);
    
    matResultVal = zeros(valK, length(setYindex));
    matResultIndex = zeros(valK, length(setYindex));
    
    for i = 1:length(setYindex)
        index = setYindex(i);
        matResultVal(:, i) = full(matSortedVal(1:valK, index));
        matResultIndex(:, i) = full(matSortedIndex(1:valK, index));
    end
    %}
end

