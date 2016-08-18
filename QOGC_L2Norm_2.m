%
% QOGC: L2Norm distiance maximization ver.2
%
% Propation Seeds Choosing Stratagy (alg_type): 
%
%   alg_type == 1:
%     1. always choose the vertex that has the most replacing reward.
%
%   alg_type == 2:
%     1. choose vertices in turn until candidate set converges. (like Gibbs Sampling)
%     2. always choose the vertex that has the most replacing reward.
%
%   alg_type == 3:
%     1. choose vertices greedily (also in turn) until candidate set converges.
%     2. always choose the vertex that has the most replacing reward.
%
%
function [matrixV] = QOGC_L2Norm(query, K, valDampFactor, alg_type, vecRel)
    global matA
    global matG 

    global N

    matrixR = spdiags(vecRel(:), 0, N, N);

    % Create identity matrix I_K
    matrixI_K = speye(K);

    % QOGC : EM algorithm
    isConverge = false;

    matrixL = zeros(2, K+1);
    matrixL2 = zeros(3, K+1);
    
    matrixV = zeros(N, K);

    if query > 0
        for i = 1:K
            matrixL(2, i+1) = query;
            matrixL2(2, i+1) = query;
            matrixV(query, i) = 1;
        end
    else
        for i = 1:K
            matrixL(2, i+1) = arrayIniLabel(i);
            matrixL2(2, i+1) = arrayIniLabel(i);
            matrixV(arrayIniLabel(i), i) = 1;
            matrixR = speye(N);
            vecRel = ones(N, 1);
        end
    end
    
    %print(matrixL(2,:), "\n");

    t = 0;

    matrixQ = sparse(matrixV);
    tmpM = valDampFactor * matrixQ;
    for i = 1:10
        matrixV = (1-valDampFactor) * (matA * matrixV) + tmpM;
    end
    clear tmpM;

    matrixM = K * matrixI_K - ones(K, K);
    matrixVdiff = matrixV * matrixM;

    % Initialize matU
    matU = zeros(N, K);
    vecU = 1/N * ones(N, 1);
    
    tmpM = valDampFactor * (vecRel .* (vecRel .* matrixV(:,1)));
    for i = 1:10
        vecU = (1-valDampFactor) * matA' * vecU + tmpM;
    end
    clear tmpM;
    
    for i = 1:K
      matU(:, i) = vecU;
    end
    clear vecU;

    %% Updating candidate label by EM algorithm
    updatelabIndex = 1;
    step = 0;
    
    % flag_strategy
    if alg_type == 1
        flag_strategy = 'GREEDY';
    elseif alg_type == 2
        flag_strategy = 'GIBBS_SAMPL';
    elseif alg_type == 3
        flag_strategy = 'GREEDY_IN_TURN';
    else
    end
    
    matrixL_batch = matrixL;
    choosen = zeros(K, 1);
    while isConverge == false
        t = t + 1;
        fprintf('\n Iter %d', t);
        vecL = matrixL(mod(t,2)+1, 2:K+1);
        
        %% Maximization Phase
        
        matVG = 1/N * ones(N, K);

        % calculate v_grad
        tmpM = (2 * valDampFactor) * (matrixR * (matrixR * matrixVdiff));
               
        for i = 1:10
            matVG = (1-valDampFactor) * matA' * matVG + tmpM;
        end
        
        clear tmpM;
        
        if query > 0
          matSol = 2*valDampFactor * matVG - 2*(K-1) * matU;
        else
          matSol = 2*valDampFactor * matVG;
        end

        arrayCandidateLab = zeros(K, 1);
        arrayCandidateVal = zeros(K, 1);
        
        for i = 1:N
            if vecRel(i ,1) == 0
               for j = 1:K
                   matSol(i, j) = -10e10;
               end
            end
        end

        for j = 1:K
            % vertex index of the j-th label
            l = matrixL(mod(t, 2) + 1, j+1);

            vec_l = (1-valDampFactor)^2 * matA * (matA * matVG(:, j));

            vecW = zeros(N, 1);

            % Approximate the vector vecW near the j-th label
            [x_idx, y_idx] = find(matG(:, l));

            for q = 1:length(x_idx)
                k = x_idx(q);
                v_k = (valDampFactor * (1-valDampFactor)) * matA(:, k) + vec_l;
                v_k(k) = v_k(k) + valDampFactor;
                vecW(k ,1) = norm(vecRel .* v_k)^2;
                clear v_k;
            end
            clear vec_l;

            if query > 0
                matSol(:, j) = matSol(:, j) + (K-1) * vecW;
            end

            clear vecW;
            
            % Avoid from assignning the labels that is the same with last run.
            for k = 1:K
                lab = matrixL(mod(t, 2) + 1 , k+1);
                matSol(lab, j) = -10e10;
            end

            [v, candidateLabel] = max(matSol(:, j));
            arrayCandidateLab(j, 1) = candidateLabel;

            vectt = matrixV(:, j);
            arrayCandidateVal(j, 1) = v - 2 * valDampFactor * matVG(l,j) + (K-1) * (vectt' * (vecRel .* (vecRel .* vectt)));
        end

        % choose a label to be updated
        if strcmp(flag_strategy, 'GIBBS_SAMPL')
            
            nextCatIndx = 0;
            for i=1:K
                if choosen(i, 1) == 0
                    nextCatIndx = i;
                    break;
                end
            end
            
            updatelabIndex = nextCatIndx;
            
            if nextCatIndx == K
                choosen = zeros(K, 1);
            else
                choosen(nextCatIndx, 1) = 1;
            end
            
        elseif strcmp(flag_strategy, 'GREEDY_IN_TURN')
            
            arrayCandidateVal = arrayCandidateVal - 10e5 * choosen;
            [bestVal, updatelabIndex] = max(arrayCandidateVal);
            choosen(updatelabIndex, 1) = 1;
            
            if sum(choosen) == K
                choosen = zeros(K, 1);
            end
            
        elseif strcmp(flag_strategy, 'GREEDY')
            [bestVal, updatelabIndex] = max(arrayCandidateVal);
        else
        end
        
        betterLabel = arrayCandidateLab(updatelabIndex);
        vecL(1, updatelabIndex) = betterLabel;

        %% Expectation Phase
        %vecV = zeros(N, 1);
        %vecV(betterLabel, 1) = 1;
        matrixV(:, updatelabIndex) = zeros(N, 1);
        matrixV(betterLabel, updatelabIndex) = 1;
        
        tmp_vecQ = valDampFactor * matrixV(:, updatelabIndex);
        for i = 1:10
          matrixV(:, updatelabIndex) = (1-valDampFactor) * (matA * matrixV(:, updatelabIndex)) + tmp_vecQ;
        end
        clear tmp_vecQ;
        %matrixV(:, updatelabIndex) = vecV;
        matrixVdiff = matrixV * matrixM;

        % update matU
        %vecU = 1/N * ones(N, 1);
        matU(:, updatelabIndex) = 1/N * ones(N, 1);
        
        tmpM = valDampFactor * (vecRel .* (vecRel .* matrixV(:, updatelabIndex)));
        for i = 1:10
            matU(:, updatelabIndex) = (1-valDampFactor) * (matA' * matU(:, updatelabIndex)) + tmpM;
        end
        clear tmpM;
        
        %matU(:, updatelabIndex) = vecU;

        %matRV = matrixR * matrixV;
        tmp = (matrixR * matrixV * matrixM);
        tmp2 = matrixV' * matrixR;
        valL = trace(tmp2 * tmp);
        matrixL(mod(t+1,2)+1, 1) = valL;
        matrixL2(mod(t+1,3)+1, 1) = valL;
        
        if mod(t, K) == 0
            matrixL_batch(mod(floor(t/k)+1,2)+1, 1) = matrixL(mod(t+1,2)+1, 1);
        end

        %% isCovergence determination
        if alg_type == 1
            % greedy
            if vecL == matrixL(mod(t+1,2)+1, 2:K+1)
                isConverge = true;
            end

            if t > 30
                isConverge = true;
            end
            
        elseif alg_type == 2
            % Gibbs Sampl => greedy
            if vecL == matrixL_batch(mod(floor(t/k)+1,2)+1, 2:K+1)
                if strcmp(flag_strategy, 'GIBBS_SAMPL')
                    flag_strategy = 'GREEDY';
                elseif strcmp(flag_strategy, 'GREEDY')
                    isConverge = true;
                else
                end
            end
            
            if t == 31
                flag_strategy = 'GREEDY';
            elseif t == 61
                isConverge = true;
            end
                      
        elseif alg_type == 3
            % Gibbs In Turn => greedy
            if vecL == matrixL_batch(mod(floor(t/k)+1,2)+1, 2:K+1)
                if strcmp(flag_strategy, 'GREEDY_IN_TURN')
                    flag_strategy = 'GREEDY';
                elseif strcmp(flag_strategy, 'GREEDY')
                    isConverge = true;
                else     
                end    
            end
            
            if t == 31
                flag_strategy = 'GREEDY';
            elseif t == 61
                isConverge = true;
            end
            
        else
        end
        
        matrixL(mod(t+1,2)+1, 2:K+1) = vecL;
        matrixL2(mod(t+1,3)+1, 2:K+1) = vecL;
        
        if mod(t, K) == 0
            matrixL_batch(mod(floor(t/k)+1,2)+1, 2:K+1) = vecL;
        end

        fprintf(' %s ', flag_strategy);
        fprintf(' %d ', matrixL(mod(t+1,2)+1, :));
        fprintf(' %d ', matrixL_batch(mod(floor(t/k)+1,2)+1, :));
    end

    %% Choose better label set
    matrixQ = zeros(N, K);
    if matrixL2(3,1) >= matrixL2(1, 1) && matrixL2(3,1) >= matrixL2(2, 1)
        for i = 1:K
            idx = matrixL2(3, i+1);
            matrixQ(idx, i) = 1;
        end
        fprintf('\n Result: ');
        fprintf(' %d ', matrixL2(3,:));
        fprintf('\n');
    elseif matrixL2(2,1) >= matrixL2(1, 1) && matrixL2(2,1) >= matrixL2(3, 1)
        for i = 1:K
            idx = matrixL2(2, i+1);
            matrixQ(idx, i) = 1;
        end
        fprintf('\n Result: ');
        fprintf(' %d ', matrixL2(2,:));
        fprintf('\n');
    else
        for i = 1:K
            idx = matrixL2(1, i+1);
            matrixQ(idx, i) = 1;
        end
        fprintf('\n Result: ');
        fprintf(' %d ', matrixL2(1,:));
        fprintf('\n');
    end
    matrixV = full(matrixQ);

    %
    % --- Expectation Phase ---
    %
    tmpM = valDampFactor * matrixQ;

    for i = 1:10
        matrixV = (1-valDampFactor) * matA' * matrixV + tmpM;
    end
end

