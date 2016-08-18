function [ vecR ] = PageRank( q, n, alpha )

    global matB
    %global matA

    if q == -1
        % PageRank
        vecP = ones(n ,1) / n;        
    else
        % Personalized PageRank
        vecP = zeros(n ,1); 
        vecP(q, 1) = 1;
    end
    
    vecR = vecP;
    
    for i = 1:50
        vecR = (1-alpha) * matB * vecR + alpha * vecP;
    end

end

