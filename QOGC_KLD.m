function [vecQ] = QOGC_KLD(query, K, vecRel)
    global matA

    global N

    vecV = rand(N, K);

    for i = 1:K
        vecV(:, i) = vecV(:, i) .* vecRel;
    end
    vecV(query, :) = 0;

    params = struct();
    [params,options] = process_options(params);
    options.Method = 'lbfgs';
    options.LS_type = 1;
    options.LS_interp = 2;
    options.MaxFunEvals = 10000;
    options.progTol = 10e-25;

    [x, f, exitflag, output] = minFunc(@ObjFunc6rel, vecV(:), options);
    params.tau = 0;
    vecV(:) = x;

    vecP_results = exp(matA * vecV);

    for i = 1:K
        vecP_results(:, i) = vecP_results(:, i) / sum(vecP_results(:, i));
    end

    vecQ = vecP_results;
    for i = 1: 30
        vecQ = 0.8 * matA * vecQ + 0.2 * vecP_results;
    end
end