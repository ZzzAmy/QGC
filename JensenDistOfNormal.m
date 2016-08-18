function [ valDist ] = JensenDistOfNormal(vecC1_mean, matC1_var, vecC2_mean, matC2_var)
    d = size(vecC1_mean, 2);
    KL_2to1 = log(det(matC2_var) / det(matC1_var)) - d + trace((matC2_var^-1) * matC1_var) + (vecC2_mean - vecC1_mean) * (matC2_var^-1) * (vecC2_mean - vecC1_mean)';
    KL_1to2 = log(det(matC1_var) / det(matC2_var)) - d + trace((matC1_var^-1) * matC2_var) + (vecC1_mean - vecC2_mean) * (matC1_var^-1) * (vecC1_mean - vecC2_mean)';
    valDist = 0.5 * (KL_2to1 + KL_1to2);
end

