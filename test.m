global matG
global matA
global matB

global cellVertexNames

global vecRel

global maxitr

global N
global E
global K
global alpha

alpha = 0.1;
t = 10;

maxitr = 200;

% 
% Test Type => 1: toy graph; 2: UIT graph; 3: UIT graph; 4: DBLP graph
%
TEST_TYPE = 1;

%
% Re-read => 0: not re-read; 1: re-read
%
REREAD = 1;


K = 3;
topN = 0;
topHop = 1;


% algorithm active
act_QOGC_Tree = 0;
act_QOGC_L2Norm = 0;
act_QOGC_KLD = 0;
act_QOGC_QGC = 1;
act_QOGC_QGC_topAll = 0;
act_QOGC_AC = 0;
act_QOGC_SC = 0;


%% Load data
if REREAD == 1
    if TEST_TYPE == 1
        %[N, E] = LoadGraphData('simpleGraph.txt');
        %[N, E] = LoadGraphData('tinyGraph.txt');
        [N, E] = LoadGraphData('simpleGraph_new.txt');
    elseif TEST_TYPE == 2
        [N, E] = LoadGraph_UIT('tags.dat', 'user_taggedartists.dat', 't', 5, 10);
    elseif TEST_TYPE == 3
        %[N, E] = LoadGraph_UIT('/Users/iankuoli/Dataset/LastFm2/artists2.txt', 'user_taggedartists.dat', 'i', 1, 5);
        [N, E] = LoadGraph_UIT('/Users/iankuoli/Dataset/LastFm2/artists2.txt', 'user_taggedartists.dat', 'i', -1, 10);
        %[N, E] = LoadGraph_UIT('C:/Dataset/LastFm2/artists2.txt', 'user_taggedartists.dat', 'i', 1, 5);
    elseif TEST_TYPE == 4
        [N, E] = LoadGraph_DBLP('/Users/iankuoli/Dataset/BigDND_DBLP/graph2012_2.txt', '/Users/iankuoli/Dataset/BigDND_DBLP/author2012_2.txt');
        %[N, E] = LoadGraph_DBLP('/C:/Users/Ian/Dropbox/dataset/BigDND_DBLP/graph2012_2.txt', 'C:/Users/Ian/Dropbox/dataset/BigDND_DBLP/author2012_2.txt');
    elseif TEST_TYPE == 5
        [N, E] = LoadGraph_IMDB('/Users/iankuoli/Dataset/IMDb/ComedyActors_(2000-2013).txt', '/Users/iankuoli/Dataset/IMDb/ComedyCountries_(2000-2013).txt', 0, 10);
        %[N, E] = LoadGraph_IMDB('C:/Users/Ian/Dropbox/dataset/IMDb/ComedyActors_(2000-2013).txt', 'C:/Users/Ian/Dropbox/dataset/IMDb/ComedyCountries_(2000-2013).txt', 0, 10);
    elseif TEST_TYPE == 6
        [N, E] = LoadGraph_Facebook('/Users/iankuoli/Dataset/sna_facebook/facebook_combined.txt');
    elseif TEST_TYPE == 7
        [N, E] = LoadGraph_ANN('/Users/iankuoli/Dataset/ANN/author-citation-network2012.txt');
    else
    end
end

%% Create symmetric normalized matrix A of G
vecTmp = sum(matG)';
[vecPosX, vecPosY, vecVal] = find(vecTmp);
matD = sparse(vecPosX, vecPosX, 1 ./ (vecVal .^ 0.5), N, N);
matA = matD * matG * matD;
matB = matG * matD * matD;

% Create identity matrix I
vecposition = zeros(N, 1);
for i = 1:N
    vecposition(i) = i;
end
matI = sparse(vecposition, vecposition, ones(1, N), N, N);

% Construct matA
xi = 0;%0.2;
matA = (1-xi) * matA + xi * matI;

clear vecTmp;
clear vecPosX;
clear vecPosY;
clear vecVal;
clear vecposition;
clear matI;
clear matD;

%% Test
if TEST_TYPE == 1
    %for Simple graph
    %queries = [1 2 3 4 5 6 7 8 9 10];
    queries = [1];
elseif TEST_TYPE == 2
    % for Last.fm tag recmmendation
    queries = [1 3 4 5 6 19 21 23 32 39 41 56 62 67 75 78 79 80 81 82 83 86 87 90 91 102 109 121 127 130 131 139 140 145 164 179 181 185 189 190 191 192 193 195 203 215 216 231 233 234 238 245 247 250 252 292 296 297 308 311 315 347 352 378 385 386 387 388 389 390 391 392 393 397 444 445 446 481 499 503 505 511 517 527 528 541 545 547 567 570 574 575 593 604 612 615 617 638 647 686 689 696 700 701 721 735 736 739 745 757 761 781 789 792 793 797 821 827 829 848 853 855 869 873 887 895 918 927 935 936 996 1005 1016 1017 1076 1095 1107 1150 1155 1160 1169 1173 1196 1206 1219 1220 1228 -1 1331 1340 1373 1424 1434 1495 1547 1683 2011];
    %queries = [5];
elseif TEST_TYPE == 3
    % for Last.fm artist recommendation
    %queries = [-1];
    queries = [25 47 89 154 159 163 203 212 215 227 229 234 238 251 278 288 289 295 298 306 ...
               329 333 344 374 390 394 399 409 412 418 420 424 425 428 429 441 444 455 458 511 ...
               517 522 533 554 562 599 609 610 611 613 630 683 703 707 709 724 728 734 735 868 ... 
               875 917 928 1117 1120 1121 1123 1130 1138 1242 1249 1262 1276 1358 1369 1372 1412 1418 1519 1814 ...
               2265 2292 2606 3357 3371 5076 7123 8843 9740 11480, 6 15 20 25 39 46 63 65 81 103 ...
               166 170 171 172 173 190 199 207 208 209 210 220 225 230 233 250 253 257 360 377 ...
               378 382 385 393 396 416 433 489 503 529 534 576 578 585 605 612 614 618 620 621 ...
               631 689 710 715 716 718 730 732 733 735 812 820 841 843 846 903 918 919 920 923 ...
               924 937 949 951 952 956 969 972 985 986 987 1001 1042 1048 1073 1083 1099 1101 1103 1104 ...
               1109 1122 1132 1150 1175 1244 1254 1259 1260 1264 1275 1298 1342 1384 1390 1394 1398 1411 1413 1469];
    %queries = [533];   % Oasis
    %queries = [47];    % Emperor
    %queries = [599];    % David Bowie
    %queries = [25];     % Cradle of Filth
    %queries = [1412];   % Led Zeppelin
    %queries = [611];
elseif TEST_TYPE == 4
    % ------ for DBLP co-authorship #2000 -----
    %queries = [128 465 567 582 613 794 825 836 945 1025 1078 2049 2148 2167 2234 2730 2851 2982 2993 3220 3322 3494 3649 3854 4095 4385 4789 5475 5599 5656 5831 5824 8585 5906 6009 6182 6323 6706 6803 6958 7166 7504 7694 7876 8070 8284 8431 8667 8865 8981 9013 9056 9156 9573 9638 9846 9895 10149 10354 10375 10444 10912 11259 11946 13009 13671 14270 14277 14968 15445 15969 16451 17022 17268 17786 18246 18715 19394 20466 21315 22078 22116 22815 23463 24759 25374 26414 27194 28611 28595 29599 30411 31807 32961 33362 33579 34135 34581 37636 39475 57647 44793 48473 106857 146114 147383 150160 47365 4785 63725 48424 90704 91463 92192 92760 93623 94583 95549 96285 97371];
    %queries = [4785];   % Ming-Syan Chen
    %queries = [63725];  % Geoffrey E. Hinton
    %queries = [48424];  % Yoshua Bengio
    %queries = [22116];  % Yann LeCun
    
    
    % ------ for DBLP co-authorship #2010 ------ 
%     queries = [35750 7176 34195 65780 10849 10847 30369 1 2 3 6 9 11 470 634 699 790 876 993 1000 ...
%                1051 1226 1243 1312 1318 1502 1528 1621 1627 1678 2055 -1 2089 2524 2705 3076 3107 3300 3614 3754 ...
%                -1 3882 4015 4026 4112 4200 4201 4351 4388 4406 4900 4961 5054 5115 5140 5160 -1 5378 5837 6033 ...
%                6569 6712 6743 7001 7059 7061 7067 7071 7505 7513 7562 7565 7566 7608 8560 9030 9329 9435 9837 -1 ...
%                -1 -1 10021 10158 10239 10265 10401 10488 10586 10974 11110 11113 11115 11116 11933 12000 12030 12071 12171 -1 ...
%                12541 13172 13182 13205 13553 13560 14691 14793 14830 14863 16090 16111 16135 16209 16373 16400 16416 16819 17315 17325];
    
    %queries = [35750];  % Ming-Syan Chen
    %queries = [7176];   % Philip S. Yu
    %queries = [34195];  % Jian Pei
    %queries = [65780];  % Jure Leskovec
    %queries = [10849];  % Geoffrey E. Hinton
    %queries = [10847];  % Yoshua Bengio
    %queries = [30369];  % Yann LeCun
    
    
    % ------ for DBLP co-authorship #2012 ------ 
    % 14445 36 668 1941 8126 1601 4689 4703 4926 302 8593
    %queries = [5144 4662 31579 6699 23680 20167  8  9220 8268 8345 8362 11585 8247 8248 333 397 8565 461 ...
    %           353 1116 1133 8571 9720 7562 7604 7609 7614 7616 7649 7646 
    queries = [7790 7861 11268 9865 8246  ...
               675 11933 4430 714 732 797 827 12263 852 855 955 957 1020 1026 11452 5300 1461 1508 1524 1530 ...
               5349 5395 1574 1593 1602 5429 1627 1657 5512 1804  1976 5572 11509 10302 2310 2328 2481 2572 ...
               2637 2656 2657 10319 2719 2855 3001 5744 3352 3381 1145 3688 8077 3753 3832 3835 3903  4330 4387 ...
               10347 4522 4523 4555 4577 4623 4662  8161 4776 4780 4837  4927 5033 8208 8242 5148 5264 ...
               1601 4689 4703 4926 302 8593 14445
               ];
           
    %queries = [5144];   % Ming-Syan Chen
    %queries = [4662];   % Philip S. Yu
    %queries = [31579];  % Jian Pei
    %queries = [6699];   % Jure Leskovec
    %queries = [23680];  % Geoffrey E. Hinton
    %queries = [20167];  % Yoshua Bengio
    %queries = [14445];  % Yann LeCun
    
elseif TEST_TYPE == 5
    %queries = [1 3 4];
    queries = [1 3 4 5 6 19 21 23 32 39 41 56 62 67 75 78 79 80 81 82 83 86 87 90 91 102 109 121 127 130 131 139 140 145 164 179 181 185 189 ...
               191 192 193 203 215 216 231 233 234 238 245 247 250 252 292 296 297 308 311 315 347 378 385 386 387 388 390 391 392 393 397 444 ...
               445 446 481 499 503 505 511 517 527 528 541 545 547 567 570 574 575 593 604 612 615 617 638 647 686 689 696 700 701 721 735 736 739 ...
               745 757 761 781 789 792 793 797 821 827 829 848 853 855 869 873 887 895 918 927 935 936 996 1005 1016 1017 1076 1095 1107 1150 1155 1160 ...
               1169 1173 1196 1206 1219 1220 1228 1330 1331 1340 1373 1424 1434 1495 1547 1683 2011 2013 2145 4564 3432 4333 5445 5678 6787 6789 7124 7178 ...
               7256 7302 7551 7552 7553 7554 7565 7694 7698 7700 7825 7846 7936 8034 8078 8163 8275 8496 8599 8743 8888 8998 9019 9111 9211 9324 9345 9574 ...
               9786 9985 10004 10038 10098 11034 11486 11756];
elseif TEST_TYPE == 6
    queries = [0 107 348 414 686 698 1684 1912 3437 3980];
    queries = [-2 -2 348 414 686 -2 946 993 1059 1078 -2 1126 1184 1185 -2 1211 1238 -2 1367 1376 1377 1390 1391 -2 1431 1471 1516 1522 1551 1554 1557 1559 1610 1612 1613 1621 1622 1663 -2 1707 1714 -2 -2 1736 1746 -2 1799 -2 -2 -2 -2 1833 1835 1839 1888 1912 1917 1918 1929 1938 -2 1943 1946 1962 1966 1971 1979 1983 1984 1985 1986 1993 2005 -2 2020 2030 2033 2037 2040 2043 2045 -2 2054 2059 2064 2069 2073 2074 2078 2081 2087 2088 2090 2093 2095 2103 2104 2108 2111 2112 2118 2121 2123 2124 2131 -2 2139 2140 -2 2150 2154 2172 2184 2188 2190 2199 2200 2201 2206 2218 2220 2229 2233 2240 2244 2253 2257 2266 2271 2275 2276 2278 2282 2283 -2 2290 2299 2308 2309 2323 2324 2326 -2 2331 2333 2339 2340 2347 2348 2352 2354 2356 2363 2369 2370 2374 2376 2381 2384 2395 2404 2408 2409 2410 2414 2423 2428 2446 2460 2464 -2 2482 2492 2500 2504 2507 2520 2526 2542 2543 2549 2550 2551 2553 2559 2560 2561 2564 2573 2578 2586 2590 2593 2598 2600 2601 2602 2604 2607 2611 2615 2619 2624 2625 2630 2638 2654 2655 -2];
    queries = queries + 1;
elseif TEST_TYPE == 7
    queries = [25 47 89 154 159 163 203 212 215 227 229 234 238 251 278 288 289 295 298 306 ...
               329 333 344 374 390 394 399 408 412 418 420 424 425 428 429 441 444 455 458 511 ...
               517 522 533 554 562 599 609 610 611 613 630 683 705 706 710 724 728 738 740 868 ... 
               875 917 928 1117 1120 1119 1125 1130 1138 1242 1249 1262 1276 1367 1369 1373 1412 1418 1519 1815 ...
               2265 2290 2606 3358 3371 5077 7121 8843 9740 11498, 6 15 20 25 39 46 63 65 81 103 ...
               166 170 171 172 173 190 199 207 208 209 210 220 225 230 233 250 253 257 361 377 ...
               378 382 385 393 396 416 433 489 503 529 534 576 578 585 605 612 614 618 620 621 ...
               631 689 710 715 716 718 730 732 743 746 812 820 841 843 846 903 918 919 920 923 ...
               926 938 949 951 952 956 969 972 985 986 987 1001 1102 1104 1109 1122 1136 1151 1175 1244 ...
               1254 1259 1260 1264 1275 1298 1342 1384 1390 1394 1398 1411 1413 1469 2195 2213 2245 2253 2288 2316];
else
end

%
% query | Tree_NCut | L2Norm_NCut | KLD_NCut | Tree_topN | L2Norm_topN | KLD_topN
%
cellResult = cell(length(queries), 19);

for i = 1:length(queries)
    %% Set the queries and their relevances
    query = queries(i);
    
    if query < 0
        continue;
    end
    
    % if the query is isolated, ignore it.
    if query > 0 && sum(matA(query, :)) == 0.2
        continue;
    end
    
    %sum(matA(445,:))

    % Get the relevance vector
    if query > 0
        vecRel = PageRank(query, N, 0.8);
        vecRel(query, 1) = 0;
    else
        vecRel = ones(N, 1)/N;
    end

    if sum(vecRel) == 0
        continue;
    end
    
    
    if topN == 0 && topHop > 0
        % assign the ?-hop expanded neighbors
        vecExpand = zeros(N, 1);
        vecExpand(query, 1) = 1;
        
        for q = 1:topHop
            vecExpand = (matG * vecExpand + vecExpand) > 0;
        end
        
        topN = sum(vecExpand) - 1;
    end
    
    %% Algorithm: QOGC_Tree (Tree Search & Merge)
    if act_QOGC_Tree == 1
        fprintf('\n Run Tree_Search algorithm...\n');
        [Tree_vecNcut_rel, Tree_vecBestPerform, Tree_NCut_BestRel] = QOGC_Tree_batch(query, K, 'jsd', 0.7);
        %[Tree_vecNcut_rel, Tree_vecBestPerform, Tree_NCut_BestRel] = QOGC_L2Norm_batch(query, K, 0.7, 1);
    else
        Tree_NCut_BestRel = 0;
        Tree_vecBestPerform = sparse(N,K);
    end
    
    %% Algorithm: QOGC2_L2Norm (L2-norm Maximization by Label Propagation)
    if act_QOGC_L2Norm == 1
        fprintf('\n Run L2 distance maximization algorithm...\n');
        [L2_vecNcut_rel, L2_vecBestPerform, L2_NCut_BestRel] = QOGC_L2Norm_batch(query, K, 0.7, 1);
    else
        L2_NCut_BestRel = 0;
        L2_vecBestPerform = sparse(N,K);
    end
    
    %% Algorithm: QOGC5_KLD (KL-divergence maximization)
    if act_QOGC_KLD == 1
        fprintf('\n Run KL divergence maximization algorithm...\n');
        [KLD_vecNcut_rel, KLD_vecBestPerform, KLD_NCut_BestRel] = QOGC_KLD_batch(query, K, K*3);
        %[KLD_vecNcut_rel, KLD_vecBestPerform, KLD_NCut_BestRel] = QOGC_L2Norm_batch(query, K, 0.7, 3);
    else
        KLD_NCut_BestRel = 0;
        KLD_vecBestPerform = sparse(N,K);
    end
    
    %% Algorithm: QOGC5_QGC (Query-oriented spectral clustering)
    if act_QOGC_QGC == 1
        fprintf('\n Run Query-oriented spectral clustering algorithm...\n');
        [iPhi, vecQ, vecPerform1, QGC_vecBestPerform, QGC_NCut_BestRel, vecPerform2, QGC_vecBestPerform2, QGC_NCut_BestRel2, vecPerform3, QGC_vecBestPerform3, QGC_NCut_BestRel3] = QOGC_QGC_batch(query, K, 3, topN);
        c_size1 = sum(full(vecPerform1))
        c_size2 = sum(full(vecPerform2))
        c_size3 = sum(full(vecPerform3))
        entropy_QGC1 = funcEntropy(c_size1');
        entropy_QGC2 = funcEntropy(c_size2');
        entropy_QGC3 = funcEntropy(c_size3');
    else
        if act_QOGC_QGC_topAll == 0
            QGC_NCut_BestRel = 0;
            QGC_vecBestPerform = sparse(N,K);
            QGC_NCut_BestRel2 = 0;
            QGC_vecBestPerform2 = sparse(N,K);
            QGC_NCut_BestRel3 = 0;
            QGC_vecBestPerform3 = sparse(N,K);
            entropy_QGC1 = 0;
            entropy_QGC2 = 0;
            entropy_QGC3 = 0;
        end
    end
    
    %% Algorithm: QOGC5_QGC (Query-oriented spectral clustering) for all M which can be represented as a vector.
    if act_QOGC_QGC_topAll == 1
        fprintf('\n Run Query-oriented spectral clustering algorithm...\n');
        topM = [0, 0, 0, 0];
        
        vecExpand = zeros(N, 1);
        vecExpand(query, 1) = 1;
        for q = 1:4
            vecExpand = (matG * vecExpand + vecExpand) > 0;
            topM(q) = sum(vecExpand) - 1;
        end
        
        [iPhi, vecQ, vecPerform1, QGC_vecBestPerform, QGC_NCut_BestRel, vecPerform2, QGC_vecBestPerform2, QGC_NCut_BestRel2, vecPerform3, QGC_vecBestPerform3, QGC_NCut_BestRel3] = QOGC_QGC_batch(query, K, 3, topM);
        entropy_QGC1 = zeros(length(topM), 1);
        entropy_QGC2 = zeros(length(topM), 1);
        entropy_QGC3 = zeros(length(topM), 1);
        for m = 1:length(topM)
            c_size1 = sum(full(vecPerform1(:,:,m)));
            c_size2 = sum(full(vecPerform2(:,:,m)));
            c_size3 = sum(full(vecPerform3(:,:,m)));
            entropy_QGC1(m) = funcEntropy(c_size1');
            entropy_QGC2(m) = funcEntropy(c_size2');
            entropy_QGC3(m) = funcEntropy(c_size3');
        end
    else
        if act_QOGC_QGC == 0
            QGC_NCut_BestRel = 0;
            QGC_vecBestPerform = sparse(N,K);
            QGC_NCut_BestRel2 = 0;
            QGC_vecBestPerform2 = sparse(N,K);
            QGC_NCut_BestRel3 = 0;
            QGC_vecBestPerform3 = sparse(N,K);
            entropy_QGC1 = 0;
            entropy_QGC2 = 0;
            entropy_QGC3 = 0;
        end
    end
    
    
    %% Algorithm: QOGC5_AC (Query-oriented agglomerative clustering)
    if act_QOGC_AC == 1 && topN > 0
        fprintf('\n Run Query-oriented agglomerative clustering algorithm...\n');
        [AC_vecNcut_rel, AC_vecBestPerform, AC_NCut_BestRel] = QOGC_AgglomClust_batch(query, K, topN * K * 1, topN);
        c_size = sum(full(AC_vecBestPerform))
        entropy_AC = funcEntropy(c_size');
    else
        AC_NCut_BestRel = 0;
        AC_vecBestPerform = sparse(N,K);
        entropy_AC = 0;
    end
    
    %% Algorithm: QOGC5_SC (Query-oriented spectral clustering)
    if act_QOGC_SC == 1 && topN > 0
        fprintf('\n Run Query-oriented spectral clustering algorithm...\n');
        [SC_vecNcut_rel, SC_vecBestPerform, SC_NCut_BestRel] = QOGC_SC_batch(query, K, topN * K * 1, topN);
        c_size = sum(full(SC_vecBestPerform))
        entropy_SC = funcEntropy(c_size');
    else
        SC_NCut_BestRel = 0;
        SC_vecBestPerform = sparse(N,K);
        entropy_SC = 0;
    end 
    
    %% Record test results
    if TEST_TYPE == 1
        matRes_Tree = zeros(topN, K);
        matRes_L2 = zeros(topN, K);
        matRes_KLD = zeros(topN, K);
        matRes_QGC = zeros(topN, K);
        matRes_QGC2 = zeros(topN, K);
        matRes_QGC3 = zeros(topN, K);
        matRes_AC = zeros(topN, K);
        matRes_SC = zeros(topN, K);
        
        Tree_tmpPerf = Tree_vecBestPerform;
        L2_tmpPerf = L2_vecBestPerform;
        KLD_tmpPerf = KLD_vecBestPerform;
        QGC_tmpPerf = QGC_vecBestPerform;
        QGC_tmpPerf2 = QGC_vecBestPerform2;
        QGC_tmpPerf3 = QGC_vecBestPerform3;
        AC_tmpPerf = AC_vecBestPerform;
        SC_tmpPerf = SC_vecBestPerform;

        for k = 1:K
            topN_Res_Tree = maxk(Tree_tmpPerf(:, k), topN);
            topN_Res_L2 = maxk(L2_tmpPerf(:, k), topN);
            topN_Res_KLD = maxk(KLD_tmpPerf(:, k), topN);
            topN_Res_QGC = maxk(QGC_tmpPerf(:, k), topN);
            topN_Res_QGC2 = maxk(QGC_tmpPerf2(:, k), topN);
            topN_Res_QGC3 = maxk(QGC_tmpPerf3(:, k), topN);
            topN_Res_AC = maxk(AC_tmpPerf(:, k), topN);
            topN_Res_SC = maxk(SC_tmpPerf(:, k), topN);

            for n = 1:topN
                if topN_Res_Tree(n) > 0
                    id_Tree = find(Tree_tmpPerf(:, k) == topN_Res_Tree(n), 1);
                    matRes_Tree(n, k) = id_Tree;
                    Tree_tmpPerf(id_Tree, k) = 0;
                end
                if topN_Res_L2(n) > 0
                    id_L2N = find(L2_tmpPerf(:, k) == topN_Res_L2(n), 1);
                    matRes_L2(n, k) = id_L2N;
                    L2_tmpPerf(id_L2N, k) = 0;
                end
                if topN_Res_KLD(n) > 0
                    id_KLD = find(KLD_tmpPerf(:, k) == topN_Res_KLD(n), 1);
                    matRes_KLD(n, k) = id_KLD;
                    KLD_tmpPerf(id_KLD, k) = 0;
                end
                if topN_Res_QGC(n) > 0
                    id_QGC = find(QGC_tmpPerf(:, k) == topN_Res_QGC(n), 1);
                    matRes_QGC(n, k) = id_QGC;
                    QGC_tmpPerf(id_QGC, k) = 0;
                end
                if topN_Res_QGC2(n) > 0
                    id_QGC2 = find(QGC_tmpPerf2(:, k) == topN_Res_QGC2(n), 1);
                    matRes_QGC2(n, k) = id_QGC2;
                    QGC_tmpPerf2(id_QGC2, k) = 0;
                end
                if topN_Res_QGC3(n) > 0
                    id_QGC3 = find(QGC_tmpPerf3(:, k) == topN_Res_QGC3(n), 1);
                    matRes_QGC3(n, k) = id_QGC3;
                    QGC_tmpPerf3(id_QGC3, k) = 0;
                end
                if topN_Res_AC(n) > 0
                    id_AC = find(AC_tmpPerf(:, k) == topN_Res_AC(n), 1);
                    matRes_AC(n, k) = id_AC;
                    AC_tmpPerf(id_AC, k) = 0;
                end  
                if topN_Res_SC(n) > 0
                    id_SC = find(SC_tmpPerf(:, k) == topN_Res_SC(n), 1);
                    matRes_SC(n, k) = id_SC;
                    SC_tmpPerf(id_SC, k) = 0;
                end
            end
        end

        cellResult{i, 1} = query;
        cellResult{i, 2} = QGC_NCut_BestRel;
        cellResult{i, 3} = entropy_QGC1;
        cellResult{i, 4} = QGC_NCut_BestRel2;
        cellResult{i, 5} = entropy_QGC2;
        cellResult{i, 6} = QGC_NCut_BestRel3;
        cellResult{i, 7} = entropy_QGC3;
        cellResult{i, 8} = AC_NCut_BestRel;
        cellResult{i, 9} = entropy_AC;
        cellResult{i, 10} = SC_NCut_BestRel;
        cellResult{i, 11} = entropy_SC;
        cellResult{i, 12} = matRes_QGC;
        cellResult{i, 13} = matRes_QGC2;
        cellResult{i, 14} = matRes_QGC3;
        cellResult{i, 15} = matRes_AC;
        cellResult{i, 16} = matRes_SC;
        
    elseif TEST_TYPE > 1 && act_QOGC_QGC_topAll == 0
        
        cellRes_Tree = cell(topN, K);
        cellRes_L2 = cell(topN, K);
        cellRes_KLD = cell(topN, K);
        cellRes_QGC = cell(topN, K);
        cellRes_QGC2 = cell(topN, K);
        cellRes_QGC3 = cell(topN, K);
        cellRes_AC = cell(topN, K);
        cellRes_SC = cell(topN, K);
        
        Tree_tmpPerf = Tree_vecBestPerform;
        L2_tmpPerf = L2_vecBestPerform;
        KLD_tmpPerf = KLD_vecBestPerform;
        QGC_tmpPerf = QGC_vecBestPerform;
        QGC_tmpPerf2 = QGC_vecBestPerform2;
        QGC_tmpPerf3 = QGC_vecBestPerform3;
        AC_tmpPerf = AC_vecBestPerform;
        SC_tmpPerf = SC_vecBestPerform;
        
        for k = 1:K
            topN_Res_Tree = maxk(Tree_tmpPerf(:, k), topN);
            topN_Res_L2 = maxk(L2_tmpPerf(:, k), topN);
            topN_Res_KLD = maxk(KLD_tmpPerf(:, k), topN);
            topN_Res_QGC = maxk(QGC_tmpPerf(:, k), topN);
            topN_Res_QGC2 = maxk(QGC_tmpPerf2(:, k), topN);
            topN_Res_QGC3 = maxk(QGC_tmpPerf3(:, k), topN);
            topN_Res_AC = maxk(AC_tmpPerf(:, k), topN);
            topN_Res_SC = maxk(SC_tmpPerf(:, k), topN);

            for n = 1:topN
                id_Tree = find(Tree_tmpPerf(:, k) == topN_Res_Tree(n), 1);
                id_L2N = find(L2_tmpPerf(:, k) == topN_Res_L2(n), 1);
                id_KLD = find(KLD_tmpPerf(:, k) == topN_Res_KLD(n), 1);
                id_QGC = find(QGC_tmpPerf(:, k) == topN_Res_QGC(n), 1);
                id_QGC2 = find(QGC_tmpPerf2(:, k) == topN_Res_QGC2(n), 1);
                id_QGC3 = find(QGC_tmpPerf3(:, k) == topN_Res_QGC3(n), 1);
                id_AC = find(AC_tmpPerf(:, k) == topN_Res_AC(n), 1);
                id_SC = find(SC_tmpPerf(:, k) == topN_Res_SC(n), 1);
                
                if Tree_tmpPerf(id_Tree, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_Tree, 1);
                    if ~isempty(id)
                        cellRes_Tree(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_Tree{n, k} = 'n / a';
                    end  
                end
                if L2_tmpPerf(id_L2N, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_L2N, 1);
                    if ~isempty(id)
                        cellRes_L2(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_L2{n, k} = 'n / a';
                    end
                end
                if KLD_tmpPerf(id_KLD, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_KLD, 1);
                    if ~isempty(id)
                        cellRes_KLD(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_KLD{n, k} = 'n / a';
                    end
                end
                if QGC_tmpPerf(id_QGC, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_QGC, 1);
                    if ~isempty(id)
                        cellRes_QGC(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_QGC{n, k} = 'n / a';
                    end
                end
                if QGC_tmpPerf2(id_QGC2, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_QGC2, 1);
                    if ~isempty(id)
                        cellRes_QGC2(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_QGC2{n, k} = 'n / a';
                    end
                end
                if QGC_tmpPerf3(id_QGC3, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_QGC3, 1);
                    if ~isempty(id)
                        cellRes_QGC3(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_QGC3{n, k} = 'n / a';
                    end
                end
                if AC_tmpPerf(id_AC, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_AC, 1);
                    if ~isempty(id)
                        cellRes_AC(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_AC{n, k} = 'n / a';
                    end
                end
                if SC_tmpPerf(id_SC, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_SC, 1);
                    if ~isempty(id)
                        cellRes_SC(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_SC{n, k} = 'n / a';
                    end
                end
                
                Tree_tmpPerf(id_Tree, k) = 0;
                L2_tmpPerf(id_L2N, k) = 0;
                KLD_tmpPerf(id_KLD, k) = 0;
                QGC_tmpPerf(id_QGC, k) = 0;
                QGC_tmpPerf2(id_QGC2, k) = 0;
                QGC_tmpPerf3(id_QGC3, k) = 0;
                AC_tmpPerf(id_AC, k) = 0;
                SC_tmpPerf(id_SC, k) = 0;
            end
        end
        if ~isempty(cellVertexNames)
            cellResult{i, 1} = cellVertexNames{1,2}(find(cellVertexNames{1,1} == query));
        else
            cellResult{i, 1} = query;
        end
        
        cellResult{i, 2} = QGC_NCut_BestRel;
        cellResult{i, 3} = entropy_QGC1;
        cellResult{i, 4} = QGC_NCut_BestRel2;
        cellResult{i, 5} = entropy_QGC2;
        cellResult{i, 6} = QGC_NCut_BestRel3;
        cellResult{i, 7} = entropy_QGC3;
        cellResult{i, 8} = AC_NCut_BestRel;
        cellResult{i, 9} = entropy_AC;
        cellResult{i, 10} = SC_NCut_BestRel;
        cellResult{i, 11} = entropy_SC;
        cellResult{i, 12} = cellRes_QGC;
        cellResult{i, 13} = cellRes_QGC2;
        cellResult{i, 14} = cellRes_QGC3;
        cellResult{i, 15} = cellRes_AC;
        cellResult{i, 16} = cellRes_SC;
        
        cellR = cell(10, 1);
        for n = 1:20
            [val , indx] = max(vecRel);
            cellR{n, 1} = cellVertexNames{1,2}(find(cellVertexNames{1,1} == indx, 1));
            vecRel(indx) = 0;
        end
        
    elseif TEST_TYPE > 1 && act_QOGC_QGC_topAll == 1
        
        cellResult{i, 1} = cellVertexNames{1,2}(find(cellVertexNames{1,1} == query));
        
        for m = 1:length(topM)
            % for each m in M, record the performance (QNCut & entropy)
            
            cellResult{i, 6 * (m-1) + 2} = QGC_NCut_BestRel(m);
            cellResult{i, 6 * (m-1) + 3} = entropy_QGC1(m);
            cellResult{i, 6 * (m-1) + 4} = QGC_NCut_BestRel2(m);
            cellResult{i, 6 * (m-1) + 5} = entropy_QGC2(m);
            cellResult{i, 6 * (m-1) + 6} = QGC_NCut_BestRel3(m);
            cellResult{i, 6 * (m-1) + 7} = entropy_QGC3(m);

            %{
            cellR = cell(10, 1);
            for n = 1:20
                [val , indx] = max(vecRel);
                cellR{n, 1} = cellVertexNames{1,2}(find(cellVertexNames{1,1} == indx, 1));
                vecRel(indx) = 0;
            end
            %}
        end
        
        % Record the item for the largest topM 
        tN = max(topM);
        cellRes_QGC = cell(tN, K);
        cellRes_QGC2 = cell(tN, K);
        cellRes_QGC3 = cell(tN, K);

        QGC_tmpPerf = QGC_vecBestPerform;
        QGC_tmpPerf2 = QGC_vecBestPerform2;
        QGC_tmpPerf3 = QGC_vecBestPerform3;

        for k = 1:K
            topN_Res_QGC = maxk(QGC_tmpPerf(:, k), tN);
            topN_Res_QGC2 = maxk(QGC_tmpPerf2(:, k), tN);
            topN_Res_QGC3 = maxk(QGC_tmpPerf3(:, k), tN);

            for n = 1:tN
                id_QGC = find(QGC_tmpPerf(:, k) == topN_Res_QGC(n), 1);
                id_QGC2 = find(QGC_tmpPerf2(:, k) == topN_Res_QGC2(n), 1);
                id_QGC3 = find(QGC_tmpPerf3(:, k) == topN_Res_QGC3(n), 1);

                if QGC_tmpPerf(id_QGC, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_QGC, 1);
                    if ~isempty(id)
                        cellRes_QGC(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_QGC{n, k} = 'n / a';
                    end
                end
                if QGC_tmpPerf2(id_QGC2, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_QGC2, 1);
                    if ~isempty(id)
                        cellRes_QGC2(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_QGC2{n, k} = 'n / a';
                    end
                end
                if QGC_tmpPerf3(id_QGC3, k) ~= 0
                    id = find(cellVertexNames{1,1} == id_QGC3, 1);
                    if ~isempty(id)
                        cellRes_QGC3(n, k) = cellVertexNames{1,2}(id);
                    else
                        cellRes_QGC3{n, k} = 'n / a';
                    end
                end

                QGC_tmpPerf(id_QGC, k) = 0;
                QGC_tmpPerf2(id_QGC2, k) = 0;
                QGC_tmpPerf3(id_QGC3, k) = 0;
            end
        end
        
        cellResult{i, 6 * m + 2} = cellRes_QGC;
        cellResult{i, 6 * m + 3} = cellRes_QGC2;
        cellResult{i, 6 * m + 4} = cellRes_QGC3;
    else 
    end
    
    fprintf('\n query: %d => ', query);
    fprintf('QGC: %f , ', QGC_NCut_BestRel);
    fprintf('QGC2: %f , ', QGC_NCut_BestRel2);
    fprintf('QGC3: %f , ', QGC_NCut_BestRel3);
    fprintf('AC: %f ,', AC_NCut_BestRel);
    fprintf('SC: %f \n', SC_NCut_BestRel);

    Tree_perf = full(Tree_vecBestPerform);
    L2_perf = full(L2_vecBestPerform);
    KLD_perf = full(KLD_vecBestPerform);
    QGC_perf = full(QGC_vecBestPerform);
    AC_perf = full(AC_vecBestPerform);
    
end