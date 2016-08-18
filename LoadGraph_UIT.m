function [ N, E ] = LoadGraph_UIT( vetex_filepath, UIT_filepath, graphtype, upperbound, lowerBound)

    % --- for tag recommendation ---
    % tag.dat: id \t tag
    % vetex_filepath = 'tags.dat'
    % UIT_filepath = 'user_taggedartists.dat'
    % upperbound = 3
    % lowerBound = 0
    %
    % --- for artist recommendation ---
    % tag.dat: id \t artist \t lastfm_web_url \t pic_url
    % vetex_filepath = 'artists2.txt'
    % UIT_filepath = 'user_taggedartists.dat'
    % upperbound = 1
    % lowerBound = 5
    
    global matG
    global cellVertexNames

    %% Read file
    
    fileVertices = fopen(vetex_filepath);
    cellVertexNames = textscan(fileVertices, '%d%s', 'delimiter', '\t');
    fclose(fileVertices);
    fileUTA = fopen(UIT_filepath);
    cellUserTagArtists = textscan(fileUTA, '%d%d%d%d%d%d', 'delimiter', '\t');
    fclose(fileUTA);

    %% Construct matrices
    
    matUTI = (cell2mat(cellUserTagArtists))';
    matTagID = cell2mat(cellVertexNames(1,1));
    
    % Size of TagID list
    tagN = double(max(matUTI(3,:)));

    % Size of ArtistID list
    itemN = double(max(matUTI(2,:)));
    
    matIT = sparse(itemN, tagN);

    for index = matUTI
        valRowIndex = index(2,1);
        valcolIndex = index(3,1);

        if find(matTagID == valcolIndex) > 0
            matIT(valRowIndex, valcolIndex) = matIT(valRowIndex, valcolIndex) + 1;
        end
    end

    if upperbound == -1
        if graphtype == 'i'
            matG = BuildKernel(matIT, lowerBound);
        elseif graphtype == 't'
            matG = BuildKernel(matIT', lowerBound);
        else
        end 
    else
        if graphtype == 'i'
            % Prune unpopular tags

            tagFreq = sum(matIT);

            tooSmall = tagFreq < lowerBound;
            tagFreq(tooSmall) = 0;
            matIT(:, tooSmall) = 0;

            % Prune too popular tags
            threshold = mean(tagFreq) + upperbound * std(tagFreq);
            for indx = 1:length(tagFreq)
                if tagFreq(indx) > threshold
                    matIT(:, indx) = sparse(itemN, 1);
                end
            end
        elseif graphtype == 't'
            % Prune unpopular items
            itemFreq = sum(matIT');

            tooSmall = itemFreq < lowerBound;
            itemFreq(tooSmall) = 0;
            matIT(tooSmall, :) = 0;

            % Prune too popular items
            threshold = mean(itemFreq) + upperbound * std(itemFreq);     
            for indx = 1:length(itemFreq)
                if itemFreq(indx) > threshold
                    matIT(indx, :) = sparse(1, tagN);
                end
            end
        else
        end


        % reweight the influence of a tag / item based on TFIDF.
        matLargeThan1 = (matIT > 0);

        if graphtype == 't'
            vecItem_pop = log2(sum(matLargeThan1, 2) + 2);
            matIT = bsxfun(@rdivide, matIT, vecItem_pop);
            matrixVertex = sparse(matIT' * matIT);
        elseif graphtype == 'i'
            vecTag_pop = log2(sum(matLargeThan1, 1) + 2);
            matIT = bsxfun(@rdivide, matIT, vecTag_pop);
            matrixVertex = sparse(matIT * matIT');
        else
        end

        matrixTrans = matrixVertex - diag(diag(matrixVertex), 0);

        matrixTrans = matrixTrans .* matrixTrans;

        [mT_i, mT_j, mT_s] = find(matrixTrans);

        if graphtype == 't'
            matG = 2*sparse(mT_i, mT_j, sigmf(mT_s, [50 0])); %for tag recommendation
        elseif graphtype == 'i'
            matG = 2*sparse(mT_i, mT_j, sigmf(mT_s, [80 0])); %for artist recommendation
        else
        end
    
    end
    
    N = length(matG);
    E = nnz(matG);
end

