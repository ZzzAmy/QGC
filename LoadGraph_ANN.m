function [ N, E ] = LoadGraph_ANN(vetex_filepath)
    % vetex_filepath = 'ComedyActors_(2000-2013).txt'
    
    global matG
    global cellVertexNames

    %% Read file
    fileVertices = fopen(vetex_filepath);
    cellCiteAuthor = textscan(fileVertices, '%s%s', 'delimiter', '\t');
    fclose(fileVertices);

    %% Use hashmap to index name <--> id
    mapAuthor2ID = java.util.HashMap;
    mapID2Author = java.util.HashMap;

    AuthorID = 1;
    s = length(cellCiteAuthor{1,1});
    for i = 1:s
        author1 = cellCiteAuthor{1,1}{i,1};
        author2 = cellCiteAuthor{1,2}{i,1};

        if (mapAuthor2ID.containsKey(author1) == false)
            mapID2Author.put(AuthorID, author1);
            mapAuthor2ID.put(author1, AuthorID);
            AuthorID = AuthorID + 1;
        end

        if (mapAuthor2ID.containsKey(author2) == false)
            mapID2Author.put(AuthorID, author2);
            mapAuthor2ID.put(author2, AuthorID);
            AuthorID = AuthorID + 1;
        end 
    end

    valAuthorNum = AuthorID - 1;

    matG = sparse(valAuthorNum, valAuthorNum);
    x_index = zeros(s,1);
    y_index = zeros(s,1);

    for i = 1:s
        author1 = cellCiteAuthor{1,1}{i,1};
        author2 = cellCiteAuthor{1,2}{i,1};

        x_index(i,1) = mapAuthor2ID.get(author1);
        y_index(i,1) = mapAuthor2ID.get(author2);
    end
    matG = sparse(x_index, y_index, ones(s,1), valAuthorNum, valAuthorNum);
    matG = matG - diag(diag(matG));
    matG = 0.5 * (matG + matG');
    %matG = matG .* matG;
    
    N = length(matG);
    E = nnz(matG);

    cellVertexNames = cell(1,2);
    cellVertexNames{1,1} = int32(zeros(valAuthorNum,1));
    cellVertexNames{1,2} = cell(valAuthorNum, 1);
    
    for i = 1:valAuthorNum
        cellVertexNames{1,1}(i) = i;
        cellVertexNames{1,2}{i} = mapID2Author.get(i);
    end
    
end
