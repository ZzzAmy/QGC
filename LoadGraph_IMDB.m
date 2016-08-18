function [ N, E ] = LoadGraph_IMDB( vetex_filepath, country_filepath, upperbound, lowerBound)
    % vetex_filepath = 'ComedyActors_(2000-2013).txt'
    % country_filepath = 'ComedyCountries_(2000-2013).txt'
    
    global matG
    global cellVertexNames

    %% Read file
    fileVertices = fopen(vetex_filepath);
    cellActorMovie = textscan(fileVertices, '%s%s', 'delimiter', '\t');
    fclose(fileVertices);

    fileCountry = fopen(country_filepath);
    cellMovieCountry = textscan(fileCountry, '%s%s', 'delimiter', '\t');
    fclose(fileCountry);

    %% Use hashmap to index name <--> id
    mapMovie2Country = java.util.HashMap;
    mapMovie2ID = java.util.HashMap;
    mapActor2ID = java.util.HashMap;
    mapID2Movie = java.util.HashMap;
    mapID2Actor = java.util.HashMap;
    mapCountry2ID = java.util.HashMap;
    mapID2Country = java.util.HashMap;

    CountryID = 1;
    s = length(cellMovieCountry{1,1});
    for i = 1:s
        movieName = cellMovieCountry{1,1}{i,1};
        countryName = cellMovieCountry{1,2}{i,1};
        mapMovie2Country.put(movieName, countryName);

        if (mapCountry2ID.containsKey(countryName) == false)
            mapCountry2ID.put(countryName, CountryID);
            mapID2Country.put(CountryID, countryName);
            CountryID = CountryID + 1;
        end
    end

    ActorID = 1;
    MovieID = 1;
    s = length(cellActorMovie{1,1});
    for i = 1:s
        actorName = cellActorMovie{1,1}{i,1};
        movieName = cellActorMovie{1,2}{i,1};

        if (mapActor2ID.containsKey(actorName) == false)
            mapID2Actor.put(ActorID, actorName);
            mapActor2ID.put(actorName, ActorID);
            ActorID = ActorID + 1;
        end

        if (mapMovie2ID.containsKey(actorName) == false)
            mapID2Movie.put(MovieID, movieName);
            mapMovie2ID.put(movieName, MovieID);
            MovieID = MovieID + 1;
        end   
    end

    valActorNum = ActorID - 1;
    valMovieNum = MovieID - 1;
    valCountryNum = CountryID - 1;

    matrixActorMovie = sparse(valActorNum, valMovieNum);
    x_index = zeros(s,1);
    y_index = zeros(s,1);

    for i = 1:s
        actorName = cellActorMovie{1,1}{i,1};
        movieName = cellActorMovie{1,2}{i,1};

        x_index(i,1) = mapActor2ID.get(actorName);
        y_index(i,1) = mapMovie2ID.get(movieName);
    end
    matrixActorMovie = sparse(x_index, y_index, ones(s,1), valActorNum, valMovieNum);
    matG = sparse(matrixActorMovie * matrixActorMovie');
    matG = matG - diag(diag(matG));
    matG = matG .* matG;
    
    [mT_i, mT_j, mT_s] = find(matG);
    
    matG = 1*sparse(mT_i, mT_j, sigmf(mT_s, [1 2]));

    N = length(matG);
    E = nnz(matG);

    cellVertexNames = cell(1,2);
    cellVertexNames{1,1} = int32(zeros(valActorNum,1));
    cellVertexNames{1,2} = cell(valActorNum, 1);
    
    for i = 1:valActorNum
        cellVertexNames{1,1}(i) = i;
        cellVertexNames{1,2}{i} = mapID2Actor.get(i);
    end
    
end
