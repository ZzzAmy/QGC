function [ y ] = myGG( x, rel )

    %
    % (mathbf{I - rr}^\top) \matthbf{L} (mathbf{I - rr}^\top)
    %
    
    global matG_SC

    tmp = x - rel * (rel' * x);
    tmp = matG_SC * tmp;
    y = tmp - rel * (rel' * tmp);

end

