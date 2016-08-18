function [ y ] = myGG2( x, rel, eta )

    %
    % (mathbf{I - rr}^\top) \matthbf{L} (mathbf{I - rr}^\top)
    %
    
    global matG_SC

    y = matG_SC * x -  eta * rel * (rel' * x);
end

