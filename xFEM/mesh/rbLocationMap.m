function [basisFunction, LM] = rbLocationMap(modes)

    % RBLOCATIONMAP it maps the shape functions local to each 
    % element to a global unknown index
    
    LM = zeros(2*modes, 2*modes);

    for i=1:2*modes
        for j=1:2*modes
            LM(i, j) = j;
        end
    end

    basisFunction = @reducedBasis;
    
end