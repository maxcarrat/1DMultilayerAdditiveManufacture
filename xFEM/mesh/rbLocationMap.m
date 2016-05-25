function [basisFunction, LM] = rbLocationMap(modes, numberOfEnrichedElements)

    % RBLOCATIONMAP it maps the shape functions local to each 
    % element to a global unknown index
    
    LM = zeros(numberOfEnrichedElements*modes, 2);

    for i=1:numberOfEnrichedElements*modes
        LM(i, 1) = i;
        LM(i, 2) = i+1;
    end

    basisFunction = @reducedBasis;
    
end