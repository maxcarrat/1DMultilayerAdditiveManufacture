function [basisFunction, LM] = locationMap(N)

    % location matrix
    % it maps the shape functions local to each element to a global unknown
    % index
    
    LM = zeros(N, 2);
    for i=1:N
        LM(i, 1) = i;
        LM(i, 2) = i+1;
    end

    basisFunction = @basisFunctions;
end