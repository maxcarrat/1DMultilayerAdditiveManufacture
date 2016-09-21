function [ layerSolution ] = getSolution( solution, layer, numberOfLayers,...
    coords, refinementDepth, numberOfRefinedElementsToBeKept )
%GETLAYERSOLUTION extract the solution coefficients of the actual layer
%   solution = global solution
%   layer = actual layer
%   numberOfLayers = number of layers
%   coords = coarse mesh coordinates

offset = floor((size(coords,2) * (layer-1))/numberOfLayers + 1);
switch numberOfRefinedElementsToBeKept
    case(1)
        layerSolution = solution(offset:end);
    case(5)
        if layer == 1
            layerSolution = solution(offset:end);
        elseif layer == 2
            layerSolution = solution(offset + 2^refinementDepth - 1:end);
        elseif layer == 3
            layerSolution = solution(offset + 2 * 2^refinementDepth - 2:end);
        elseif layer == 4
            layerSolution = solution(offset + 3 * 2^refinementDepth - 3:end);
        else
            layerSolution = solution(offset + 4 * 2^refinementDepth - 4:end);
        end
    otherwise
        disp('Case not implemented yet!!!');
end
end

