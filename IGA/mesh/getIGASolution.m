function [ layerSolution ] = getIGASolution( solution, layer, numberOfLayers,...
    numberOfControlPoints, refinementDepth, numberOfRefinedElementsToBeKept )
%GETLAYERSOLUTION extract the solution coefficients of the actual layer
%   solution = global solution
%   layer = actual layer
%   numberOfLayers = number of layers
%   coords = coarse mesh coordinates

offset = layer;
switch numberOfRefinedElementsToBeKept
    case(1)
        layerSolution = solution(offset:end);
    otherwise
        disp('Case not implemented yet!!!');
end
end

