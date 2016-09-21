function [ layerSolution ] = getIGASolution( solution, layer, problem,...
    numberOfRefinedElementsToBeKept )
%GETLAYERSOLUTION extract the solution coefficients of the actual layer
%   solution = global solution
%   layer = actual layer
%   numberOfLayers = number of layers
%   coords = coarse mesh coordinates

offset = layer;
switch numberOfRefinedElementsToBeKept
    case(1)
        if layer == 1 || problem.p == 1
           layerSolution = solution(offset:end);
        else
           layerSolution = solution(offset+(problem.p-1):end);
        end
    otherwise
        disp('Case not implemented yet!!!');
end
end

