function [ layerSolution ] = getLayerSolution( solution, layer, numberOfLayers, coords )
%GETLAYERSOLUTION extract the solution coefficients of the actual layer
%   solution = global solution
%   layer = actual layer
%   numberOfLayers = number of layers
%   coords = coarse mesh coordinates

offset = floor((size(coords,2) * (layer-1))/numberOfLayers + 1);
layerSolution = solution(offset:end);

end

