function [ activeCoords, numberOfElementsPerLayer ] = getLayerActiveCoords( coords, layer, numberOfLayers )
%[activeCoords] = GETLAYERACTIVECOORDINATES: get the active coords at a given
%layer of the process
%   coords = 1D mesh
%   layer = actual layer
%   numberOfLayers = number of layers

numberOfElementsPerLayer = floor(size(coords,2)/numberOfLayers);

offset = ceil((size(coords,2) * layer)/numberOfLayers);
activeCoords = coords(1, 1:offset);

end

