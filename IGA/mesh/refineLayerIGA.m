function [ refinedMesh ] = refineLayerIGA( mesh, refinementDepth, layer, numberOfLayers, p )
%REFINELAYER refine the elements layer of the initial mesh up to a given depth
%   mesh = original coarse mesh
%   refinementDepth = depth of the refinement
%   layer = actual layer
%   numberOfLayers = number of layers

%coarsening of the previous active configuration
mesh = linspace(mesh(1), mesh(end), size(mesh, 2));

%generate refined mesh for the new active configuration
[activeMesh, numberOfElementsPerLayer] = getLayerActiveCoords(mesh, layer, numberOfLayers);
refinedActiveElement = linspace(activeMesh(end-numberOfElementsPerLayer), activeMesh(end), numberOfElementsPerLayer*2^(refinementDepth)+ 1);
refinedMesh = [mesh(1:((layer-1) * numberOfElementsPerLayer)) refinedActiveElement];
  
        
end

